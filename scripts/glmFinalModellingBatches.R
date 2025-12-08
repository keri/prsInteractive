#!/usr/bin/env Rscript

# Set CRAN mirror first
#options(repos = c(CRAN = "https://cran.rstudio.com/"))
#
#library(glmnet)
#
## Conditional installation
#if (!require(data.table, quietly = TRUE)) {
# install.packages("data.table")
# library(data.table)
#}
#
#library(dplyr)
#library(doMC)
#library(pROC)
#library(DescTools)
#library(stringr)
#
#if (!require(argparse, quietly = TRUE)) {
# install.packages("argparse")
# library(argparse)
#}

library(glmnet)
library(data.table)
library(dplyr)
library(doMC)
library(pROC)
library(DescTools)
library(stringr)
library(argparse)

registerDoMC(cores = 60)

# Simple timeout function using setTimeLimit
run_glmnet_with_timeout <- function(X, y, timeout_hours = 1, nfolds = 3, 
                                   penalty.factor = NULL, parallel = TRUE) {
                                      
  timeout_seconds <- timeout_hours * 3600
                                      
  tryCatch({
    cat("Starting cv.glmnet with", timeout_hours, "hour timeout...\n")
    
    # Set the timeout
    setTimeLimit(cpu = timeout_seconds, elapsed = timeout_seconds, transient = TRUE)
    
    # Record start time
    start_time <- Sys.time()
    
    # Run cv.glmnet
    result <- cv.glmnet(X, y, 
                       nfolds = nfolds, 
                       family = "binomial", 
                       penalty.factor = penalty.factor, 
                       parallel = parallel,
                       thresh = 1e-4,  # Less strict convergence
                       maxit = 100000) # More iterations allowed
    
    # Reset timeout
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    
    end_time <- Sys.time()
    duration <- round(difftime(end_time, start_time, units = "mins"), 2)
    cat("cv.glmnet completed in", duration, "minutes\n")
    
    return(list(model = result, success = TRUE, timed_out = FALSE))
    
  }, error = function(e) {
    # Reset timeout in case of error
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
      cat("TIMEOUT: cv.glmnet exceeded", timeout_hours, "hours\n")
      return(list(model = NULL, success = FALSE, timed_out = TRUE))
    } else {
      cat("ERROR in cv.glmnet:", e$message, "\n")
      return(list(model = NULL, success = FALSE, timed_out = FALSE, error = e$message))
    }
  })
}

get_main_cardio_features <- function(feature_file_path, results_path) {

  
  cat("=== PROCESSING MAIN CARDIO FEATURES ===\n")
  cat("Feature file:", feature_file_path, "\n")
  cat("Results path:", results_path, "\n")
  
  # Step 1: Read the feature file
  if (!file.exists(feature_file_path)) {
    stop(paste("Feature file not found:", feature_file_path))
  }
  
  feature_data <- fread(feature_file_path)
  cat("Loaded feature data with", nrow(feature_data), "rows and", ncol(feature_data), "columns\n")
  cat("Columns:", paste(names(feature_data), collapse = ", "), "\n")
  
  # Step 2: Filter for main_E features (main_E == 1)
  main_E_features <- feature_data[main_E == 1]
  cat("Found", nrow(main_E_features), "main_E features\n")
  
  if (nrow(main_E_features) == 0) {
    stop("No main_E features found in the data")
  }
  
  # Step 3: Extract unique environmental feature names
  main_cardio_feature_list <- unique(main_E_features$envFeature)
  main_cardio_feature_list <- main_cardio_feature_list[!is.na(main_cardio_feature_list)]
  
  cat("Main cardio features to extract:\n")
  for(i in 1:min(10, length(main_cardio_feature_list))) {
    cat(sprintf("  %d. %s\n", i, main_cardio_feature_list[i]))
  }
  if(length(main_cardio_feature_list) > 10) {
    cat(sprintf("  ... and %d more\n", length(main_cardio_feature_list) - 10))
  }
  
  # Step 4: Read participant environment data
  participant_env_file <- file.path(results_path, "participant_environment.csv")
  
  if (!file.exists(participant_env_file)) {
    stop(paste("Participant environment file not found:", participant_env_file))
  }
  
  cat("\nReading participant environment data from:", participant_env_file, "\n")
  participant_env <- fread(participant_env_file)
  
  # Handle different possible column names for participant ID
  if ("Participant ID" %in% names(participant_env)) {
    setnames(participant_env, "Participant ID", "IID")
  } else if (!"IID" %in% names(participant_env)) {
    # Check for other common ID column names
    id_cols <- c("ID", "id", "ParticipantID", "participant_id")
    found_id_col <- intersect(id_cols, names(participant_env))
    if (length(found_id_col) > 0) {
      setnames(participant_env, found_id_col[1], "IID")
      cat("Renamed", found_id_col[1], "to IID\n")
    } else {
      stop("No IID or participant ID column found in environment data")
    }
  }
  
  # 

  # Ensure IID is integer type for consistent merging
  participant_env[, IID := as.integer(IID)]
  
  cat("Participant environment data shape:", nrow(participant_env), "x", ncol(participant_env), "\n")
  cat("Available columns:", paste(head(names(participant_env), 15), collapse = ", "), 
      if(ncol(participant_env) > 15) "..." else "", "\n")
  
  # Step 5: Check which main cardio features are available in participant data
  available_features <- intersect(main_cardio_feature_list, names(participant_env))
  missing_features <- setdiff(main_cardio_feature_list, names(participant_env))
  
  cat("\nFeature availability check:\n")
  cat("Available features:", length(available_features), "/", length(main_cardio_feature_list), "\n")
  
  if (length(missing_features) > 0) {
    cat("Missing features:\n")
    for(i in 1:min(10, length(missing_features))) {
      cat(sprintf("  %s\n", missing_features[i]))
    }
    if(length(missing_features) > 10) {
      cat(sprintf("  ... and %d more\n", length(missing_features) - 10))
    }
  }
  
  if (length(available_features) == 0) {
    stop("No main cardio features found in participant environment data")
  }
  
  # Remove features with too much missing data (>5%)
  high_quality_features <- c()
  for(feature in available_features) {
    missing_rate <- sum(is.na(participant_env[[feature]])) / nrow(participant_env)
    
    if(missing_rate <= 0.05) {
      high_quality_features <- c(high_quality_features, feature)
    } else {
      cat(sprintf("Excluding %s: Participant Data %.1f%%, (>5%% threshold)\n", 
                  feature, missing_rate*100))
    }
  }
  
  cat("\nFeatures passing 5% missing data threshold:", length(high_quality_features), "/", length(available_features), "\n")
  
  
  # Update datasets to include only high-quality features
  final_columns <- c("IID", high_quality_features)
  cardio_env_data <- participant_env[, ..final_columns]
  
  cat("\nExtracted cardio environment data with", ncol(cardio_env_data), "columns\n")
  
  return(cardio_env_data)
}



calculate_nagelkerke_rsquared <- function(yTest,yProba) {
  
  # Fit the null model (intercept only)
  null_model <- glm(yTest ~ 1, family = binomial)
  
  # Extract log-likelihood of the null model
  logLik_null <- logLik(null_model)
  
  
  # Manually calculate the log-likelihood of the fitted model
  logLik_fit <- sum(yTest * log(yProba) + (1 - yTest) * log(1 - yProba))
  
  
  # Calculate Cox-Snell's R^2
  R2_cox_snell <- 1 - exp((2 * (logLik_null - logLik_fit)) / length(yTest))
  print('cox snell R squared =')
  print(R2_cox_snell)
  
  # Calculate Nagelkerke's R^2
  R2_nagelkerke <- R2_cox_snell / (1 - exp(2 * logLik_null / length(yTest)))
  print('Nagelkerke r squared = ')
  print(R2_nagelkerke)
  
  return(R2_nagelkerke)
  
}

# Fixed version of download_covariate_data to ensure consistent data types

download_covariate_data <- function(covar_file) {
  # Load required package
  if (!require(data.table, quietly = TRUE)) {
    stop("data.table package is required")
  }
  
  # Check if file exists
  if (!file.exists(covar_file)) {
    stop(paste("File not found:", covar_file))
  }
  
  cat('Getting covariate data .....\n')
  
  # Determine file type and separator
  file_ext <- tolower(tools::file_ext(covar_file))
  separator <- if (file_ext == "csv") "," else " "
  
  cat("Reading file as", if (file_ext == "csv") "CSV" else "space-separated", "format...\n")
  
  # Try to read the file and see what columns are available
  tryCatch({
    # First, peek at the file to see available columns
    df_peek <- fread(covar_file, sep = separator, nrows = 1)
    available_cols <- names(df_peek)
    cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
    
    # Define preferred column names (in order of preference)
    column_preferences <- list(
      IID = c("IID", "id", "ID"),
      age = c("age", "Age", "AGE"),
      SEX = c("SEX", "Sex", "sex", "Gender"),
      PC1 = "PC1", PC2 = "PC2", PC3 = "PC3", PC4 = "PC4", PC5 = "PC5",
      PC6 = "PC6", PC7 = "PC7", PC8 = "PC8", PC9 = "PC9", PC10 = "PC10"
    )
    
    # Find which columns exist
    columns_to_read <- c()
    column_mapping <- list()
    
    for (standard_name in names(column_preferences)) {
      possible_names <- column_preferences[[standard_name]]
      found_col <- intersect(possible_names, available_cols)[1]  # Take first match
      
      if (!is.na(found_col)) {
        columns_to_read <- c(columns_to_read, found_col)
        column_mapping[[found_col]] <- standard_name
        cat("Using '", found_col, "' for ", standard_name, "\n", sep = "")
      } else {
        cat("Column for", standard_name, "not found\n")
      }
    }
    
    if (length(columns_to_read) == 0) {
      stop("No matching columns found in the file")
    }
    
    # Read the file with found columns
    df <- fread(covar_file, sep = separator, select = columns_to_read)
    
    # Rename columns to standard names
    for (old_name in names(column_mapping)) {
      if (old_name %in% names(df)) {
        setnames(df, old_name, column_mapping[[old_name]])
      }
    }
    
    # Handle missing IID - create if needed
    if (!"IID" %in% names(df)) {
      cat("IID column not found, creating from row numbers...\n")
      df[, IID := paste0("ID_", 1:.N)]
    }
    
    # CRITICAL: Ensure IID is integer type for consistent merging
    df[, IID := as.integer(IID)]
    cat("IID column set to integer type\n")
    
    # Handle Sex column conversion
    if ("Sex" %in% names(df)) {
      cat("Processing Sex column...\n")
      if (is.character(df$Sex) || is.factor(df$Sex)) {
        df[, SEX := ifelse(tolower(as.character(Sex)) %in% c("female", "f","Female",'F'), 0, 1)]
      }
      df[, SEX := as.numeric(SEX)]
    }
    
    # Convert age to numeric if present
    if ("age" %in% names(df)) {
      df[, age := as.numeric(age)]
    }
    
    # Convert PC columns to numeric
    pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
    if (length(pc_cols) > 0) {
      for (col in pc_cols) {
        df[, (col) := as.numeric(get(col))]
      }
    }
    
    # Remove any rows with all missing values
    df <- df[rowSums(is.na(df)) < ncol(df)]
    
    # Set IID as key
    if ("IID" %in% names(df)) {
      setkey(df, IID)
    }
    
    cat("Successfully loaded", nrow(df), "rows and", ncol(df), "columns\n")
    cat("Final columns:", paste(names(df), collapse = ", "), "\n")
    cat("IID data type:", class(df$IID), "\n")
    
    return(df)
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    cat("Attempting fallback approach...\n")
    
    # Fallback: try reading all columns and manually select
    tryCatch({
      df_all <- fread(covar_file, sep = separator)
      cat("Fallback successful. All columns:", paste(names(df_all), collapse = ", "), "\n")
      
      # Return the data.table as is, but ensure we have an ID column
      if (!"IID" %in% names(df_all) && !"ID" %in% names(df_all)) {
        df_all[, IID := paste0("ID_", 1:.N)]
      } else if ("ID" %in% names(df_all) && !"IID" %in% names(df_all)) {
        setnames(df_all, "ID", "IID")
      }
      
      # Ensure IID is integer
      if ("IID" %in% names(df_all)) {
        df_all[, IID := as.integer(IID)]
        setkey(df_all, IID)
      }
      
      return(df_all)
      
    }, error = function(e2) {
      stop(paste("Both main and fallback approaches failed:", e2$message))
    })
  })
}


create_epi_df <- function(epiDf, pairList,combo="sum") {
  # Convert epiDf to data.table if it is not already
  epiDf <- as.data.table(epiDf)
  
  # Initialize an empty data.table for the final result
  epiArrayFinal <- data.table()
  
  # Loop through each pair in pairList
  for (pair in pairList) {
    
    snps <- strsplit(pair, ",")[[1]]
    

    
    if (combo == "sum") {
      snps_combined <- rowSums(epiDf[, ..snps, drop = FALSE])
    } else {
      snps_combined <- apply(epiDf[, ..snps, drop = FALSE], 1, prod)
    }
    # Create a data.table for the multiplied SNPs with column name as pair
#   temp_dt <- data.table(temp = snps_product)
    temp_dt <- data.table(temp = snps_combined)
    setnames(temp_dt, "temp", pair)
    
    # Concatenate with the final result
    epiArrayFinal <- cbind(epiArrayFinal, temp_dt)
  }
  
  # Bind the IID column to array
  epiArrayFinal <- cbind(epiDf$IID, epiArrayFinal)
  
  # Rename the column 'V1' to 'IID'
  setnames(epiArrayFinal, "V1", "IID")
  
  return(epiArrayFinal)
}


#get_geno_dataset <- function(data_file,data_path,columns_to_get){
# 
# start_time <- Sys.time()
# 
# # Adding 'IID' and 'PHENOTYPE' to the list of columns to get
# columns_to_get <- c('IID', 'PHENOTYPE', columns_to_get)
# 
# # Load the main dataset
# main_data <- fread(data_file, sep=" ", select = columns_to_get)  
# 
# # Load the file containing the values to match (assuming no header)
# withdrawal_path = paste0(data_path,"/withdrawals.csv")
# withdrawn <- fread(withdrawal_path, header = FALSE)
# 
# # Rename the column in values_to_remove to 'IID' (or whatever the matching column is in main_data)
# setnames(withdrawn, "V1", "IID")
# 
# # Remove withdrawn individuals from the DataFrame
# main_data <- main_data[!main_data$IID %in% withdrawn$IID]
# 
# # Modify the 'PHENOTYPE' column in place
# main_data[PHENOTYPE == 1, PHENOTYPE := 0]
# main_data[PHENOTYPE == 2, PHENOTYPE := 1]
# 
# end_time <- Sys.time()
# print(paste('time it took to download entire dataset is', round(difftime(end_time, start_time, units = "mins"), 2), 'minutes'))
# 
# return(main_data)
# 
#}

get_geno_dataset <- function(data_file, data_path, columns_to_get){
    
    start_time <- Sys.time()
    
    # Adding 'IID' and 'PHENOTYPE' to the list of columns to get
    columns_to_get <- c('IID', 'PHENOTYPE', columns_to_get)
    
    # Check how many columns we're actually reading
    print(paste("Reading", length(columns_to_get), "columns from file"))
    
    # Load the main dataset with progress tracking
    main_data <- fread(data_file, 
                        sep = " ", 
                        select = columns_to_get,
                        showProgress = TRUE,
                        verbose = TRUE) 
    
    print(paste("Dataset loaded. Dimensions:", nrow(main_data), "x", ncol(main_data)))
    gc()  # Force garbage collection after loading
    
    # Load withdrawals
    withdrawal_path = paste0(data_path, "/withdrawals.csv")
    withdrawn <- fread(withdrawal_path, header = FALSE)
    setnames(withdrawn, "V1", "IID")
    
    print(paste("Removing", nrow(withdrawn), "withdrawn individuals"))
    
    # Remove withdrawn individuals (more memory efficient)
    main_data <- main_data[!IID %in% withdrawn$IID]
    rm(withdrawn)  # Remove withdrawn data from memory
    gc()
    
    print(paste("After withdrawals. Dimensions:", nrow(main_data), "x", ncol(main_data)))
    
    # Modify the 'PHENOTYPE' column in place
    main_data[PHENOTYPE == 1, PHENOTYPE := 0]
    main_data[PHENOTYPE == 2, PHENOTYPE := 1]
    
    end_time <- Sys.time()
    print(paste('Time taken:', round(difftime(end_time, start_time, units = "mins"), 2), 'minutes'))
    
    return(main_data)
  }


get_geno_read_table <- function(data_file, data_path, columns_to_get, test_mode = FALSE, test_rows = 1000) {
  
  start_time <- Sys.time()
  print("=== USING read.table (NO MEMORY MAPPING) ===")
  
  if (test_mode) {
    print(paste("TEST MODE: Reading only", test_rows, "rows"))
  } else {
    print("PRODUCTION MODE: Reading entire file")
  }
  
  # Load withdrawals
  withdrawal_path <- paste0(data_path, "/withdrawals.csv")
  withdrawn <- fread(withdrawal_path, header = FALSE)
  withdrawn_set <- as.integer(withdrawn$V1)
  print(paste("Loaded", length(withdrawn_set), "withdrawal IDs"))
  rm(withdrawn); gc()
  
  # Get header
  con <- file(data_file, "r")
  header_line <- readLines(con, n = 1)
  close(con)
  
  header_cols <- trimws(strsplit(header_line, " ")[[1]])
  header_cols <- header_cols[header_cols != ""]
  
  print(paste("File has", length(header_cols), "columns"))
  
  # Find column positions
  base_cols <- c('IID', 'PHENOTYPE')
  col_positions <- match(base_cols, header_cols)
  
  if (any(is.na(col_positions))) {
    stop("IID or PHENOTYPE columns not found")
  }
  
  print(paste("IID at position", col_positions[1], ", PHENOTYPE at position", col_positions[2]))
  
  # Use read.table with colClasses to read only specific columns
  print(paste("Reading first", test_rows, "rows using read.table..."))
  
  # Create colClasses vector - NULL for unwanted columns, appropriate type for wanted
  col_classes <- rep("NULL", length(header_cols))  # Skip all columns by default
  col_classes[col_positions[1]] <- "integer"       # IID as integer
  col_classes[col_positions[2]] <- "integer"       # PHENOTYPE as integer
  
  # Add requested columns
  if (length(columns_to_get) > 0) {
    requested_positions <- match(columns_to_get, header_cols)
    if (any(is.na(requested_positions))) {
      missing <- columns_to_get[is.na(requested_positions)]
      stop("Requested columns not found: ", paste(missing, collapse = ", "))
    }
    
    # In production mode, add ALL requested columns
    if (test_mode) {
      # Test mode: just first few columns
      test_snp_cols <- columns_to_get[1:min(5, length(columns_to_get))]
      test_positions <- match(test_snp_cols, header_cols)
      col_classes[test_positions] <- "numeric"
      print(paste("TEST MODE: Will read", length(test_snp_cols), "SNP columns"))
    } else {
      # Production mode: ALL requested columns
      col_classes[requested_positions] <- "numeric"
      print(paste("PRODUCTION MODE: Will read", length(columns_to_get), "SNP columns"))
    }
  }
  
  # Read using read.table (avoids memory mapping)
  main_data <- read.table(data_file, 
                         header = TRUE,
                         sep = " ",
                         colClasses = col_classes,
                         nrows = if (test_mode) test_rows else -1,  # Full file in production
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote = "",
                         check.names = FALSE)  # CRITICAL: Preserve original column names
  
  print(paste("Successfully loaded:", nrow(main_data), "x", ncol(main_data)))
  
  # Convert to data.table
  main_data <- as.data.table(main_data)
  
  # Filter withdrawals
  original_rows <- nrow(main_data)
  main_data <- main_data[!IID %in% withdrawn_set]
  print(paste("After withdrawal filtering:", nrow(main_data), "rows kept"))
  
  # Modify phenotype
  main_data[PHENOTYPE == 1, PHENOTYPE := 0L]
  main_data[PHENOTYPE == 2, PHENOTYPE := 1L]
  
  print("Sample data:")
  print(names(main_data)[1:10])
  print(paste("IID type:", class(main_data$IID), "- good for merging!"))
  
  end_time <- Sys.time()
  print(paste('Time taken:', round(difftime(end_time, start_time, units = "secs"), 2), 'seconds'))
  
  return(main_data)
}

# PRODUCTION WRAPPERS - Full file processing
get_geno_production_read_table <- function(data_file, data_path, columns_to_get) {
  print("=== PRODUCTION RUN WITH read.table ===")
  print(paste("Processing", length(columns_to_get), "columns from entire file"))
  
  result <- get_geno_read_table(
    data_file = data_file,
    data_path = data_path,
    columns_to_get = columns_to_get,
    test_mode = FALSE  # Full file
  )
  
  return(result)
}

# UPDATED PRODUCTION WRAPPERS with column name fixing
get_geno_production_fixed <- function(data_file, data_path, columns_to_get) {
  print("=== PRODUCTION RUN WITH COLUMN NAME FIXING ===")
  print(paste("Processing", length(columns_to_get), "columns from entire file"))
  
  result <- get_geno_read_table_fixed(
    data_file = data_file,
    data_path = data_path,
    columns_to_get = columns_to_get,
    test_mode = FALSE  # Full file
  )
  
  return(result)
}

test_read_table_fixed <- function(data_file, data_path, columns_to_get) {
  print("=== TEST RUN WITH COLUMN NAME FIXING ===")
  
  result <- get_geno_read_table_fixed(
    data_file = data_file,
    data_path = data_path,
    columns_to_get = columns_to_get,
    test_mode = TRUE,
    test_rows = 1000
  )
  
  return(result)
}

# Helper function to check column name changes
check_column_names <- function(original_names, loaded_data_names) {
  print("=== COLUMN NAME VALIDATION ===")
  
  # Check for mismatches
  mismatches <- list()
  
  for (orig_name in original_names) {
    if (!orig_name %in% loaded_data_names) {
      # Look for similar names (R might have converted characters)
      possible_matches <- loaded_data_names[grepl(gsub("[^A-Za-z0-9]", ".", orig_name), loaded_data_names)]
      
      if (length(possible_matches) > 0) {
        mismatches[[orig_name]] <- possible_matches[1]
        print(paste("MISMATCH:", orig_name, "->", possible_matches[1]))
      } else {
        print(paste("MISSING:", orig_name, "not found in loaded data"))
      }
    }
  }
  
  if (length(mismatches) > 0) {
    print(paste("Found", length(mismatches), "column name changes"))
    return(mismatches)
  } else {
    print("✓ All column names preserved correctly")
    return(NULL)
  }
}

# Fixed version that handles column name changes
get_geno_read_table_fixed <- function(data_file, data_path, columns_to_get, test_mode = FALSE, test_rows = 1000) {
  
  start_time <- Sys.time()
  print("=== USING read.table WITH COLUMN NAME FIXING ===")
  
  if (test_mode) {
    print(paste("TEST MODE: Reading only", test_rows, "rows"))
  } else {
    print("PRODUCTION MODE: Reading entire file")
  }
  
  # Load withdrawals
  withdrawal_path <- paste0(data_path, "/withdrawals.csv")
  withdrawn <- fread(withdrawal_path, header = FALSE)
  withdrawn_set <- as.integer(withdrawn$V1)
  print(paste("Loaded", length(withdrawn_set), "withdrawal IDs"))
  rm(withdrawn); gc()
  
  # Get original header
  con <- file(data_file, "r")
  header_line <- readLines(con, n = 1)
  close(con)
  
  original_header_cols <- trimws(strsplit(header_line, " ")[[1]])
  original_header_cols <- original_header_cols[original_header_cols != ""]
  
  print(paste("File has", length(original_header_cols), "columns"))
  
  # Find column positions using original names
  base_cols <- c('IID', 'PHENOTYPE')
  col_positions <- match(base_cols, original_header_cols)
  
  if (any(is.na(col_positions))) {
    stop("IID or PHENOTYPE columns not found")
  }
  
  # Create colClasses vector
  col_classes <- rep("NULL", length(original_header_cols))
  col_classes[col_positions[1]] <- "integer"  # IID
  col_classes[col_positions[2]] <- "integer"  # PHENOTYPE
  
  # Add requested columns
  if (length(columns_to_get) > 0) {
    requested_positions <- match(columns_to_get, original_header_cols)
    if (any(is.na(requested_positions))) {
      missing <- columns_to_get[is.na(requested_positions)]
      stop("Requested columns not found: ", paste(missing, collapse = ", "))
    }
    
    if (test_mode) {
      test_snp_cols <- columns_to_get[1:min(5, length(columns_to_get))]
      test_positions <- match(test_snp_cols, original_header_cols)
      col_classes[test_positions] <- "numeric"
      print(paste("TEST MODE: Will read", length(test_snp_cols), "SNP columns"))
    } else {
      col_classes[requested_positions] <- "numeric"
      print(paste("PRODUCTION MODE: Will read", length(columns_to_get), "SNP columns"))
    }
  }
  
  # Read using read.table with column name preservation
  main_data <- read.table(data_file, 
                         header = TRUE,
                         sep = " ",
                         colClasses = col_classes,
                         nrows = if (test_mode) test_rows else -1,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote = "",
                         check.names = FALSE)  # Preserve original names
  
  print(paste("Successfully loaded:", nrow(main_data), "x", ncol(main_data)))
  
  # Check if column names were preserved
  loaded_col_names <- names(main_data)
  expected_cols <- if (test_mode) {
    c(base_cols, columns_to_get[1:min(5, length(columns_to_get))])
  } else {
    c(base_cols, columns_to_get)
  }
  
  # Validate column names
  name_changes <- check_column_names(expected_cols, loaded_col_names)
  
  # Fix any column name changes
  if (!is.null(name_changes)) {
    for (orig_name in names(name_changes)) {
      changed_name <- name_changes[[orig_name]]
      if (changed_name %in% loaded_col_names) {
        setnames(main_data, changed_name, orig_name)
        print(paste("FIXED:", changed_name, "->", orig_name))
      }
    }
  }
  
  # Convert to data.table
  main_data <- as.data.table(main_data)
  
  # Final column name check
  print("Final column names:")
  print(names(main_data)[1:10])
  
  # Filter withdrawals
  original_rows <- nrow(main_data)
  main_data <- main_data[!IID %in% withdrawn_set]
  print(paste("After withdrawal filtering:", nrow(main_data), "rows kept"))
  
  # Modify phenotype
  main_data[PHENOTYPE == 1, PHENOTYPE := 0L]
  main_data[PHENOTYPE == 2, PHENOTYPE := 1L]
  
  print("Sample data:")
  print(head(main_data[, 1:3], 3))
  print(paste("IID type:", class(main_data$IID), "- good for merging!"))
  
  end_time <- Sys.time()
  print(paste('Time taken:', round(difftime(end_time, start_time, units = "mins"), 2), 'minutes'))
  
  return(main_data)
}


# TEST WRAPPERS - Limited rows for validation
test_read_table <- function(data_file, data_path, columns_to_get) {
  print("=== TEST RUN WITH read.table ===")
  
  result <- get_geno_read_table(
    data_file = data_file,
    data_path = data_path,
    columns_to_get = columns_to_get,
    test_mode = TRUE,
    test_rows = 1000
  )
  
  return(result)
}

test_scan_method <- function(data_file, data_path, columns_to_get) {
  print("=== TEST RUN WITH scan() ===")
  
  result <- get_geno_scan_method(
    data_file = data_file,
    data_path = data_path,
    columns_to_get = columns_to_get,
    test_mode = TRUE,
    test_rows = 1000
  )
  
  return(result)
}

get_epi_snps <- function(epiFeatures) {
  # Initialize an empty vector to store SNPs
  epiSnps <- c()
  
  # Loop through each pair in epiFeatures
  for (pair in epiFeatures) {
    # Split the pair into individual SNPs and append to epiSnps
    snps <- unlist(strsplit(pair, ","))
    epiSnps <- c(epiSnps, snps)
  }
  
  # Get unique SNPs
  epiSnps <- unique(epiSnps)
  
  # Return the list of unique SNPs
  return(epiSnps)
}

#get_important_features <- function(feature_file){
#
# #file has 3 columns : [feature,data_type:(main,epi)]
## feature_file = paste0(feature_pathway,'/importantFeaturesForAssociationAnalysis.csv')
# #feature_file = paste0(feature_pathway,'/importantFeaturesPostShap.csv')
# important_features = read.csv(feature_file)
#
# epi_main_features = c(important_features$feature)
#
# main_features = important_features$feature[important_features$data_type == "main"]
#
# epi_features = important_features$feature[important_features$data_type == "epi"]
#
# return(list(epi_main_features=epi_main_features,main_features=main_features,epi_features=epi_features))
#}

get_important_features <- function(feature_file){
  
  print("=== LOADING AND CLEANING FEATURE LIST ===")
  
  # Load feature file
  important_features = read.csv(feature_file)
  
  print(paste("Original features loaded:", nrow(important_features)))
  
  # Check for duplicates before processing
  duplicate_features <- important_features$feature[duplicated(important_features$feature)]
  if (length(duplicate_features) > 0) {
    print(paste("Found", length(duplicate_features), "duplicate features"))
    cat("First few duplicates:", paste(head(duplicate_features, 10), collapse = ", "), "\n")
    
    # Remove duplicates - keep first occurrence
    important_features <- important_features[!duplicated(important_features$feature), ]
    print(paste("After removing duplicates:", nrow(important_features), "unique features"))
  } else {
    print("No duplicates found")
  }
  
  # Create the feature lists with unique values
  epi_main_features = unique(important_features$feature)  # All unique features
  
  main_features = unique(important_features$feature[important_features$data_type == "main"])
  
  epi_features = unique(important_features$feature[important_features$data_type == "epi"])
  
  # Print summary
  cat("Feature summary:\n")
  cat("- Total unique features:", length(epi_main_features), "\n")
  cat("- Main features:", length(main_features), "\n")
  cat("- Epi features:", length(epi_features), "\n")
  
  # Verify no overlaps if you expect them to be mutually exclusive
  overlap <- intersect(main_features, epi_features)
  if (length(overlap) > 0) {
    print(paste("WARNING: Found", length(overlap), "features in both main and epi"))
    cat("First few overlapping features:", paste(head(overlap, 5), collapse = ", "), "\n")
  }
  
  return(list(
    epi_main_features = epi_main_features,
    main_features = main_features,
    epi_features = epi_features
  ))
}

# Also add this function to check data types before merging:
check_merge_compatibility <- function(df1, df2, merge_col = "IID") {
  cat("\n=== MERGE COMPATIBILITY CHECK ===\n")
  cat("Dataset 1 (", deparse(substitute(df1)), "):\n")
  cat("  Rows:", nrow(df1), "\n")
  cat("  ", merge_col, "type:", class(df1[[merge_col]]), "\n")
  cat("  Sample", merge_col, "values:", paste(head(df1[[merge_col]], 3), collapse = ", "), "\n")
  
  cat("Dataset 2 (", deparse(substitute(df2)), "):\n")
  cat("  Rows:", nrow(df2), "\n")
  cat("  ", merge_col, "type:", class(df2[[merge_col]]), "\n")
  cat("  Sample", merge_col, "values:", paste(head(df2[[merge_col]], 3), collapse = ", "), "\n")
  
  # Check for common values
  common_ids <- length(intersect(df1[[merge_col]], df2[[merge_col]]))
  cat("Common", merge_col, "values:", common_ids, "\n")
  
  if (common_ids == 0) {
    cat("WARNING: No common", merge_col, "values found!\n")
  }
  
  # Check data types
  if (class(df1[[merge_col]]) != class(df2[[merge_col]])) {
    cat("WARNING: Data type mismatch for", merge_col, "column!\n")
    return(FALSE)
  }
  
  cat("Merge compatibility: OK\n\n")
  return(TRUE)
}


################################ GLOBAL VARIABLES  ########################
# 
parser <- ArgumentParser()
parser$add_argument("--results_path", required = TRUE)
parser$add_argument("--data_path", required = TRUE)
parser$add_argument("--hla_file", required = TRUE) 
parser$add_argument("--covar_file", required = TRUE)
parser$add_argument("--pheno_data", required = TRUE)
parser$add_argument("--training_file", required = TRUE)
parser$add_argument("--test_file", required = TRUE)
parser$add_argument("--training_env_gen_file", required = TRUE)
parser$add_argument("--test_env_gen_file", required = TRUE)
parser$add_argument("--feature_model_file" ,required = TRUE)
# 
args <- parser$parse_args()

results_path <- args$results_path
data_path <- args$data_path
pheno_path <- args$pheno_data
hla_file <- args$hla_file
training_file <- args$training_file
test_file <- args$test_file
training_env_gen_file <- args$training_env_gen_file
test_env_gen_file <- args$test_env_gen_file
covar_file <- args$covar_file
feature_model_file <- args$feature_model_file


#results_path <- '/Users/kerimulterer/prsInteractive/results'
#pheno_path <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
#data_path <- '/Users/kerimulterer/prsInteractive/data'
#hla_file <- '/Users/kerimulterer/prsInteractive/results/participant_hla.csv'
## training_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw'
#training_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw'
#test_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/testCombined.raw'
## test_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/testCombined_final.raw'
#training_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTraining.csv'
#test_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTest.csv'
#covar_file='/Users/kerimulterer/prsInteractive/results/covar.csv'
#feature_model_file='/Users/kerimulterer/prsInteractive/results/type2Diabetes/scores/importantFeaturesPostShap.csv'

scores_path = paste0(pheno_path,'/scores')
#covar_pathway = paste0(results_path,'/covar.txt')


############################### FEATURES TO BE USED IN EACH MODEL ###############

#combined epi and main snps after filtering for FDR
all_features = get_important_features(feature_model_file)

epi_main_features = all_features$epi_main_features
main_features = all_features$main_features
epi_features = all_features$epi_features
columns_to_get = get_epi_snps(epi_main_features)


####################################################################################
#    CREATE FILES FOR MODEL SCORES, FEATURE SCORES, AND PREDICTION PROBABILITIES   #
####################################################################################

###### MODEL SCORE FILE ##########
model_scores_file = paste0(scores_path,'/modelScoresReducedFinalModelBatches.csv')
model_score_colnames = c("auc","ci_lower","ci_upper","nagelkerke_rquared","model","batch")
# Check if the file exists
if (!file.exists(model_scores_file)) {
  
  # If the file doesn't exist, create it with the column headings
  write.table(
    x = data.frame(matrix(ncol = length(model_score_colnames), nrow = 0)),
    file = model_scores_file,
    col.names = model_score_colnames,
    row.names = FALSE,
    sep = ",",
    quote = FALSE
  )
  
  cat("Model score file created successfully.\n")
} else {
  cat("Model score file already exists.\n")
}

###### PREDICTION FILE ##########

predictions_file = paste0(scores_path,'/predictProbsReducedFinalModelBatches.csv')
predictions_colnames = c("yProba","model","batch")

# Check if the file exists
if (!file.exists(predictions_file)) {
  # If the file doesn't exist, create it with the column headings
  write.table(
    x = data.frame(matrix(ncol = length(predictions_colnames), nrow = 0)),
    file = predictions_file,
    col.names = predictions_colnames,
    row.names = FALSE,
    sep = ",",
    quote = FALSE
  )
  
  cat("prediction probability file created successfully.\n")
} else {
  cat("prediction probability file already exists.\n")
}

# ###### FEATURE SCORE FILE ##########
#feature_scores_file = paste0(scores_path,'/importantFeaturesForAssociationAnalysis.csv')
feature_scores_file = paste0(scores_path,'/featureScoresReducedFinalModelBatches.csv')
feature_scores_colnames = c("coefs","model","feature","batch")

# Check if the file exists
if (!file.exists(feature_scores_file)) {
  # If the file doesn't exist, create it with the column headings
  write.table(
    x = data.frame(matrix(ncol = length(feature_scores_colnames), nrow = 0)),
    file = feature_scores_file,
    col.names = feature_scores_colnames,
    row.names = FALSE,
    sep = ",",
    quote = FALSE
  )
  
  cat("feature scores file created successfully.\n")
} else {
  cat("feature scores file already exists.\n")
}




####################################################################################
#                   DOWNLOAD THE TRAINING AND TEST GENO AND COVARIATE DATA         #
####################################################################################

######################### DATASETS TO BE USED IN EACH MODEL #####################

################### COVARIATE DATA ######################
#covariate data with IID + first 10 PCs, SEX, and AGE
covariate_data = download_covariate_data(covar_file)

################### HLA DATA ######################


# hlaData = fread(paste0(machine_path,'/ukbiobank/tanigawaData/HLAImputationCleaned_participant.csv'),sep=",")
hlaData = fread(hla_file,sep=",")
setnames(hlaData, "Participant ID", "IID")

# CRITICAL: Ensure IID is integer type
hlaData[, IID := as.integer(IID)]
cat("HLA IID column set to integer type\n")

hla_columns = setdiff(names(hlaData), "IID")

################## CARDIO METABOLIC DATA ##########

# epi_cardio_participant_feature_list = download_important_cardiometabolic_features(feature_pathway,env_type,machine_path)
envTraining <- fread(training_env_gen_file)
envTest <- fread(test_env_gen_file)

# CRITICAL: Ensure IID is integer type in both
envTraining[, IID := as.integer(IID)]
envTest[, IID := as.integer(IID)]
cat("Environment data IID columns set to integer type\n")

epi_cardio_features = setdiff(names(envTraining), "IID")

############## MAIN CARDIO FEATURES TO MODEL ##########

main_cardio_file=paste0(pheno_path,'/scores/','cardioMetabolicimportantFeaturesPostShap.csv')
main_cardio_df = get_main_cardio_features(main_cardio_file,results_path)
main_cardio_features = setdiff(names(main_cardio_df),'IID')


################### ARRAY TYPE DATA #####################
#get the genotype batch : column Genotype measurement batch
# array_path = paste0(machine_path,'/ukbiobank/tanigawaData/participant.csv')
#array_path = paste0(results_path,'/participant.csv')
#array_data <- fread(array_path, sep=",", select = c('Participant ID','Genotype measurement batch'))
#colnames(array_data) = c("IID","array")
#change batch to 0 and 1 if UK in name for the 2 types of arrays
#array_data$array <- ifelse(grepl("UK", array_data$array), 0, 1)


################# TRAINING data #####################
# trainingDf = get_geno_dataset(paste0(training_path,'/data/',training_file),machine_path,columns_to_get)
trainingDf = get_geno_production_fixed(training_file,data_path,columns_to_get)

#phenotype
yTraining = trainingDf$PHENOTYPE

#combined main SNPs and epi pairs in which haplotypes are added for snp1 and snp2 of pair
trainingDf = create_epi_df(trainingDf,epi_main_features)

#merge covariate data to geno data
# Before merging with covariate data:
# check_merge_compatibility(trainingDf, covariate_data)
trainingDf = merge(trainingDf,covariate_data, by = "IID", all.x = TRUE)

#merge HLA data to geno data
# Before merging with covariate data:
# check_merge_compatibility(trainingDf, hla_data)
trainingDf = merge(trainingDf,hlaData, by = "IID", all.x = TRUE)

#combine cardio and geno features and scale data after combined
#shape will be length(epi geno features to combine list X # people in test)

trainingDf = merge(trainingDf,envTraining, by = "IID", all.x = TRUE)

#merge with cardio main features
trainingDf = merge(trainingDf,main_cardio_df, by= "IID", all.x = TRUE)

################# TEST data #####################
testDf = get_geno_production_fixed(test_file,data_path,columns_to_get)

#phenotype
yTest = testDf$PHENOTYPE

#combined main SNPs and epi pairs in which haplotypes are added for snp1 and snp2 of pair
testDf = create_epi_df(testDf,epi_main_features)


#merge covariate data to geno data
testDf = merge(testDf,covariate_data, by = "IID", all.x = TRUE)

#merge HLA data to geno data
testDf = merge(testDf,hlaData, by = "IID", all.x = TRUE)

testDf = merge(testDf,envTest, by = "IID", all.x = TRUE)

#merge with cardio main features
testDf = merge(testDf,main_cardio_df, by= "IID", all.x = TRUE)


#####################################################################################
#                        PROCESS DATA FOR MODELLING                                 #
#                        1) impute missing data with the mean                       #
#####################################################################################

# Loop through each column except IID and impute missing values with the mean
numeric_cols <- names(trainingDf)[sapply(trainingDf, is.numeric) & names(trainingDf) != "IID"]

cat("Processing", length(numeric_cols), "numeric columns for imputation\n")
cat("Sample columns:", paste(head(numeric_cols, 5), collapse = ", "), "\n")

# Check for missing data before imputation
missing_before_train <- sum(is.na(trainingDf))
missing_before_test <- sum(is.na(testDf))
cat("Missing values before imputation - Training:", missing_before_train, "Test:", missing_before_test, "\n")

# More robust approach - calculate means one by one and apply immediately
cat("Calculating training means and applying imputation...\n")

# Apply imputation column by column
for (col in numeric_cols) {
  # Calculate mean from training data only
  train_mean <- mean(trainingDf[[col]], na.rm = TRUE)
  
  # Check if we have a valid mean
  if (is.na(train_mean) || !is.finite(train_mean)) {
    cat("Warning: Invalid mean for column", col, "- using 0\n")
    train_mean <- 0
  }
  
  # Apply to training data
  missing_count_train <- sum(is.na(trainingDf[[col]]))
  if (missing_count_train > 0) {
    #cat("Imputing", missing_count_train, "missing values in training", col, "with mean", round(train_mean, 4), "\n")
    trainingDf[is.na(get(col)), (col) := train_mean]
  }
  
  # Apply to test data (if column exists)
  if (col %in% names(testDf)) {
    missing_count_test <- sum(is.na(testDf[[col]]))
    if (missing_count_test > 0) {
      #cat("Imputing", missing_count_test, "missing values in test", col, "with mean", round(train_mean, 4), "\n")
      testDf[is.na(get(col)), (col) := train_mean]
    }
  } else {
    cat("Column", col, "not found in test data - skipping\n")
  }
}

# Verify imputation worked
missing_after_train <- sum(is.na(trainingDf))
missing_after_test <- sum(is.na(testDf))

cat("Imputation completed successfully\n")
cat("Missing values after imputation - Training:", missing_after_train, "Test:", missing_after_test, "\n")

# Final check - ensure no missing values remain in key columns
if (missing_after_train > 0) {
  cat("Warning: Some missing values remain in training data\n")
  # Show which columns still have missing values
  missing_cols <- sapply(trainingDf, function(x) sum(is.na(x)))
  missing_cols <- missing_cols[missing_cols > 0]
  if (length(missing_cols) > 0) {
    cat("Columns with remaining missing values:\n")
    print(missing_cols)
  }
}

if (missing_after_test > 0) {
  cat("Warning: Some missing values remain in test data\n")
  # Show which columns still have missing values
  missing_cols <- sapply(testDf, function(x) sum(is.na(x)))
  missing_cols <- missing_cols[missing_cols > 0]
  if (length(missing_cols) > 0) {
    cat("Columns with remaining missing values:\n")
    print(missing_cols)
  }
}


####################################################################################################
#                                LOOP THROUGH ALL DATASETS FOR MODELLING                           #
####################################################################################################

#GET COVARIATE ONLY DATA FOR COVARIATE ONLY MODEL + ARRAY

#get covariate columns to set penalty values to 0
covariate_columns = setdiff(names(covariate_data),'IID')

all_features = setdiff(names(trainingDf),'IID')
#list of lists include the data_type and features in data_type

#get a list of epi cardio features to use in modelling
# epi_cardio_features = setdiff(names(cardioTraining), "IID")
# epi_cardio_features = setdiff(names(envTestDf), "IID")

dataset_list = list(
  list('main',c(main_features,covariate_columns)),
  list('epi',c(epi_features,covariate_columns)),
  list('epi+main',c(epi_main_features,covariate_columns)),
  list('cardio',c(epi_cardio_features,covariate_columns)),
  list('hla',c(covariate_columns,hla_columns)),
  list('cardio_main',c(main_cardio_features,covariate_columns)),
  list('covariate',c(covariate_columns)),
  list('all',c(epi_cardio_features,epi_main_features,hla_columns,covariate_columns)),
  list('all+cardio_main',c(epi_cardio_features,epi_main_features,hla_columns,main_cardio_features,covariate_columns))
  
)

#####################################################################################
#                               ASSIGN PENALTIES                                    #
#####################################################################################

########################. START LOOP THROUGH EACH DATASET. #########################

for (dataset in dataset_list) {
  
  data_type = dataset[[1]]
  cat(paste0('data type being modelled .... ',data_type,'\n'))
  
  all_columns = dataset[[2]]
  
  # Ensure the columns exist in your data and exclude IID
  all_columns = intersect(all_columns, names(testDf))
  all_columns = setdiff(all_columns, 'IID')
  
  # Separate non-covariate features from covariates
  non_covariate_features = setdiff(all_columns, covariate_columns)
  available_covariates = intersect(all_columns, covariate_columns)
  
  cat(paste0('Total features for ', data_type, ': ', length(non_covariate_features), 
             ' + ', length(available_covariates), ' covariates\n'))
  
  # If only covariates, process as before
  if (length(non_covariate_features) == 0 || data_type == "covariate") {
    
    columns_in_model = available_covariates
    
    cat(paste0('Using ', length(columns_in_model), ' covariate columns for ', data_type, ' model\n'))
    
    ########## set the penalties ############
    p.fac <- rep(0, length(columns_in_model))  # All covariates get penalty factor 0
    
    #####################################################################################
    #                        RUN THE GLM LASSO MODEL                                    #
    #####################################################################################
    
    # Prepare the design matrix
    Xtraining <- as.matrix(trainingDf[, ..columns_in_model])
    Xtest <- as.matrix(testDf[, ..columns_in_model])
    
    # Use GLM for covariate-only model
    cvModel <- glm(yTraining ~ ., data = as.data.frame(Xtraining), family = binomial())
    coefs <- as.data.frame(as.matrix(coef(cvModel)))
    
    ############## GET THE COEFS ###############
    colnames(coefs) <- c("coefs")
    coefs$model = data_type
    coefs$feature = rownames(coefs)
    coefs$batch = "_batch_1"
    rownames(coefs) = NULL
    
    # Save coefficients
    if (!file.exists(feature_scores_file)) {
      cat('creating feature scores file ....\n')
      write.table(coefs, file = feature_scores_file, row.names = FALSE, sep = ",")
    } else {
      write.table(coefs, file = feature_scores_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      cat('appending to feature scores file ....\n')
    }
    
    #############################################################################
    #                        MAKE PREDICTIONS                                   #
    #############################################################################
    
    yProba <- predict(cvModel, newdata = as.data.frame(Xtest), type = "response")
    yHat <- ifelse(yProba > 0.5, 1, 0)
    
    yProba_vector = as.numeric(yProba)
    yProba_df <- data.frame(yProba = yProba_vector)
    yProba_df$model = data_type
    yProba_df$batch = "_batch_1"
    yProba_df$IID = testDf$IID
    
    # Save predictions
    if (!file.exists(predictions_file)) {
      cat('creating predictions file .... \n')
      write.table(yProba_df, file = predictions_file, row.names = FALSE, sep = ",")
    } else {
      write.table(yProba_df, file = predictions_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      cat('appending data to predictions file ... \n')
    }
    
    ############################################################################
    #                               CALCULATE MODEL PERFORMANCE                #
    ############################################################################
    
    # Calculate the AUC
    roc_curve <- roc(yTest, yProba_vector)
    auc_value <- auc(roc_curve)
    test_result <- roc_curve$auc.p.value
    ci_intervals = ci.auc(roc_curve)
    
    # Calculate nagelkerke r squared
    n_rsquared = calculate_nagelkerke_rsquared(yTest, yProba_vector)
    
    # Convert CI results to a data frame
    model_scores <- data.frame(
      auc = c(auc_value[1]),
      ci_lower = c(ci_intervals[1]),
      ci_upper = c(ci_intervals[3]),
      nagelkerke_rsquared = c(n_rsquared),
      model =  data_type,
      batch = "_batch_1"
    )
    
    # Save model scores
    if (!file.exists(model_scores_file)) {
      print('creating model scores file ....')
      write.table(model_scores, file = model_scores_file, row.names = FALSE, sep = ",")
      cat('creating model scores file with scored features....\n')
    } else {
      write.table(model_scores, file = model_scores_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      cat('writing to model scores file ......\n')
    }
    
  } else {
    
    # Process features in batches of 1000
    batch_size = 100
    n_batches = ceiling(length(non_covariate_features) / batch_size)
    
    cat(paste0('Processing ', length(non_covariate_features), ' non-covariate features in ', 
               n_batches, ' batches of up to ', batch_size, ' features each\n'))
    
    for (batch_num in 1:n_batches) {
      
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat("Starting batch", batch_num, "of", n_batches, "for", data_type, "\n")
      
      # Calculate batch indices
      start_idx = (batch_num - 1) * batch_size + 1
      end_idx = min(batch_num * batch_size, length(non_covariate_features))
      batch_features = non_covariate_features[start_idx:end_idx]
      columns_in_model = c(batch_features, available_covariates)
      
      cat("Batch features:", length(batch_features), "+ covariates:", length(available_covariates), 
          "= total:", length(columns_in_model), "\n")
      
      # Set penalty factors
      p.fac <- rep(1, length(columns_in_model))
      p.fac[columns_in_model %in% covariate_columns] <- 0
      
      # Prepare matrices
      Xtraining <- as.matrix(trainingDf[, ..columns_in_model])
      Xtest <- as.matrix(testDf[, ..columns_in_model])
      
      # Run with timeout
      batch_start_time <- Sys.time()
      
      glm_result <- run_glmnet_with_timeout(
        X = Xtraining, 
        y = yTraining, 
        timeout_hours = 1,  # Adjust as needed
        nfolds = 3,
        penalty.factor = p.fac,
        parallel = TRUE
      )
      
      # Handle timeout
      if (glm_result$timed_out) {
        cat("TIMEOUT: Skipping batch", batch_num, "after 1 hours\n")
        
        # Log the timeout
        timeout_entry <- data.frame(
          model = paste0(data_type, "_batch_", batch_num),
          status = "TIMEOUT_1_HOUR",
          n_features = length(batch_features),
          timestamp = as.character(Sys.time()),
          stringsAsFactors = FALSE
        )
        
        timeout_log_file <- paste0(scores_path, '/batch_timeouts.csv')
        if (!file.exists(timeout_log_file)) {
          write.csv(timeout_entry, timeout_log_file, row.names = FALSE)
        } else {
          write.table(timeout_entry, timeout_log_file, row.names = FALSE, 
                      col.names = FALSE, sep = ",", append = TRUE)
        }
        
        # Clean up and move to next batch
        rm(Xtraining, Xtest)
        gc()
        next
      }
      
      # Handle other failures
      if (!glm_result$success) {
        cat("FAILED: Batch", batch_num, "failed with error:", 
            ifelse(is.null(glm_result$error), "Unknown error", glm_result$error), "\n")
        rm(Xtraining, Xtest)
        gc()
        next
      }
      
      # Success - continue with normal processing
      cvModel <- glm_result$model
      cat("SUCCESS: Batch", batch_num, "completed\n")
      
      # Get coefficients
      coefs <- as.data.frame(as.matrix(coef(cvModel, s = "lambda.min")))
      colnames(coefs) <- c("coefs")
      coefs$model = data_type
      coefs$batch = batch_num
      coefs$feature = rownames(coefs)
      rownames(coefs) = NULL
      
      # Save coefficients
      if (!file.exists(feature_scores_file)) {
        write.table(coefs, file = feature_scores_file, row.names = FALSE, sep = ",")
      } else {
        write.table(coefs, file = feature_scores_file, row.names = FALSE, 
                    col.names = FALSE, sep = ",", append = TRUE)
      }
      
      # Make predictions
      yProba <- predict(cvModel, newx = Xtest, type = "response", s="lambda.min")
      yProba_vector = as.numeric(yProba)
      yProba_df <- data.frame(yProba = yProba_vector)
      yProba_df$model = data_type
      yProba_df$batch = batch_num
      yProba_df$IID = testDf$IID
      
      # Save predictions
      if (!file.exists(predictions_file)) {
        write.table(yProba_df, file = predictions_file, row.names = FALSE, sep = ",")
      } else {
        write.table(yProba_df, file = predictions_file, row.names = FALSE, 
                    col.names = FALSE, sep = ",", append = TRUE)
      }
      
      # Calculate performance metrics
      roc_curve <- roc(yTest, yProba_vector)
      auc_value <- auc(roc_curve)
      ci_intervals = ci.auc(roc_curve)
      n_rsquared = calculate_nagelkerke_rsquared(yTest, yProba_vector)
      
      model_scores <- data.frame(
        auc = c(auc_value[1]),
        ci_lower = c(ci_intervals[1]),
        ci_upper = c(ci_intervals[3]),
        nagelkerke_rsquared = c(n_rsquared),
        model = data_type,
        batch = batch_num
      )
      
      # Save model scores
      if (!file.exists(model_scores_file)) {
        write.table(model_scores, file = model_scores_file, row.names = FALSE, sep = ",")
      } else {
        write.table(model_scores, file = model_scores_file, row.names = FALSE, 
                    col.names = FALSE, sep = ",", append = TRUE)
      }
      
      # Report batch completion
      batch_end_time <- Sys.time()
      batch_duration <- round(difftime(batch_end_time, batch_start_time, units = "mins"), 2)
      cat("Batch", batch_num, "total time:", batch_duration, "minutes\n")
      
      # Clean up memory
      rm(cvModel, coefs, yProba, yProba_df, model_scores, roc_curve, Xtraining, Xtest)
      gc()
      
    } # End batch loop
    
    
  } # End else (non-covariate processing)
  
} # End dataset loop

  