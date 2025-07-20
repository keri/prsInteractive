# Set CRAN mirror first
options(repos = c(CRAN = "https://cran.rstudio.com/"))

library(glmnet)

# Conditional installation
if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
  library(data.table)
}

library(dplyr)
library(doMC)
library(pROC)
library(DescTools)
library(stringr)

if (!require(argparse, quietly = TRUE)) {
  install.packages("argparse")
  library(argparse)
}

#library(glmnet)
#library(data.table)
#library(dplyr)
#library(doMC)
#library(pROC)
#library(DescTools)
#library(stringr)
#library(argparse)

registerDoMC(cores = 18)



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
    
    # Handle SEX column conversion
    if ("Sex" %in% names(df)) {
      cat("Processing SEX column...\n")
      if (is.character(df$SEX) || is.factor(df$SEX)) {
        df[, SEX := ifelse(tolower(as.character(SEX)) %in% c("fnemale", "f"), 0, 1)]
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
      
      if ("IID" %in% names(df_all)) {
        setkey(df_all, IID)
      }
      
      return(df_all)
      
    }, error = function(e2) {
      stop(paste("Both main and fallback approaches failed:", e2$message))
    })
  })
}


create_epi_df <- function(epiDf, pairList, cardio_df=data.table()) {
  # Convert epiDf to data.table if it is not already
  epiDf <- as.data.table(epiDf)
  
  # Initialize an empty data.table for the final result
  epiArrayFinal <- data.table()
  
  # Loop through each pair in pairList
  for (pair in pairList) {
    
    snps <- strsplit(pair, ",")[[1]]
    # Multiply the columns corresponding to the SNPs (element-wise product)
#   snps_product <- apply(epiDf[, ..snps], 1, prod)
    snps_sum <- rowSums(epiDf[, ..snps])
    
    # Create a data.table for the multiplied SNPs with column name as pair
#   temp_dt <- data.table(temp = snps_product)
    temp_dt <- data.table(temp = snps_sum)
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


# get_geno_dataset <- function(file_path,data_path,columns_to_get){
#   
#   start_time <- Sys.time()
#   
#   # Adding 'IID' and 'PHENOTYPE' to the list of columns to get
#   columns_to_get <- c('IID', 'PHENOTYPE', columns_to_get)
#   # columns_to_get <- c('IID', 'PHENOTYPE')
#   
#   # Load the main dataset
#   main_data <- fread(file_path, sep=" ", select = columns_to_get)  
#   
#   # Load the file containing the values to match (assuming no header)
#   withdrawal_path = paste0(data_path,"/withdrawals.csv")
#   withdrawn <- fread(withdrawal_path, header = FALSE)
#   
#   # Rename the column in values_to_remove to 'IID' (or whatever the matching column is in main_data)
#   setnames(withdrawn, "V1", "IID")
#   
#   # Remove withdrawn individuals from the DataFrame
#   main_data <- main_data[!main_data$IID %in% withdrawn$IID]
#   
#   # Modify the 'PHENOTYPE' column in place
#   main_data[PHENOTYPE == 1, PHENOTYPE := 0]
#   main_data[PHENOTYPE == 2, PHENOTYPE := 1]
#   
#   end_time <- Sys.time()
#   print(paste('time it took to download entire dataset is', round(difftime(end_time, start_time, units = "mins"), 2), 'minutes'))
#   
#   return(main_data)
#   
# }

get_geno_dataset <- function(file_path, data_path, columns_to_get) {
  
  library(data.table)
  start_time <- Sys.time()
  
  cat("Reading file with mixed separators...\n")
  
  # Read the file manually since header and data have different separators
  all_lines <- readLines(file_path)
  
  # Parse header (tab-separated)
  header_line <- all_lines[1]
  col_names <- unlist(strsplit(header_line, "\t"))
  cat("Found", length(col_names), "columns from tab-separated header\n")
  cat("First 10 columns:", paste(head(col_names, 10), collapse = ", "), "\n")
  
  # Parse data lines (space-separated)
  data_lines <- all_lines[-1]
  cat("Processing", length(data_lines), "data rows...\n")
  
  # Parse first data line to check field count
  first_data <- unlist(strsplit(data_lines[1], "\\s+"))
  cat("First data line has", length(first_data), "fields\n")
  cat("First 10 data values:", paste(head(first_data, 10), collapse = ", "), "\n")
  
  # Adjust column names if needed
  if (length(first_data) != length(col_names)) {
    cat("Adjusting column count from", length(col_names), "to", length(first_data), "\n")
    if (length(first_data) < length(col_names)) {
      col_names <- col_names[1:length(first_data)]
    } else {
      # Add extra column names if data has more fields
      extra_cols <- paste0("extra_", 1:(length(first_data) - length(col_names)))
      col_names <- c(col_names, extra_cols)
    }
  }
  
  cat("Using", length(col_names), "column names\n")
  
  # Create a temporary file with consistent separators
  temp_file <- tempfile()
  cat("Creating temporary file with consistent formatting...\n")
  
  # Write header with tabs
  writeLines(paste(col_names, collapse = "\t"), temp_file)
  
  # Process and write data lines
  cat("Converting data lines to tab-separated format...\n")
  
  # Process in chunks for better memory management
  chunk_size <- 100
  for (i in seq(1, length(data_lines), chunk_size)) {
    end_idx <- min(i + chunk_size - 1, length(data_lines))
    chunk_lines <- data_lines[i:end_idx]
    
    # Convert each line from space-separated to tab-separated
    processed_lines <- character(length(chunk_lines))
    for (j in 1:length(chunk_lines)) {
      fields <- unlist(strsplit(chunk_lines[j], "\\s+"))
      # Pad or truncate to match column count
      if (length(fields) < length(col_names)) {
        fields <- c(fields, rep(NA, length(col_names) - length(fields)))
      } else if (length(fields) > length(col_names)) {
        fields <- fields[1:length(col_names)]
      }
      processed_lines[j] <- paste(fields, collapse = "\t")
    }
    
    # Append to temp file
    write.table(processed_lines, temp_file, append = TRUE, quote = FALSE, 
                row.names = FALSE, col.names = FALSE, sep = "\t")
    
    if (i %% 500 == 1) cat("Processed", min(end_idx, length(data_lines)), "of", length(data_lines), "rows\n")
  }
  
  cat("Reading processed file...\n")
  
  # Now read the properly formatted file
  main_data <- fread(temp_file, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  
  # Clean up temp file
  unlink(temp_file)
  
  cat("Successfully loaded data:", nrow(main_data), "x", ncol(main_data), "\n")
  
  # Check column names
  cat("Column names (first 15):", paste(head(names(main_data), 15), collapse = ", "), "\n")
  
  # Find SNP columns
  found_snps <- intersect(columns_to_get, names(main_data))
  missing_snps <- setdiff(columns_to_get, names(main_data))
  
  cat("Found", length(found_snps), "of", length(columns_to_get), "requested SNP columns\n")
  if (length(missing_snps) > 0 && length(missing_snps) <= 10) {
    cat("Missing SNPs:", paste(missing_snps, collapse = ", "), "\n")
  }
  
  # Select columns
  essential_cols <- c("IID", "PHENOTYPE")
  final_columns <- c(essential_cols, found_snps)
  
  # Ensure all selected columns exist
  final_columns <- intersect(final_columns, names(main_data))
  main_data <- main_data[, final_columns, with = FALSE]
  
  cat("Selected", ncol(main_data), "columns\n")
  
  # Convert data types
  cat("Converting data types...\n")
  if ("PHENOTYPE" %in% names(main_data)) {
    main_data[, PHENOTYPE := as.numeric(PHENOTYPE)]
  }
  
  # Convert SNP columns to numeric
  snp_cols <- setdiff(names(main_data), c("IID", "PHENOTYPE"))
  for (col in snp_cols) {
    main_data[, (col) := as.numeric(get(col))]
  }
  
  # Handle withdrawals
  if ("IID" %in% names(main_data) && nrow(main_data) > 0) {
    withdrawal_path <- paste0(data_path, "/withdrawals.csv")
    if (file.exists(withdrawal_path)) {
      cat("Processing withdrawals...\n")
      withdrawn <- fread(withdrawal_path, header = FALSE)
      setnames(withdrawn, "V1", "IID")
      
      initial_count <- nrow(main_data)
      main_data <- main_data[!IID %in% withdrawn$IID]
      removed_count <- initial_count - nrow(main_data)
      cat("Removed", removed_count, "withdrawn participants\n")
    }
  }
  
  # Handle phenotype recoding
  if ("PHENOTYPE" %in% names(main_data) && nrow(main_data) > 0) {
    cat("Recoding phenotypes (1->0, 2->1)...\n")
    pheno_before <- table(main_data$PHENOTYPE, useNA = "ifany")
    cat("Before recoding:\n")
    print(pheno_before)
    
    main_data[PHENOTYPE == 1, PHENOTYPE := 0]
    main_data[PHENOTYPE == 2, PHENOTYPE := 1]
    
    pheno_after <- table(main_data$PHENOTYPE, useNA = "ifany")
    cat("After recoding (0=control, 1=case):\n")
    print(pheno_after)
  }
  
  end_time <- Sys.time()
  time_taken <- round(difftime(end_time, start_time, units = "mins"), 2)
  cat("Processing completed in", time_taken, "minutes\n")
  
  # Final summary
  cat("\n=== FINAL SUMMARY ===\n")
  cat("Dimensions:", nrow(main_data), "x", ncol(main_data), "\n")
  if (nrow(main_data) > 0) {
    cat("Sample data (first 3 rows, first 5 columns):\n")
    print(main_data[1:min(3, nrow(main_data)), 1:min(5, ncol(main_data))])
  }
  
  return(main_data)
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

get_important_features <- function(feature_pathway){

  #file has 3 columns : [feature,data_type:(main,epi)]
# feature_file = paste0(feature_pathway,'/importantFeaturesForAssociationAnalysis.csv')
  feature_file = paste0(feature_pathway,'/importantFeaturesPostShap.csv')
  important_features = read.csv(feature_file)

  epi_main_features = c(important_features$feature)

  main_features = important_features$feature[important_features$data_type == "main"]

  epi_features = important_features$feature[important_features$data_type == "epi"]

  return(list(epi_main_features=epi_main_features,main_features=main_features,epi_features=epi_features))
}


################################ GLOBAL VARIABLES  ########################

#parser <- ArgumentParser()
#parser$add_argument("--results_path", required = TRUE)
#parser$add_argument("--data_path", required = TRUE)
#parser$add_argument("--hla_file", required = TRUE) 
#parser$add_argument("--covar_file", required = TRUE)
#parser$add_argument("--pheno_path", required = TRUE)
#parser$add_argument("--training_file", required = TRUE)
#parser$add_argument("--test_file", required = TRUE)
#parser$add_argument("--training_env_gen_file", required = TRUE)
#parser$add_argument("--test_env_gen_file", required = TRUE)
#
#args <- parser$parse_args()

#results_path <- args$results_path
#data_path <- args$data_path
#pheno_path <- args$pheno_path
#hla_file <- args$hla_file
#training_file <- args$training_file
#test_file <- args$test_file
#training_env_gen_file <- args$training_env_gen_file
#test_env_gen_file <- args$test_env_gen_file
#covar_file <- args$covar_file


results_path <- '/Users/kerimulterer/prsInteractive/results'
data_path <- '/Users/kerimulterer/prsInteractive/data'
pheno_path <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
hla_file <- '/Users/kerimulterer/prsInteractive/results/participant_hla.csv'
training_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw'
test_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/testCombined.raw'
training_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTraining.csv'
test_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTest.csv'
covar_file='/Users/kerimulterer/prsInteractive/results/covar.csv'

scores_path = paste0(pheno_path,'/scores')
#covar_pathway = paste0(results_path,'/covar.txt')


############################### FEATURES TO BE USED IN EACH MODEL ###############

#combined epi and main snps after filtering for FDR
all_features = get_important_features(scores_path)

epi_main_features = all_features$epi_main_features
main_features = all_features$main_features
epi_features = all_features$epi_features
columns_to_get = get_epi_snps(epi_main_features)


####################################################################################
#    CREATE FILES FOR MODEL SCORES, FEATURE SCORES, AND PREDICTION PROBABILITIES   #
####################################################################################

###### MODEL SCORE FILE ##########
model_scores_file = paste0(scores_path,'/modelScoresReducedFinalModel.csv')
model_score_colnames = c("auc","ci_lower","ci_upper","nagelkerke_rquared","model")
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

predictions_file = paste0(scores_path,'/predictProbsReducedFinalModel.csv')
predictions_colnames = c("yProba","model")

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

feature_scores_file = paste0(scores_path,'/featureScoresReducedFinalModel.csv')
feature_scores_colnames = c("coefs","model","feature")

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
hla_columns = setdiff(names(hlaData), "IID")

################## CARDIO METABOLIC DATA ##########

# epi_cardio_participant_feature_list = download_important_cardiometabolic_features(feature_pathway,env_type,machine_path)
envTraining <- fread(training_env_gen_file)
envTest <- fread(test_env_gen_file)
epi_cardio_features = setdiff(names(envTraining), "IID")

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
trainingDf = get_geno_dataset(training_file,data_path,columns_to_get)

#phenotype
yTraining = trainingDf$PHENOTYPE

#combined main SNPs and epi pairs in which haplotypes are added for snp1 and snp2 of pair
trainingDf = create_epi_df(trainingDf,epi_main_features)

#merge covariate data to geno data
trainingDf = merge(trainingDf,covariate_data, by = "IID", all.x = TRUE)

#merge HLA data to geno data
trainingDf = merge(trainingDf,hlaData, by = "IID", all.x = TRUE)

#combine cardio and geno features and scale data after combined
#shape will be length(epi geno features to combine list X # people in test)

trainingDf = merge(trainingDf,envTraining, by = "IID", all.x = TRUE)

################# TEST data #####################
testDf = get_geno_dataset(test_file,data_path,columns_to_get)

#phenotype
yTest = testDf$PHENOTYPE

#combined main SNPs and epi pairs in which haplotypes are added for snp1 and snp2 of pair
testDf = create_epi_df(testDf,epi_main_features)


#merge covariate data to geno data
testDf = merge(testDf,covariate_data, by = "IID", all.x = TRUE)

#merge HLA data to geno data
testDf = merge(testDf,hlaData, by = "IID", all.x = TRUE)

testDf = merge(testDf,envTest, by = "IID", all.x = TRUE)

#####################################################################################
#                        PROCESS DATA FOR MODELLING                                 #
#                        1) impute missing data with the mean                       #
#####################################################################################

# Loop through each column except IID and impute missing values with the mean
numeric_cols <- names(trainingDf)[sapply(trainingDf, is.numeric) & names(trainingDf) != "IID"]

################ TRAINING. ################

trainingDf[, (numeric_cols) := lapply(.SD, function(x)  {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}), .SDcols = numeric_cols]



############# TEST  ######################

testDf[, (numeric_cols) := lapply(.SD, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}), .SDcols = numeric_cols]


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
  list('main',c(main_features,covariate_columns,hla_columns)),
  list('epi',c(epi_features,covariate_columns,hla_columns)),
  list('epi+main',c(epi_main_features,covariate_columns,hla_columns)),
  list('cardio',c(epi_cardio_features,covariate_columns,hla_columns)),
  list('all',c(all_features,covariate_columns)),
  # list('cardio_main',c(main_cardio_features,covariate_columns)),
  list('covariate',c(covariate_columns))
  
)

#####################################################################################
#                               ASSIGN PENALTIES                                    #
#####################################################################################



########################. START LOOP THROUGH EACH DATASET. #########################

for (dataset in dataset_list) {

  data_type = dataset[[1]]
  cat(paste0('data type being modelled .... ',data_type,'\n'))
  
  
  columns_in_model = dataset[[2]]
  
  # Instead, ensure the columns exist in your data and exclude IID
  columns_in_model = intersect(columns_in_model, names(testDf))
  columns_in_model = setdiff(columns_in_model, 'IID')
  
  cat(paste0('Using ', length(columns_in_model), ' columns for ', data_type, ' model\n'))
  
  ########## set the penalties ############
  
  #Set penalty factor for all columns to 1
  p.fac <- rep(1, length(columns_in_model))
  
  # Update values to 0 for columns in the covariate data
  p.fac[columns_in_model %in% covariate_columns] <- 0
  
  
  
  
  #####################################################################################
  #                        RUN THE GLM LASSO MODEL                                    #
  #####################################################################################
  
  # Prepare the design matrix (excluding the 'IID' column)
  Xtraining <- as.matrix(trainingDf[, ..columns_in_model])
  Xtest <- as.matrix(testDf[,..columns_in_model])
  
  if (data_type != "covariate") {
    
    cvModel <- cv.glmnet(Xtraining, yTraining, nfolds=5,family = "binomial", penalty.factor = p.fac, parallel = TRUE)
    coefs <- as.data.frame(as.matrix(coef(cvModel, s = "lambda.min")))
    
  } else {
    
    cvModel <- glm(yTraining ~ ., data = as.data.frame(Xtraining), family = binomial())
    coefs <- as.data.frame(as.matrix(coef(cvModel)))
  }
  
  
  ############## GET THE COEFS ###############
  #capture the data into a table to add onto
  
  
  colnames(coefs) <- c("coefs")
  coefs$model = rep(data_type,length(coefs))
  coefs$feature = rownames(coefs)
  rownames(coefs) = NULL
  
  # Check if the file exists
  if (!file.exists(feature_scores_file)) {
    
    cat('creating feature scores file ....\n')
    
    # If the file doesn't exist, create it with column names
    write.table(coefs, file = feature_scores_file, row.names = FALSE, sep = ",")
  } else {
    # If the file exists, append the new data without column names
    write.table(coefs, file = feature_scores_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
    
    cat('appending to feature scores file ....\n')
  }
  
  #############################################################################
  #                        MAKE PREDICTIONS                                   #
  #############################################################################
  
  if (data_type != "covariate") {
    
    yHat = predict(cvModel,Xtest,type='class',s="lambda.min")
    
    # Predict probabilities for the training data
    yProba <- predict(cvModel, newx = Xtest, type = "response",s="lambda.min")
    
    
  } else {
    
    yProba <- predict(cvModel, newdata = as.data.frame(Xtest), type = "response")
    yHat <- ifelse(yProba > 0.5, 1, 0)
    
  }
  
  yProba_vector = as.numeric(yProba)
  yProba <- data.frame(yProba_vector)
  colnames(yProba) = 'yProba'
  yProba$model = rep(data_type,length(yProba))
  
  #Check if the file exists
  if (!file.exists(predictions_file)) {
    
    cat('creating predictions file .... \n')
    
    # If the file doesn't exist, create it with column names
    write.table(yProba, file = predictions_file, row.names = FALSE, sep = ",")
  } else {
    # If the file exists, append the new data without column names
    write.table(yProba, file = predictions_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
    
    cat('appending data to predictions file ... \n')
    
  }
  
  
  ############################################################################
  #                               CALCULATE MODEL PERFORMANCE                #
  ############################################################################
  
  model_scores = data.table()
  # Calculate the AUC
  roc_curve <- roc(yTest, yProba_vector)
  auc_value <- auc(roc_curve)
  
# Perform test directly on ROC curve
  test_result <- roc_curve$auc.p.value
  ci_intervals = ci.auc(roc_curve)
  
  #calculate nagelkerke r squared
  n_rsquared = calculate_nagelkerke_rsquared(yTest,yProba_vector)
  
  # Convert CI results to a data frame
  model_scores <- data.frame(
    auc = c(auc_value[1]),
    ci_lower = c(ci_intervals[1]),
    ci_upper = c(ci_intervals[3]),
    nagelkerke_rsquared = c(n_rsquared),
    model = c(data_type)
  )
  
  # Check if the file exists
  #columns : c("auc","ci_lower","ci_upper","nagelkerke_rquared","model")
  if (!file.exists(model_scores_file)) {
    
    print('creating model scores file ....')
    
    # If the file doesn't exist, create it with column names
    write.table(model_scores, file = model_scores_file, row.names = FALSE, sep = ",")
    cat('creating model scores file with scored features....\n')
    
  } else {
    # If the file exists, append the new data without column names
    write.table(model_scores, file = model_scores_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
    
    cat('writing to model scores file ......\n')
  }
  
}



