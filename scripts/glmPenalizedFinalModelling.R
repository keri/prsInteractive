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

registerDoMC(cores = 20)

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


get_geno_dataset <- function(data_file,data_path,columns_to_get){
  
  start_time <- Sys.time()
  
  # Adding 'IID' and 'PHENOTYPE' to the list of columns to get
  columns_to_get <- c('IID', 'PHENOTYPE', columns_to_get)
  
  # Load the main dataset
  main_data <- fread(data_file, sep=" ", select = columns_to_get)  
  
  # Load the file containing the values to match (assuming no header)
  withdrawal_path = paste0(data_path,"/withdrawals.csv")
  withdrawn <- fread(withdrawal_path, header = FALSE)
  
  # Rename the column in values_to_remove to 'IID' (or whatever the matching column is in main_data)
  setnames(withdrawn, "V1", "IID")
  
  # Remove withdrawn individuals from the DataFrame
  main_data <- main_data[!main_data$IID %in% withdrawn$IID]
  
  # Modify the 'PHENOTYPE' column in place
  main_data[PHENOTYPE == 1, PHENOTYPE := 0]
  main_data[PHENOTYPE == 2, PHENOTYPE := 1]
  
  end_time <- Sys.time()
  print(paste('time it took to download entire dataset is', round(difftime(end_time, start_time, units = "mins"), 2), 'minutes'))
  
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

get_important_features <- function(feature_file){

  #file has 3 columns : [feature,data_type:(main,epi)]
# feature_file = paste0(feature_pathway,'/importantFeaturesForAssociationAnalysis.csv')
  #feature_file = paste0(feature_pathway,'/importantFeaturesPostShap.csv')
  important_features = read.csv(feature_file)

  epi_main_features = c(important_features$feature)

  main_features = important_features$feature[important_features$data_type == "main"]

  epi_features = important_features$feature[important_features$data_type == "epi"]

  return(list(epi_main_features=epi_main_features,main_features=main_features,epi_features=epi_features))
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
# parser <- ArgumentParser()
# parser$add_argument("--results_path", required = TRUE)
# parser$add_argument("--data_path", required = TRUE)
# parser$add_argument("--hla_file", required = TRUE) 
# parser$add_argument("--covar_file", required = TRUE)
# parser$add_argument("--pheno_path", required = TRUE)
# parser$add_argument("--training_file", required = TRUE)
# parser$add_argument("--test_file", required = TRUE)
# parser$add_argument("--training_env_gen_file", required = TRUE)
# parser$add_argument("--test_env_gen_file", required = TRUE)
# parser$add_argument("--feature_model_file" ,required = TRUE)
# 
# args <- parser$parse_args()
# 
# results_path <- args$results_path
# data_path <- args$data_path
# pheno_path <- args$pheno_path
# hla_file <- args$hla_file
# training_file <- args$training_file
# test_file <- args$test_file
# training_env_gen_file <- args$training_env_gen_file
# test_env_gen_file <- args$test_env_gen_file
# covar_file <- args$covar_file
# feature_model_file <- args$feature_model_file


results_path <- '/Users/kerimulterer/prsInteractive/results'
pheno_path <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
data_path <- '/Users/kerimulterer/prsInteractive/data'
hla_file <- '/Users/kerimulterer/prsInteractive/results/participant_hla.csv'
# training_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw'
training_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw'
test_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/testCombined.raw'
# test_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/testCombined_final.raw'
training_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTraining.csv'
test_env_gen_file <- '/Users/kerimulterer/prsInteractive/results/type2Diabetes/geneEnvironmentTest.csv'
covar_file='/Users/kerimulterer/prsInteractive/results/covar.csv'
feature_model_file='/Users/kerimulterer/prsInteractive/results/type2Diabetes/scores/importantFeaturesPostShap.csv'

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
#feature_scores_file = paste0(scores_path,'/importantFeaturesForAssociationAnalysis.csv')
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
trainingDf = get_geno_dataset(training_file,data_path,columns_to_get)

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
  list('main',c(main_features,covariate_columns,hla_columns)),
  list('epi',c(epi_features,covariate_columns,hla_columns)),
  list('epi+main',c(epi_main_features,covariate_columns,hla_columns)),
  list('cardio',c(epi_cardio_features,covariate_columns,hla_columns)),
  list('all',c(epi_cardio_features,epi_main_features,hla_columns,covariate_columns)),
  list('cardio_main',c(main_cardio_features,covariate_columns)),
  list('all+main_cardio',c(epi_cardio_features,epi_main_features,hla_columns,main_cardio_features,covariate_columns)),
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
  yProba$IID = trainingDf$IID
  
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



