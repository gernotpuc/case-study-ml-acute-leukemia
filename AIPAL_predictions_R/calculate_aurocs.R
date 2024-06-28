############################################################
# Calculate AUROC and CI 95% using the pretrained AIPAL model
# The following processing steps are performed:
# 1. Define input list of laboratory values
# 2. Assign paths to CSV-files of AL types or differential diagnosis
# 3. Preprocess and merge datasets
# 3. Load pretrained AIPAL-model
# 4. Perform predictions
# 5. Perform bootstrapping for calculating AUROC with CI
############################################################


# Load required libraries
library(dplyr)
library(tidyr)
library(pROC)

# The following inputs are expected by the model.
# The labels in the column aipal_labels must correspond to these expected inputs.
input <- c(
  'MCV_fL',
  'PT_percent',
  'LDH_UI_L',
  'MCHC_g_L',
  'WBC_G_L',
  'Fibrinogen_g_L',
  'age',
  'Monocytes_G_L',
  'Platelets_G_L',
  'Lymphocytes_G_L',
  'Monocytes_percent'
)

# Read the CSV file

data_all <- read.csv("/PATH_TO_CSV/R_aipal_ALL.csv")
data_aml <- read.csv("/PATH_TO_CSV/R_aipal_AML.csv")
data_apl <- read.csv("/PATH_TO_CSV/R_aipal_APL.csv")
# Assign diagnosis classes if not inlcuded in the datasets already.
# Replace with differential diagnosis eventually.
data_all$diagnosis <- 'ALL'
data_aml$diagnosis <- 'AML'
data_apl$diagnosis <- 'APL'

# Merge datasets and convert values of labortatory values to numeric
data <- bind_rows(data_all, data_aml, data_apl)
data$value_quantity <- as.numeric(data$value_quantity)

# Group by encounter_reference and spread the value_quantity based on the input labels of laboratory values
data_wide <- data %>%
  group_by(encounter_reference) %>%
  spread(key = aipal_labels, value = value_quantity)

# Extract only the necessary columns based on the input list
data_filtered <- data_wide %>%
  select(encounter_reference, one_of(input), diagnosis) %>%
  group_by(encounter_reference) %>%
  summarize_all(~first(na.omit(.)))

# Convert age to numeric
data_filtered$age <- as.numeric(data_filtered$age)

# Load the model
res_list <- readRDS("data/221003_Final_model_res_list.rds")
model <- res_list$final_model

# Initialize lists to store predictions and true labels for each class
predictions_list <- list(ALL = numeric(), AML = numeric(), APL = numeric())
true_labels_list <- list(ALL = numeric(), AML = numeric(), APL = numeric())


# Loop over each encounter_reference to make predictions
for (encounter_id in data_filtered$encounter_reference) {
  # Filter the data for the current encounter_id
  input_df <- data_filtered %>%
    filter(encounter_reference == encounter_id) %>%
    select(-encounter_reference)

  # Extract the true label for each class
  for (class in c("ALL", "AML", "APL")) {
    true_label <- ifelse(input_df$diagnosis == class, 1, 0)
    true_labels_list[[class]] <- c(true_labels_list[[class]], true_label)
  }

  # Remove the true label from input data
  input_df <- input_df %>%
    select(-diagnosis)

  # Check if all required columns from the input list are present
  for (col in input) {
    if (!col %in% colnames(input_df)) {
      input_df[[col]] <- input[[col]]
    }
  }

  # Convert to data frame
  input_df <- as.data.frame(input_df)

  # Make prediction for the current encounter_id
  prediction <- predict(model, newdata = input_df, type = "prob", na.action = "na.pass")

  # Store the prediction with the encounter_id for each class
  for (class in c("ALL", "AML", "APL")) {
    predictions_list[[class]] <- c(predictions_list[[class]], prediction[[class]])
  }
}


# Function to calculate AUROC with bootstrapping and return the AUROC and confidence intervals
calculate_auroc_with_ci <- function(true_labels, probabilities, num_bootstrap_samples = 2000, ci_level = 0.95) {
  options(warn = -1)

  # Calculate AUROC on original data
  suppressMessages(roc_obj <- roc(true_labels, probabilities))
  auroc_original <- auc(roc_obj)

  # Perform bootstrapping
  auroc_bootstrapped <- numeric(num_bootstrap_samples)
  for (i in 1:num_bootstrap_samples) {
    # Resample with replacement
    indices <- sample(length(true_labels), replace = TRUE)
    bootstrapped_true_labels <- true_labels[indices]
    bootstrapped_probabilities <- probabilities[indices]
    # Calculate AUROC on bootstrapped data
    suppressMessages(bootstrapped_roc_obj <- roc(bootstrapped_true_labels, bootstrapped_probabilities))
    auroc_bootstrapped[i] <- auc(bootstrapped_roc_obj)
  }

  # Compute confidence intervals
  ci_lower <- quantile(auroc_bootstrapped, (1 - ci_level) / 2)
  ci_upper <- quantile(auroc_bootstrapped, 1 - (1 - ci_level) / 2)
  return(list(auroc_original = auroc_original, ci_lower = ci_lower, ci_upper = ci_upper))
}


# Calculate AUROC for each class with confidence intervals
for (class in c("ALL", "AML", "APL")) {
  # Extract probabilities for the current class
  probabilities <- predictions_list[[class]]
  true_labels <- true_labels_list[[class]]

  # Calculate AUROC with confidence intervals
  result <- calculate_auroc_with_ci(true_labels, probabilities)

  # Display the AUROC and confidence intervals for the current class
  cat("AUROC for", class, "class:", result$auroc_original, "\n")
  cat("95% Confidence Interval:", result$ci_lower, "-", result$ci_upper, "\n")
}





