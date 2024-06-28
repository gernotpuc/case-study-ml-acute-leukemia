############################################################
# Calculate AUROC and CI 95% using the pretrained AIPAL model
# The following processing steps are performed:
# 1. Define input list of laboratory values
# 2. Assign paths to CSV-files of AL types or differential diagnosis
# 3. Preprocess and merge datasets
# 3. Load pretrained AIPAL-model
# 4. Perform predictions
# 5. Define confidence cutoffs
# 6. Calculate ans summarize rates of confident positive and negative diagnoses
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

# Initialize an empty list to store predictions and true labels
predictions_list <- list()
true_labels <- numeric()

# Loop over each encounter_reference to make predictions
for (encounter_id in data_filtered$encounter_reference) {
  # Filter the data for the current encounter_id
  input_df <- data_filtered %>%
    filter(encounter_reference == encounter_id) %>%
    select(-encounter_reference)

  # Extract the true label
  true_label <- ifelse(input_df$diagnosis == 'APL', 1, 0)
  true_labels <- c(true_labels, true_label)

  # Remove the diagnosis from input data
  input_df <- input_df %>%
    select(-diagnosis)

  # Ensure all required columns from the input list are present
  for (col in input) {
    if (!col %in% colnames(input_df)) {
      input_df[[col]] <- input[[col]]
    }
  }

  # Convert to data frame
  input_df <- as.data.frame(input_df)

  # Make prediction for the current encounter_id
  prediction <- predict(model, newdata = input_df, type = "prob", na.action = "na.pass")

  # Store the prediction with the encounter_id
  predictions_list[[encounter_id]] <- prediction
}

library(dplyr)

# Convert predictions_list into a dataframe
predictions <- do.call(rbind, lapply(names(predictions_list), function(id) {
  pred <- predictions_list[[id]]
  data.frame(encounter_reference = id, pred)
}))

# These are the original cutoff values for confidence as defined by AIPAL
cutoffs <- data.frame(
  category = c("ALL", "AML", "APL"),
  PPV = c(0.9593636, 0.9448341, 0.7487954),
  NPV = c(0.04462764, 0.02712490, 0.04154142),
  ACC = c(0.4528432, 0.5001304, 0.3770665)
)

# Merge predictions with diagnosis column
diagnosis <- data_filtered[, c("encounter_reference", "diagnosis")]
predictions_diagnosis <- merge(predictions, diagnosis, by = "encounter_reference", all.x = TRUE)

# Add confidence columns for each AL type based on cutoff values
predictions_confidence <- predictions_diagnosis %>%
  mutate(

    confident_APL = ifelse(APL >= cutoffs$PPV[cutoffs$category == "APL"], TRUE, FALSE),
    confident_not_APL = ifelse(APL <= cutoffs$NPV[cutoffs$category == "APL"], TRUE, FALSE),
    confident_AML = ifelse(AML >= cutoffs$PPV[cutoffs$category == "AML"], TRUE, FALSE),
    confident_not_AML = ifelse(AML <= cutoffs$NPV[cutoffs$category == "AML"], TRUE, FALSE),
    confident_ALL = ifelse(ALL >= cutoffs$PPV[cutoffs$category == "ALL"], TRUE, FALSE),
    confident_not_ALL = ifelse(ALL <= cutoffs$NPV[cutoffs$category == "ALL"], TRUE, FALSE),
    # Determine if any confident prediction was made
    confident = confident_APL | confident_not_APL | confident_AML | confident_not_AML | confident_ALL | confident_not_ALL
  ) %>%
  # Determine uncertain predictions
  mutate(predicted_class = case_when(
    confident_not_APL ~ "APL",
    confident_not_AML ~ "AML",
    confident_not_ALL ~ "ALL",
    confident_APL ~ "APL",
    confident_AML ~ "AML",
    confident_ALL ~ "ALL",
    TRUE ~ "Uncertain"
  ))

# Summarize the confidence results grouped by AL type
result <- predictions_confidence %>%
  group_by(diagnosis) %>%
  summarise(
    confident_APL = mean(confident_APL) * 100,
    confident_not_APL = mean(confident_not_APL) * 100,
    confident_AML = mean(confident_AML) * 100,
    confident_not_AML = mean(confident_not_AML) * 100,
    confident_ALL = mean(confident_ALL) * 100,
    confident_not_ALL = mean(confident_not_ALL) * 100,
    uncertain = mean(predicted_class == "Uncertain") * 100
  )

print(result)

