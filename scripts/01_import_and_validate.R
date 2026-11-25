# Confirm the source files are available and preview their structure.

# Set paths to the raw files.
raw_dir <- "data/raw"

patient_file <- file.path(raw_dir, "data_clinical_patient.txt")
sample_file <- file.path(raw_dir, "data_clinical_sample.txt")
expression_file <- file.path(raw_dir, "data_mrna_illumina_microarray.txt")

# Confirm that all three files are present.
input_files <- c(patient_file, sample_file, expression_file)
file.exists(input_files)

# Preview the clinical file headers.
readLines(patient_file, n = 6)
readLines(sample_file, n = 6)

# Import the clinical files
# The first four lines describe the fields; column names start on line 5.

clinical_patient <- read.delim(
  patient_file,
  skip = 4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clinical_sample <- read.delim(
  sample_file,
  skip = 4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check the imported tables.
dim(clinical_patient)
dim(clinical_sample)

names(clinical_patient)[1:6]
names(clinical_sample)[1:6]

head(clinical_sample[, c("PATIENT_ID", "SAMPLE_ID")])

# Check whether each identifier is unique.
sum(duplicated(clinical_patient$PATIENT_ID))
sum(duplicated(clinical_sample$SAMPLE_ID))
sum(duplicated(clinical_sample$PATIENT_ID))

# Check whether both clinical files cover the same patients.
setequal(clinical_patient$PATIENT_ID, clinical_sample$PATIENT_ID)

# Find ER-related fields, including any IHC field.
grep("(^ER|IHC)", names(clinical_patient), value = TRUE, ignore.case = TRUE)
grep("(^ER|IHC)", names(clinical_sample), value = TRUE, ignore.case = TRUE)

# Preview a few rows of the mRNA expression matrix.
expression_preview <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  nrows = 3,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(expression_preview)
names(expression_preview)[1:6]
head(expression_preview[, 1:6])

# Read only the gene-symbol column.
gene_symbols <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  colClasses = c("character", rep("NULL", ncol(expression_preview) - 1)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Find the row labeled ESR1.
esr1_row <- which(gene_symbols$Hugo_Symbol == "ESR1")

esr1_row
gene_symbols[esr1_row, ]

# Read the ESR1 row from the expression matrix.
esr1_data <- read.delim(
  expression_file,
  header = FALSE,
  sep = "\t",
  skip = esr1_row,
  nrows = 1,
  col.names = names(expression_preview),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check the gene identifiers and a few sample values.
esr1_data[, 1:8]

# Put ESR1 expression values into a sample-level table.
esr1_by_sample <- data.frame(
  SAMPLE_ID = names(esr1_data)[3:ncol(esr1_data)],
  ESR1_expression = as.numeric(esr1_data[1, 3:ncol(esr1_data)]),
  stringsAsFactors = FALSE)

# Check a few sample IDs and expression values.
head(esr1_by_sample)
length(esr1_by_sample$SAMPLE_ID)
sum(esr1_by_sample$SAMPLE_ID %in% clinical_sample$SAMPLE_ID)
sum(duplicated(esr1_by_sample$SAMPLE_ID))

# Match ESR1 expression values to sample-level clinical information.
expression_clinical <- merge(
  esr1_by_sample,
  clinical_sample,
  by = "SAMPLE_ID")

# Add patient-level clinical information.
analysis_data <- merge(
  expression_clinical,
  clinical_patient,
  by = "PATIENT_ID")

# Check the joined dataset
dim(analysis_data)
head(analysis_data[, c("PATIENT_ID", "SAMPLE_ID", "ESR1_expression", "ER_IHC")])

# Check ER_IHC categories before correcting the typo.
table(analysis_data$ER_IHC, useNA = "ifany")

# Correct the source spelling.
analysis_data$ER_IHC[analysis_data$ER_IHC == "Positve"] <- "Positive"

# Count blank ER_IHC entries.
sum(analysis_data$ER_IHC == "", na.rm = TRUE)

# Treat blank entries as missing values.
analysis_data$ER_IHC[analysis_data$ER_IHC == ""] <- NA

# Check the ER_IHC counts, including missing values.
table(analysis_data$ER_IHC, useNA = "ifany")

names(analysis_data)

# The expression and clinical records have been joined, and ER_IHC values have
# been cleaned. This section checks missingness in the response and planned
# predictor pool. These counts will help us assess which variables are suitable
# for the final models; they do not mean every candidate will be included.

candidate_variables <- c(
  "ESR1_expression",
  "ER_IHC",
  "AGE_AT_DIAGNOSIS",
  "GRADE",
  "TUMOR_SIZE",
  "LYMPH_NODES_EXAMINED_POSITIVE",
  "CELLULARITY",
  "COHORT"
)

count_missing <- function(x) {
  if (is.character(x)) {
    sum(is.na(x) | trimws(x) == "")
  } else {
    sum(is.na(x))
  }
}

sapply(analysis_data[candidate_variables], count_missing)
str(analysis_data[candidate_variables])

# Treat blank cellularity entries as missing.
analysis_data$CELLULARITY[
  !is.na(analysis_data$CELLULARITY) &
    trimws(analysis_data$CELLULARITY) == ""
] <- NA

# Check categories and their counts.
table(analysis_data$GRADE, useNA = "ifany")
table(analysis_data$CELLULARITY, useNA = "ifany")
table(analysis_data$COHORT, useNA = "ifany")

# Summarize numeric candidate variables.
summary(analysis_data[c(
  "AGE_AT_DIAGNOSIS",
  "TUMOR_SIZE",
  "LYMPH_NODES_EXAMINED_POSITIVE"
)])

table(analysis_data$GRADE, useNA = "ifany")
table(analysis_data$CELLULARITY, useNA = "ifany")
table(analysis_data$COHORT, useNA = "ifany")
sum(complete.cases(analysis_data[candidate_variables]))

# Create the complete-case modeling dataset.
model_data <- analysis_data[
  complete.cases(analysis_data[candidate_variables]),
  candidate_variables]

# Treat categorical predictors as factors.
model_data$ER_IHC <- factor(model_data$ER_IHC)
model_data$GRADE <- factor(model_data$GRADE)
model_data$CELLULARITY <- factor(model_data$CELLULARITY)
model_data$COHORT <- factor(model_data$COHORT)

# Confirm the size and variable types.
dim(model_data)
str(model_data)

# Save the cleaned data for use in later scripts.
saveRDS(
  analysis_data,
  file.path("data", "processed", "metabric_esr1_joined.rds"))

saveRDS(
  model_data,
  file.path("data", "processed", "metabric_esr1_complete_cases.rds"))