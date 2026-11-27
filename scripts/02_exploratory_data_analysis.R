# This script explores ESR1 expression and its relationships with the
# planned clinical and tumor predictors.

# Load the prepared datasets from script 1.
analysis_data <- readRDS(file.path(
  "data", "processed", "metabric_esr1_joined.rds"))

model_data <- readRDS(file.path(
  "data", "processed", "metabric_esr1_complete_cases.rds"))

# Confirm the dataset sizes.
dim(analysis_data)
dim(model_data)

# Describe the complete-case analysis sample.
summary(model_data)

# Examine ESR1 expression across all profiled samples.
summary(analysis_data$ESR1_expression)

# Save a histogram of ESR1 expression.
png(
  filename = file.path("figures", "esr1_expression_distribution.png"),
  width = 1200,
  height = 900,
  res = 150
)

hist(
  analysis_data$ESR1_expression,
  main = "Distribution of ESR1 Expression",
  xlab = "ESR1 mRNA expression (log2 intensity)",
  ylab = "Number of samples",
  col = "lightblue",
  border = "white"
)

dev.off()

# Compare ESR1 expression across the IHC-based ER groups.
er_counts <- table(model_data$ER_IHC)
er_counts

# Summarize ESR1 expression separately for each ER-IHC group.
tapply(model_data$ESR1_expression, model_data$ER_IHC, summary)

# Add group sample sizes to the plot labels.
er_labels <- paste0(names(er_counts), "\n(n = ", as.integer(er_counts), ")")

png(
  filename = file.path("figures", "esr1_by_er_ihc.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  ESR1_expression ~ ER_IHC,
  data = model_data,
  main = "ESR1 Expression by IHC ER Status",
  xlab = "IHC ER status",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue")
)

# Add ER-group labels and sample sizes without axis tick marks.
axis(1, at = 1:2, labels = er_labels, tick = FALSE)

dev.off()

# Compare ESR1 expression across the PR-status groups.
pr_counts <- table(model_data$PR_STATUS)
pr_counts

# Summarize ESR1 expression separately for each PR-status group.
tapply(model_data$ESR1_expression, model_data$PR_STATUS, summary)

# Add group sample sizes to the plot labels.
pr_labels <- paste0(names(pr_counts), "\n(n = ", as.integer(pr_counts), ")")

# Save a boxplot of ESR1 expression by PR status.
png(
  filename = file.path("figures", "esr1_by_pr_status.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  ESR1_expression ~ PR_STATUS,
  data = model_data,
  main = "ESR1 Expression by PR Status",
  xlab = "PR status",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue")
)

# Add group labels and sample sizes without tick marks.
axis(1, at = 1:2, labels = pr_labels, tick = FALSE)

dev.off()

# Compare ESR1 expression across the HER2-status groups.
her2_counts <- table(model_data$HER2_STATUS)
her2_counts

# Summarize ESR1 expression separately for each HER2-status group.
tapply(model_data$ESR1_expression, model_data$HER2_STATUS, summary)

# Add group sample sizes to the plot labels.
her2_labels <- paste0(
  names(her2_counts),
  "\n(n = ",
  as.integer(her2_counts),
  ")"
)

# Save a boxplot of ESR1 expression by HER2 status.
png(
  filename = file.path("figures", "esr1_by_her2_status.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  ESR1_expression ~ HER2_STATUS,
  data = model_data,
  main = "ESR1 Expression by HER2 Status",
  xlab = "HER2 status",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue")
)

# Add group labels and sample sizes without tick marks.
axis(1, at = 1:2, labels = her2_labels, tick = FALSE)

dev.off()

# Compare ESR1 expression across the tumor-grade groups.
grade_counts <- table(model_data$GRADE)
grade_counts

# Summarize ESR1 expression separately for each grade.
tapply(model_data$ESR1_expression, model_data$GRADE, summary)

# Add group sample sizes to the plot labels.
grade_labels <- paste0(
  names(grade_counts),
  "\n(n = ",
  as.integer(grade_counts),
  ")"
)

# Save a boxplot of ESR1 expression by tumor grade.
png(
  filename = file.path("figures", "esr1_by_grade.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  ESR1_expression ~ GRADE,
  data = model_data,
  main = "ESR1 Expression by Tumor Grade",
  xlab = "Neoplasm histologic grade",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue", "lightpink")
)

# Add grade labels and sample sizes without tick marks.
axis(1, at = 1:3, labels = grade_labels, tick = FALSE)

dev.off()

# Compare ESR1 expression across the cellularity groups.
cellularity_counts <- table(model_data$CELLULARITY)
cellularity_counts

# Summarize ESR1 expression separately for each cellularity group.
tapply(model_data$ESR1_expression, model_data$CELLULARITY, summary)

# Add group sample sizes to the plot labels.
cellularity_labels <- paste0(
  names(cellularity_counts),
  "\n(n = ",
  as.integer(cellularity_counts),
  ")"
)

# Save a boxplot of ESR1 expression by cellularity.
png(
  filename = file.path("figures", "esr1_by_cellularity.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  ESR1_expression ~ CELLULARITY,
  data = model_data,
  main = "ESR1 Expression by Tumor Cellularity",
  xlab = "Tumor content",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue", "lightpink")
)

# Add group labels and sample sizes without tick marks.
axis(1, at = 1:3, labels = cellularity_labels, tick = FALSE)

dev.off()

# Save a scatterplot of ESR1 expression by age at diagnosis.
png(
  filename = file.path("figures", "esr1_by_age.png"),
  width = 1200,
  height = 900,
  res = 150
)

plot(
  ESR1_expression ~ AGE_AT_DIAGNOSIS,
  data = model_data,
  main = "ESR1 Expression by Age at Diagnosis",
  xlab = "Age at diagnosis",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  pch = 16,
  cex = 0.5,
  col = "steelblue"
)

dev.off()

# Save a scatterplot of ESR1 expression by tumor size.
png(
  filename = file.path("figures", "esr1_by_tumor_size.png"),
  width = 1200,
  height = 900,
  res = 150
)

plot(
  ESR1_expression ~ TUMOR_SIZE,
  data = model_data,
  main = "ESR1 Expression by Tumor Size",
  xlab = "Tumor size",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  pch = 16,
  cex = 0.5,
  col = "steelblue"
)

dev.off()

# Save a scatterplot of ESR1 expression by positive lymph-node count.
png(
  filename = file.path("figures", "esr1_by_positive_lymph_nodes.png"),
  width = 1200,
  height = 900,
  res = 150
)

plot(
  ESR1_expression ~ LYMPH_NODES_EXAMINED_POSITIVE,
  data = model_data,
  main = "ESR1 Expression by Positive Lymph-Node Count",
  xlab = "Number of positive lymph nodes",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  pch = 16,
  cex = 0.5,
  col = "steelblue"
)

dev.off()

# Compare ER-IHC status with PR status.
er_pr_counts <- table(model_data$ER_IHC, model_data$PR_STATUS)
er_pr_counts
round(prop.table(er_pr_counts, margin = 1), 3)

# Compare ER-IHC status with HER2 status.
er_her2_counts <- table(model_data$ER_IHC, model_data$HER2_STATUS)
er_her2_counts
round(prop.table(er_her2_counts, margin = 1), 3)

# Check correlations among the numeric candidate predictors.
numeric_predictors <- model_data[c(
  "AGE_AT_DIAGNOSIS",
  "TUMOR_SIZE",
  "LYMPH_NODES_EXAMINED_POSITIVE"
)]

round(cor(numeric_predictors), 2)

# Mark samples as included or excluded using the model-data variables.
# Classify each row by whether it has complete values for all model variables.
case_status <- ifelse(
  complete.cases(analysis_data[names(model_data)]),
  "Included",
  "Excluded"
)

# Count samples and summarize ESR1 expression in each group.
case_counts <- table(case_status)
case_counts

tapply(analysis_data$ESR1_expression, case_status, summary)

# Add sample sizes to the plot labels.
case_labels <- paste0(
  names(case_counts),
  "\n(n = ",
  as.integer(case_counts),
  ")"
)

# Save a boxplot comparing ESR1 expression by complete-case status.
png(
  filename = file.path("figures", "esr1_by_case_status.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  analysis_data$ESR1_expression ~ case_status,
  main = "ESR1 Expression by Complete-Case Status",
  xlab = "Model analysis sample",
  ylab = "ESR1 mRNA expression (log2 intensity)",
  xaxt = "n",
  col = c("lightgray", "lightblue")
)

# Add group labels and sample sizes without tick marks.
axis(1, at = 1:2, labels = case_labels, tick = FALSE)

dev.off()