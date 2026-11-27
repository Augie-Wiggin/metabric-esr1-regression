# Load the complete-case data saved by Script 1.
model_data <- readRDS(
  file.path("data", "processed", "metabric_esr1_complete_cases.rds")
)

# Check the number of samples, variables, and their types.
dim(model_data)
names(model_data)
str(model_data)

# Fit the reference model with all eight predictors.
reference_model <- lm( ESR1_expression ~
    ER_IHC + PR_STATUS + HER2_STATUS +
    AGE_AT_DIAGNOSIS + GRADE + TUMOR_SIZE + 
    LYMPH_NODES_EXAMINED_POSITIVE + CELLULARITY,
    data = model_data
)

# Review the model results.
summary(reference_model)

# Use the stepAIC() function from the MASS package for backward selection.
aic_model <- MASS::stepAIC(
  reference_model,
  direction = "backward"
)

# Review the model selected by AIC.
summary(aic_model)

# Compare the AIC values of the two models.
AIC(reference_model, aic_model)

# Save the reference and AIC-selected models for later scripts.
saveRDS(
  reference_model,
  file.path("data", "processed", "metabric_esr1_reference_model.rds")
)

saveRDS(
  aic_model,
  file.path("data", "processed", "metabric_esr1_aic_model.rds")
)

# Load glmnet, the package used for LASSO regression.
library(glmnet)

# Convert the predictors into numeric columns for LASSO.
lasso_x <- model.matrix(ESR1_expression ~ ., data = model_data)[, -1]

# Store ESR1 expression as the response.
lasso_y <- model_data$ESR1_expression

# Check the LASSO inputs.
dim(lasso_x)
length(lasso_y)

# Set a seed so the fold assignment can be reproduced.
set.seed(123)
lasso_folds <- sample(rep(1:10, length.out = nrow(lasso_x)))

# Fit LASSO and use ten-fold cross-validation to choose its penalty.
lasso_model <- cv.glmnet(
  lasso_x,
  lasso_y,
  alpha = 1,
  foldid = lasso_folds,
  standardize = TRUE,
  type.measure = "mse"
)

# Display the penalty with the lowest cross-validation MSE.
lasso_model$lambda.min

# Display the LASSO coefficients at the selected penalty.
lasso_coefficients <- coef(lasso_model, s = "lambda.min")
lasso_coefficients

# Save the LASSO cross-validation plot.
png(
  filename = file.path("figures", "esr1_lasso_cv.png"),
  width = 1200,
  height = 900,
  res = 150
)

plot(lasso_model)
title("LASSO Cross-Validation", line = 3)

dev.off()

# Save the fitted LASSO model for later use.
saveRDS(
  lasso_model,
  file.path("data", "processed", "metabric_esr1_lasso_model.rds")
)