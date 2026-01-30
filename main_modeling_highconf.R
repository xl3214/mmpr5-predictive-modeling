# ==============================================================================
# RESISTANCE PREDICTION MODELING: BASELINE VS ENHANCED
# ==============================================================================
#
# Title:        Comparison of Baseline and Structurally-Enhanced Models for
#               Bedaquiline Resistance Prediction in mmpR5 Mutants
# Author:       Xuan Lu
# Date:         01/30/2026
#
# Description:  This script compares baseline (mutation category + frameshift)
#               and enhanced (+ structural ΔΔG features) models for predicting
#               bedaquiline resistance. Both logistic regression and random
#               forest approaches are evaluated.
#
# Input:        data_highconf_singleocc.csv (high-confidence dataset)
#
# Output:       - Tables: model_comparison_logistic.csv
#                         model_comparison_rf.csv
#                         model_comparison_all.csv
#                         coefficient_table.csv
#                         variable_importance.csv
#               - Figures: fig_logistic_coefficients.png
#                          fig_logistic_predictions.png
#                          fig_rf_importance.png
#                          fig_roc_comparison.png
#                          fig_method_comparison.png
#
# Methods:      - Logistic regression with likelihood ratio test
#               - Random forest (500 trees) with permutation importance
#               - ROC/AUC analysis for model comparison
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(car)
library(pdp)
library(broom)
library(patchwork)

# Set seed for reproducibility
set.seed(2025)
set.seed(2025)

# Define file paths
INPUT_FILE <- "data_highconf_singleocc.csv"
OUTPUT_DIR <- "."

# Define feature groups
RAV_FEATURE <- "mmpR5"
STRUCTURAL_FEATURES <- c("MAESTRO", "mCSM_DNA", "mCSM_protein", 
                         "mCSM_ligand", "d_hydrophobicity", "g_score")
OUTCOME <- "S.R"

# ------------------------------------------------------------------------------
# 1. DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------

cat("=== 1. DATA LOADING AND PREPROCESSING ===\n\n")

# Load data
df <- read.csv(INPUT_FILE)
colnames(df) <- make.names(colnames(df))

# Convert outcome to ordered factor
df$S.R <- factor(df$S.R, levels = c("S", "R"), ordered = TRUE)

# Prepare mmpR5 as factor with all levels
all_mmpR5_levels <- unique(df[[RAV_FEATURE]])
df[[RAV_FEATURE]] <- factor(df[[RAV_FEATURE]], levels = all_mmpR5_levels)

# Remove missing outcome values
df <- df[!is.na(df$S.R), ]

cat("Data loaded:", nrow(df), "samples\n")
cat("Unique mmpR5 variants:", length(all_mmpR5_levels), "\n")

# ------------------------------------------------------------------------------
# 2. MUTATION DOMAIN ANNOTATION
# ------------------------------------------------------------------------------

cat("\n=== 2. MUTATION DOMAIN ANNOTATION ===\n\n")

# Define mutations by structural domain
dna_binding <- c(
  '132insG', '132insGG', '138insG', '140insGG', '141insG', '192insG', '211insG', 
  '234insG', 'A59V', 'A62T', 'A62V', 'A84V', 'C46G', 'C46R', 'C46Y', 'E55D', 
 'G41V', 'G65R', 'I67L', 'I67M', 'I67S', 'I67V', 'I80M', 'I80S', 'L40F', 
  'L40M', 'L40S', 'L40V', 'L74M', 'L74V', 'L83F', 'L83P', 'M73I', 'N70D', 
  'P48H', 'P48L', 'Q51K', 'Q51R', 'Q76E', 'R50Q', 'R72W', 'S52F', 'S52P', 
  'S53P', 'S63G', 'S63R', 'S68G', 'S68I', 'S68N', 'S68R', 'T58P', 'V85A', 
  'W42C', 'W42R'
)

dna_non_binding <- c(
  '107insG', '274insG', '28insG', '324insG', '418insG', '464insG', '465insG', 
  '492insGG', 'A101T', 'A102V', 'A118T', 'A36T', 'A36V', 'A86T', 'A86V', 
  'A99P', 'D15G', 'D5G', 'D5N', 'D88G', 'E104G', 'E113K', 'E147D', 'E147G', 
  'E21D', 'E21K', 'F100Y', 'F19S', 'F79S', 'F93L', 'F93S', 'G103S', 'G121R', 
  'G162E', 'G25C', 'G87A', 'G87R', 'I108S', 'I108T', 'I108V', 'L114P', 'L117R', 
  'L125R', 'L136P', 'L142R', 'L154R', 'M10I', 'M139I', 'M146T', 'M17V', 'M23V', 
  'N4T', 'N98D', 'P14L', 'P14S', 'Q22R', 'R105C', 'R107C', 'R109P', 'R109Q', 
  'R123K', 'R134G', 'R135W', 'R137P', 'R156Q', 'R90C', 'R90L', 'R94Q', 'R94W', 
  'R96G', 'R96L', 'R96Q', 'R96T', 'R96W', 'S158R', 'S2I', 'S2R', 'S31R', 
  'T33S', 'V20A', 'V20G', 'Y145H', 'Y26C', 'L32S', 'Y92D', 'M139T', 'D8G', 
  'A101E', 'L39S'
)

# Annotate mutations
df$dna_bind <- NA
df$dna_bind[df[[RAV_FEATURE]] %in% dna_binding] <- "DNA_binding"
df$dna_bind[df[[RAV_FEATURE]] %in% dna_non_binding & is.na(df$dna_bind)] <- "DNA_non_binding"
df$dna_bind <- factor(df$dna_bind, levels = c("DNA_non_binding", "DNA_binding"))

# Create frameshift indicator
df$FS <- ifelse(grepl("ins", df[[RAV_FEATURE]]), 1, 0)

cat("Mutation categories assigned.\n")

# ------------------------------------------------------------------------------
# 3. DATASET SUMMARY
# ------------------------------------------------------------------------------

cat("\n=== 3. DATASET SUMMARY ===\n\n")

cat("Total samples:", nrow(df), "\n")
cat("Frameshift mutations:", sum(df$FS == 1), "\n")
cat("Point mutations:", sum(df$FS == 0), "\n")
cat("Resistant (R):", sum(df[[OUTCOME]] == "R"), "\n")
cat("Susceptible (S):", sum(df[[OUTCOME]] == "S"), "\n")
cat("Resistance rate:", round(mean(df[[OUTCOME]] == "R"), 3), "\n")

# Category distribution
category_summary <- df %>%
  group_by(dna_bind, FS) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))

cat("\nDistribution by mutation category:\n")
print(category_summary)

# Prepare complete cases
complete_vars <- c(OUTCOME, RAV_FEATURE, STRUCTURAL_FEATURES, "dna_bind", "FS")
complete_cases <- complete.cases(df[, complete_vars])
df_complete <- df[complete_cases, ]

cat("\nComplete cases for modeling:", nrow(df_complete), "/", nrow(df), "\n")

# ------------------------------------------------------------------------------
# 4. MODEL SPECIFICATIONS
# ------------------------------------------------------------------------------

cat("\n=== 4. MODEL SPECIFICATIONS ===\n\n")

# Define model formulas
baseline_formula <- as.formula("S.R ~ dna_bind + FS")
enhanced_formula <- as.formula(paste("S.R ~ dna_bind + FS +", 
                                     paste(STRUCTURAL_FEATURES, collapse = " + ")))

cat("Baseline model: S.R ~ dna_bind + FS\n")
cat("Enhanced model: S.R ~ dna_bind + FS +", paste(STRUCTURAL_FEATURES, collapse = " + "), "\n")

# ------------------------------------------------------------------------------
# 5. LOGISTIC REGRESSION
# ------------------------------------------------------------------------------

cat("\n=== 5. LOGISTIC REGRESSION ===\n\n")

# --- 5.1 Fit Models ---
null_model <- glm(S.R ~ 1, data = df_complete, family = binomial())
baseline_logit <- glm(baseline_formula, data = df_complete, family = binomial())
enhanced_logit <- glm(enhanced_formula, data = df_complete, family = binomial())

cat("Models fitted successfully.\n")

# --- 5.2 Model Fit Statistics ---
model_comparison_logit <- data.frame(
  Model = c("Baseline", "Enhanced"),
  AIC = c(AIC(baseline_logit), AIC(enhanced_logit)),
  BIC = c(BIC(baseline_logit), BIC(enhanced_logit)),
  LogLik = c(as.numeric(logLik(baseline_logit)), as.numeric(logLik(enhanced_logit))),
  Deviance = c(deviance(baseline_logit), deviance(enhanced_logit)),
  df = c(length(coef(baseline_logit)), length(coef(enhanced_logit))),
  McFadden_R2 = c(
    as.numeric(1 - (logLik(baseline_logit) / logLik(null_model))),
    as.numeric(1 - (logLik(enhanced_logit) / logLik(null_model)))
  )
)

model_comparison_logit$Delta_AIC <- model_comparison_logit$AIC - min(model_comparison_logit$AIC)
model_comparison_logit$Delta_BIC <- model_comparison_logit$BIC - min(model_comparison_logit$BIC)

cat("\n--- Logistic Regression Model Comparison ---\n")
print(model_comparison_logit)

# --- 5.3 Likelihood Ratio Test ---
cat("\n--- Likelihood Ratio Test ---\n")
lr_test <- anova(baseline_logit, enhanced_logit, test = "Chisq")
print(lr_test)

# --- 5.4 Coefficient Analysis ---
enhanced_coef <- coef(enhanced_logit)
enhanced_ci <- confint(enhanced_logit)
enhanced_summary <- summary(enhanced_logit)

coef_table <- data.frame(
  Term = names(enhanced_coef),
  Estimate = enhanced_coef,
  SE = enhanced_summary$coefficients[, "Std. Error"],
  z_value = enhanced_summary$coefficients[, "z value"],
  p_value = enhanced_summary$coefficients[, "Pr(>|z|)"],
  CI_lower = enhanced_ci[, 1],
  CI_upper = enhanced_ci[, 2],
  OR = exp(enhanced_coef),
  OR_lower = exp(enhanced_ci[, 1]),
  OR_upper = exp(enhanced_ci[, 2]),
  row.names = NULL
)

cat("\n--- Enhanced Model Coefficients ---\n")
print(coef_table, digits = 3)

# --- 5.5 Prediction Performance ---

# Helper function for metrics
calc_metrics <- function(cm, actual_levels) {
  accuracy <- sum(diag(cm)) / sum(cm)
  sensitivity <- if ("R" %in% rownames(cm) && "R" %in% colnames(cm)) {
    cm["R", "R"] / sum(cm[, "R"])
  } else 0
  specificity <- if ("S" %in% rownames(cm) && "S" %in% colnames(cm)) {
    cm["S", "S"] / sum(cm[, "S"])
  } else 0
  list(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity)
}

# Generate predictions
baseline_pred <- predict(baseline_logit, type = "response")
enhanced_pred <- predict(enhanced_logit, type = "response")

baseline_class <- factor(ifelse(baseline_pred > 0.5, "R", "S"), levels = c("S", "R"))
enhanced_class <- factor(ifelse(enhanced_pred > 0.5, "R", "S"), levels = c("S", "R"))

# Confusion matrices
cat("\n--- Baseline Logistic Confusion Matrix ---\n")
baseline_cm <- table(Predicted = baseline_class, Actual = df_complete$S.R)
print(baseline_cm)

cat("\n--- Enhanced Logistic Confusion Matrix ---\n")
enhanced_cm <- table(Predicted = enhanced_class, Actual = df_complete$S.R)
print(enhanced_cm)

# Calculate metrics
baseline_metrics <- calc_metrics(baseline_cm, levels(df_complete$S.R))
enhanced_metrics <- calc_metrics(enhanced_cm, levels(df_complete$S.R))

# ROC analysis
baseline_roc <- roc(df_complete$S.R, baseline_pred, quiet = TRUE)
enhanced_roc <- roc(df_complete$S.R, enhanced_pred, quiet = TRUE)

# Performance summary
performance_logit <- data.frame(
  Model = c("Baseline", "Enhanced"),
  Accuracy = c(baseline_metrics$accuracy, enhanced_metrics$accuracy),
  Sensitivity = c(baseline_metrics$sensitivity, enhanced_metrics$sensitivity),
  Specificity = c(baseline_metrics$specificity, enhanced_metrics$specificity),
  AUC = c(auc(baseline_roc), auc(enhanced_roc))
)

cat("\n--- Logistic Regression Performance ---\n")
print(performance_logit, digits = 3)

# ------------------------------------------------------------------------------
# 6. RANDOM FOREST
# ------------------------------------------------------------------------------

cat("\n=== 6. RANDOM FOREST ===\n\n")

# --- 6.1 Fit Models ---
rf_baseline <- randomForest(
  baseline_formula,
  data = df_complete,
  ntree = 500,
  importance = TRUE,
  proximity = TRUE
)

rf_enhanced <- randomForest(
  enhanced_formula,
  data = df_complete,
  ntree = 500,
  importance = TRUE,
  proximity = TRUE
)

cat("--- Baseline RF ---\n")
print(rf_baseline)

cat("\n--- Enhanced RF ---\n")
print(rf_enhanced)

# --- 6.2 Prediction Performance ---
pred_baseline_rf <- predict(rf_baseline, type = "prob")[, "R"]
pred_enhanced_rf <- predict(rf_enhanced, type = "prob")[, "R"]

class_baseline_rf <- predict(rf_baseline)
class_enhanced_rf <- predict(rf_enhanced)

# Confusion matrices
cat("\n--- Baseline RF Confusion Matrix ---\n")
cm_baseline_rf <- table(Predicted = class_baseline_rf, Actual = df_complete$S.R)
print(cm_baseline_rf)

cat("\n--- Enhanced RF Confusion Matrix ---\n")
cm_enhanced_rf <- table(Predicted = class_enhanced_rf, Actual = df_complete$S.R)
print(cm_enhanced_rf)

# Calculate metrics
baseline_metrics_rf <- calc_metrics(cm_baseline_rf, levels(df_complete$S.R))
enhanced_metrics_rf <- calc_metrics(cm_enhanced_rf, levels(df_complete$S.R))

# ROC analysis
roc_baseline_rf <- roc(df_complete$S.R, pred_baseline_rf, quiet = TRUE)
roc_enhanced_rf <- roc(df_complete$S.R, pred_enhanced_rf, quiet = TRUE)

# Performance summary
performance_rf <- data.frame(
  Model = c("Baseline", "Enhanced"),
  Accuracy = c(baseline_metrics_rf$accuracy, enhanced_metrics_rf$accuracy),
  Sensitivity = c(baseline_metrics_rf$sensitivity, enhanced_metrics_rf$sensitivity),
  Specificity = c(baseline_metrics_rf$specificity, enhanced_metrics_rf$specificity),
  AUC = c(auc(roc_baseline_rf), auc(roc_enhanced_rf))
)

cat("\n--- Random Forest Performance ---\n")
print(performance_rf, digits = 3)

# --- 6.3 Variable Importance ---
importance_enhanced <- importance(rf_enhanced)
imp_sorted <- sort(importance_enhanced[, "MeanDecreaseAccuracy"], decreasing = TRUE)

cat("\n--- Variable Importance (Mean Decrease Accuracy) ---\n")
print(round(imp_sorted, 3))

# Importance table
importance_table <- data.frame(
  Variable = rownames(importance_enhanced),
  MeanDecreaseAccuracy = importance_enhanced[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance_enhanced[, "MeanDecreaseGini"],
  row.names = NULL
) %>%
  arrange(desc(MeanDecreaseAccuracy))

# ------------------------------------------------------------------------------
# 7. COMPREHENSIVE MODEL COMPARISON
# ------------------------------------------------------------------------------

cat("\n=== 7. COMPREHENSIVE MODEL COMPARISON ===\n\n")

model_comparison_all <- data.frame(
  Method = c("Logistic_Baseline", "Logistic_Enhanced", 
             "RF_Baseline", "RF_Enhanced"),
  Accuracy = c(baseline_metrics$accuracy, enhanced_metrics$accuracy,
               baseline_metrics_rf$accuracy, enhanced_metrics_rf$accuracy),
  Sensitivity = c(baseline_metrics$sensitivity, enhanced_metrics$sensitivity,
                  baseline_metrics_rf$sensitivity, enhanced_metrics_rf$sensitivity),
  Specificity = c(baseline_metrics$specificity, enhanced_metrics$specificity,
                  baseline_metrics_rf$specificity, enhanced_metrics_rf$specificity),
  AUC = c(auc(baseline_roc), auc(enhanced_roc),
          auc(roc_baseline_rf), auc(roc_enhanced_rf)),
  McFadden_R2 = c(model_comparison_logit$McFadden_R2[1], 
                  model_comparison_logit$McFadden_R2[2], NA, NA),
  AIC = c(model_comparison_logit$AIC[1], 
          model_comparison_logit$AIC[2], NA, NA)
)

cat("--- All Models Comparison ---\n")
print(model_comparison_all, digits = 3)

# Performance improvements
cat("\n--- Performance Improvements (Enhanced vs Baseline) ---\n")

cat("\nLogistic Regression:\n")
cat("  Accuracy:    +", round(enhanced_metrics$accuracy - baseline_metrics$accuracy, 3), "\n")
cat("  Sensitivity: +", round(enhanced_metrics$sensitivity - baseline_metrics$sensitivity, 3), "\n")
cat("  Specificity: +", round(enhanced_metrics$specificity - baseline_metrics$specificity, 3), "\n")
cat("  AUC:         +", round(auc(enhanced_roc) - auc(baseline_roc), 3), "\n")

cat("\nRandom Forest:\n")
cat("  Accuracy:    +", round(enhanced_metrics_rf$accuracy - baseline_metrics_rf$accuracy, 3), "\n")
cat("  Sensitivity: +", round(enhanced_metrics_rf$sensitivity - baseline_metrics_rf$sensitivity, 3), "\n")
cat("  Specificity: +", round(enhanced_metrics_rf$specificity - baseline_metrics_rf$specificity, 3), "\n")
cat("  AUC:         +", round(auc(roc_enhanced_rf) - auc(roc_baseline_rf), 3), "\n")

# Best model
best_idx <- which.max(model_comparison_all$AUC)
cat("\n--- Recommended Model ---\n")
cat("Model:", model_comparison_all$Method[best_idx], "\n")
cat("  AUC:", round(model_comparison_all$AUC[best_idx], 3), "\n")
cat("  Accuracy:", round(model_comparison_all$Accuracy[best_idx], 3), "\n")
cat("  Sensitivity:", round(model_comparison_all$Sensitivity[best_idx], 3), "\n")
cat("  Specificity:", round(model_comparison_all$Specificity[best_idx], 3), "\n")

# ------------------------------------------------------------------------------
# 8. EXPORT TABLES
# ------------------------------------------------------------------------------

cat("\n=== 8. EXPORTING TABLES ===\n\n")

write.csv(model_comparison_logit, 
          file = file.path(OUTPUT_DIR, "model_comparison_logistic.csv"), 
          row.names = FALSE)

write.csv(performance_rf, 
          file = file.path(OUTPUT_DIR, "model_comparison_rf.csv"), 
          row.names = FALSE)

write.csv(model_comparison_all, 
          file = file.path(OUTPUT_DIR, "model_comparison_all.csv"), 
          row.names = FALSE)

write.csv(coef_table, 
          file = file.path(OUTPUT_DIR, "coefficient_table.csv"), 
          row.names = FALSE)

write.csv(importance_table, 
          file = file.path(OUTPUT_DIR, "variable_importance.csv"), 
          row.names = FALSE)

cat("Tables exported successfully.\n")

# ------------------------------------------------------------------------------
# 9. VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("\n=== 9. GENERATING FIGURES ===\n\n")

# Store predictions in dataframe
df_complete$pred_logit <- enhanced_pred
df_complete$pred_rf <- pred_enhanced_rf

# --- 9.1 Coefficient Forest Plot ---
coef_data <- tidy(enhanced_logit, conf.int = TRUE) %>%
  filter(term != "(Intercept)")

p_coef <- ggplot(coef_data, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3) +
  labs(title = "Enhanced Logistic Regression Coefficients",
       subtitle = "95% confidence intervals",
       x = "Coefficient Estimate (log-odds)", 
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "fig_logistic_coefficients.png"), 
       plot = p_coef, width = 10, height = 6, dpi = 300)

# --- 9.2 Logistic Prediction Plots ---
p1 <- ggplot(df_complete, aes(x = dna_bind, y = pred_logit)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(title = "By DNA Binding Category", 
       y = "Predicted P(R)", x = NULL) +
  theme_minimal()

p2 <- ggplot(df_complete, aes(x = factor(FS, labels = c("Point", "Frameshift")), 
                              y = pred_logit)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(title = "By Mutation Type", 
       y = "Predicted P(R)", x = NULL) +
  theme_minimal()

df_struct_long <- df_complete %>%
  select(pred_logit, all_of(STRUCTURAL_FEATURES[1:4])) %>%
  pivot_longer(cols = -pred_logit, names_to = "Feature", values_to = "Value")

p3 <- ggplot(df_struct_long, aes(x = Value, y = pred_logit)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 0.8) +
  facet_wrap(~Feature, scales = "free_x", nrow = 1) +
  labs(title = "By Structural Features", 
       x = "ΔΔG Score", y = "Predicted P(R)") +
  theme_minimal()

p_logit_combined <- (p1 | p2) / p3 +
  plot_annotation(title = "Enhanced Logistic Regression Predictions",
                  theme = theme(plot.title = element_text(face = "bold")))

ggsave(file.path(OUTPUT_DIR, "fig_logistic_predictions.png"), 
       plot = p_logit_combined, width = 12, height = 8, dpi = 300)

# --- 9.3 Random Forest Variable Importance ---
p_importance <- ggplot(importance_table, 
                       aes(x = reorder(Variable, MeanDecreaseAccuracy), 
                           y = MeanDecreaseAccuracy)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Variable Importance",
       subtitle = "Enhanced model",
       x = NULL, 
       y = "Mean Decrease in Accuracy") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "fig_rf_importance.png"), 
       plot = p_importance, width = 10, height = 6, dpi = 300)

# --- 9.4 ROC Curve Comparison ---
png(file.path(OUTPUT_DIR, "fig_roc_comparison.png"), 
    width = 8, height = 8, units = "in", res = 300)

plot(baseline_roc, col = "lightblue", lty = 2, lwd = 2,
     main = "ROC Curve Comparison", legacy.axes = TRUE)
plot(enhanced_roc, col = "pink", lty = 2, lwd = 2, add = TRUE)
plot(roc_baseline_rf, col = "blue", lwd = 2, add = TRUE)
plot(roc_enhanced_rf, col = "red", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(
         paste0("RF Enhanced (AUC = ", round(auc(roc_enhanced_rf), 3), ")"),
         paste0("RF Baseline (AUC = ", round(auc(roc_baseline_rf), 3), ")"),
         paste0("Logistic Enhanced (AUC = ", round(auc(enhanced_roc), 3), ")"),
         paste0("Logistic Baseline (AUC = ", round(auc(baseline_roc), 3), ")")
       ),
       col = c("red", "blue", "pink", "lightblue"), 
       lty = c(1, 1, 2, 2), lwd = 2, cex = 0.9)

dev.off()

# --- 9.5 Method Comparison ---
pred_comparison <- df_complete %>%
  select(dna_bind, pred_logit, pred_rf) %>%
  pivot_longer(cols = c(pred_logit, pred_rf), 
               names_to = "Method", 
               values_to = "Prediction") %>%
  mutate(Method = recode(Method, 
                         "pred_logit" = "Logistic Regression",
                         "pred_rf" = "Random Forest"))

p_method <- ggplot(pred_comparison, aes(x = dna_bind, y = Prediction, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Logistic Regression" = "steelblue", 
                               "Random Forest" = "coral")) +
  labs(title = "Prediction Comparison: Logistic Regression vs Random Forest",
       subtitle = "Enhanced models",
       x = "DNA Binding Category", 
       y = "Predicted P(R)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "fig_method_comparison.png"), 
       plot = p_method, width = 10, height = 6, dpi = 300)

cat("Figures exported successfully.\n")

# ------------------------------------------------------------------------------
# 10. SESSION INFORMATION
# ------------------------------------------------------------------------------

cat("\n=== 10. SESSION INFORMATION ===\n\n")
sessionInfo()
