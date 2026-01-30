# ==============================================================================
# DATA FILTER TO HIGH-CONFIDENCE mmpR5 MUTANTS
# STATISTICAL ASSESSMENT OF LOW-CONFIDENCE mmpR5 MUTANTS
# ==============================================================================
#
# Title:        Binomial Testing for Phenotype Classification of Low-Confidence
#               mmpR5 Mutants
# Author:       Xuan Lu
# Date:         01/03/2026
#
# Description:  This script identifies mmpR5 mutants with conflicting resistance
#               phenotypes (low-confidence) and applies exact binomial tests to 
#               assess whether each mutant shows a statistically significant 
#               tendency toward resistance (R) or susceptibility (S).
#
# Input:        GscoreMMPR5mt.xlsx (sheet: "Gscore")
# Output:       - data_highconf_singleocc.csv (filtered dataset for modeling)
#               - low_confidence_summary.csv (statistical results)
#               - fig_forest_low_confidence.png (forest plot visualization)
#
# Methods:      Exact binomial test with Benjamini-Hochberg FDR correction
#               Null hypothesis: P(Resistant) = 0.5
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(readxl)
library(pROC)        # Note: masks ggplot2::margin()
library(ggplot2)

# Set seed for reproducibility
set.seed(2025)

# Define file paths (modify as needed)
INPUT_FILE  <- "GscoreMMPR5mt.xlsx"
OUTPUT_DIR  <- "."

# ------------------------------------------------------------------------------
# 1. DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------

# Load data
df <- read_excel(INPUT_FILE, sheet = "Gscore")
colnames(df) <- make.names(colnames(df))

# Recode intermediate (I) phenotype as resistant (R)
# Rationale: Conservative approach treating intermediate MIC as resistant
df$S.R <- ifelse(df$S.I.R == "I", "R", df$S.I.R)
df$S.R <- factor(df$S.R, levels = c("S", "R"), ordered = TRUE)
df <- df %>% select(-"S.I.R")

# Remove samples with missing phenotype
df <- df[!is.na(df$S.R), ]

# Prepare mutation factor with all observed levels
df$mmpR5 <- factor(df$mmpR5, levels = unique(df$mmpR5))

cat("Data loaded:", nrow(df), "samples with", n_distinct(df$mmpR5), "unique mutants\n")

# ------------------------------------------------------------------------------
# 2. MUTATION DOMAIN ANNOTATION
# ------------------------------------------------------------------------------

# Define mutations by structural domain location
# Based on mmpR5 protein structure analysis

dna_binding <- c(

  # Frameshift insertions in DNA-binding region
  '132insG', '132insGG', '138insG', '140insGG', '141insG', '192insG', '211insG', '234insG',
  # Point mutations in DNA-binding domain
 'A59V', 'A62T', 'A62V', 'A84V', 'C46G', 'C46R', 'C46Y', 'E55D', 'G41V', 'G65R', 
  'I67L', 'I67M', 'I67S', 'I67V', 'I80M', 'I80S', 'L40F', 'L40M', 'L40S', 'L40V', 
  'L74M', 'L74V', 'L83F', 'L83P', 'M73I', 'N70D', 'P48H', 'P48L', 'Q51K', 'Q51R', 
  'Q76E', 'R50Q', 'R72W', 'S52F', 'S52P', 'S53P', 'S63G', 'S63R', 'S68G', 'S68I', 
  'S68N', 'S68R', 'T58P', 'V85A', 'W42C', 'W42R'
)

dna_non_binding <- c(
  # Frameshift insertions outside DNA-binding region
  '107insG', '274insG', '28insG', '324insG', '418insG', '464insG', '465insG', '492insGG',
  # Point mutations outside DNA-binding domain
  'A101T', 'A102V', 'A118T', 'A36T', 'A36V', 'A86T', 'A86V', 'A99P', 'D15G', 'D5G', 
  'D5N', 'D88G', 'E104G', 'E113K', 'E147D', 'E147G', 'E21D', 'E21K', 'F100Y', 'F19S', 
  'F79S', 'F93L', 'F93S', 'G103S', 'G121R', 'G162E', 'G25C', 'G87A', 'G87R', 'I108S', 
  'I108T', 'I108V', 'L114P', 'L117R', 'L125R', 'L136P', 'L142R', 'L154R', 'M10I', 
  'M139I', 'M146T', 'M17V', 'M23V', 'N4T', 'N98D', 'P14L', 'P14S', 'Q22R', 'R105C', 
  'R107C', 'R109P', 'R109Q', 'R123K', 'R134G', 'R135W', 'R137P', 'R156Q', 'R90C', 
  'R90L', 'R94Q', 'R94W', 'R96G', 'R96L', 'R96Q', 'R96T', 'R96W', 'S158R', 'S2I', 
  'S2R', 'S31R', 'T33S', 'V20A', 'V20G', 'Y145H', 'Y26C', 'L32S', 'Y92D', 'M139T', 
  'D8G', 'A101E', 'L39S'
)

# Annotate mutations
df <- df %>%
  mutate(
    FS = ifelse(grepl("ins", mmpR5), 1, 0),
    dna_bind = case_when(
      mmpR5 %in% dna_binding ~ "DNA_binding",
      mmpR5 %in% dna_non_binding ~ "DNA_non_binding",
      TRUE ~ NA_character_
    ),
    dna_bind = factor(dna_bind, levels = c("DNA_non_binding", "DNA_binding"))
  )

# ------------------------------------------------------------------------------
# 3. CONFIDENCE LEVEL CLASSIFICATION
# ------------------------------------------------------------------------------

# Classify mutants based on phenotype consistency across samples
# - High confidence: Multiple samples, all same phenotype
# - Low confidence:  Multiple samples, mixed phenotypes (requires statistical assessment)
# - Single occurrence: Only one sample observed

df_labeled <- df %>%
  group_by(mmpR5) %>%
  mutate(
    n_occurrences = n(),
    n_unique_phenotypes = n_distinct(S.R),
    confidence = case_when(
      n_occurrences == 1 ~ "single_occurrence",
      n_unique_phenotypes == 1 ~ "high_confidence",
      n_unique_phenotypes > 1 ~ "low_confidence"
    )
  ) %>%
  ungroup()

# Summary of confidence levels
confidence_summary <- df_labeled %>%
  group_by(confidence) %>%
  summarise(
    n_samples = n(),
    n_unique_mutants = n_distinct(mmpR5),
    .groups = "drop"
  )

cat("\n--- Confidence Level Summary ---\n")
print(confidence_summary)

# Phenotype distribution by confidence level
phenotype_by_conf <- df_labeled %>%
  group_by(confidence, S.R) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = S.R, values_from = count, values_fill = 0)

cat("\n--- Phenotype Distribution by Confidence Level ---\n")
print(phenotype_by_conf)

# ------------------------------------------------------------------------------
# 4. EXPORT HIGH-CONFIDENCE DATASET FOR MODELING
# ------------------------------------------------------------------------------

# Filter to high-confidence and single-occurrence samples
df_model <- df_labeled %>%
  filter(confidence %in% c("high_confidence", "single_occurrence"))

cat("\n--- Modeling Dataset ---\n")
cat("Samples retained:", nrow(df_model), "\n")
cat("Unique mutants:", n_distinct(df_model$mmpR5), "\n")

# Export
write.csv(df_model, file = file.path(OUTPUT_DIR, "data_highconf_singleocc.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. STATISTICAL ASSESSMENT OF LOW-CONFIDENCE MUTANTS
# ------------------------------------------------------------------------------

# Calculate overall resistance rate as reference
overall_p_R <- mean(df_labeled$S.R == "R")
cat("\n--- Statistical Assessment ---\n")
cat("Overall resistance rate:", round(overall_p_R, 3), "\n")

# Summarize low-confidence mutants
low_conf_stats <- df_labeled %>%
  filter(confidence == "low_confidence") %>%
  group_by(mmpR5) %>%
  summarise(
    n_R = sum(S.R == "R"),
    n_S = sum(S.R == "S"),
    n_total = n(),
    prop_R = n_R / n_total,
    .groups = "drop"
  )

# Apply exact binomial test to each mutant
# H0: P(Resistant) = 0.5
low_conf_stats <- low_conf_stats %>%
  rowwise() %>%
  mutate(
    p_value = binom.test(n_R, n_total, p = 0.5)$p.value,
    ci_lower = binom.test(n_R, n_total)$conf.int[1],
    ci_upper = binom.test(n_R, n_total)$conf.int[2]
  ) %>%
  ungroup()

# Apply Benjamini-Hochberg FDR correction
low_conf_stats <- low_conf_stats %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Classify tendency based on adjusted p-values
low_conf_stats <- low_conf_stats %>%
  mutate(
    classification = case_when(
      p_adj < 0.05 & prop_R > 0.5 ~ "Likely_R",
      p_adj < 0.05 & prop_R < 0.5 ~ "Likely_S",
      TRUE ~ "Uncertain"
    )
  )

# Format and display results
low_conf_summary <- low_conf_stats %>%
  arrange(p_adj) %>%
  select(
    Mutant = mmpR5,
    n_R, n_S, n_total,
    Proportion_R = prop_R,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    p_value,
    p_adj,
    Classification = classification
  )

cat("\n--- Low-Confidence Mutant Assessment ---\n")
print(low_conf_summary, n = Inf)

# Classification summary
cat("\n--- Classification Summary ---\n")
print(table(low_conf_stats$classification))

# Export results
write.csv(low_conf_summary, file = file.path(OUTPUT_DIR, "low_confidence_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. VISUALIZATION: FOREST PLOT
# ------------------------------------------------------------------------------

# Prepare plot data
plot_data <- low_conf_stats %>%
  mutate(
    mmpR5 = fct_reorder(mmpR5, prop_R),
    label = paste0(n_R, "/", n_total)
  )

# Create forest plot
p_forest <- ggplot(plot_data, aes(x = prop_R, y = mmpR5, color = classification)) +
  
  # Reference lines
  geom_vline(xintercept = 0.5, linetype = "solid", color = "red", linewidth = 0.5) +
  geom_vline(xintercept = overall_p_R, linetype = "dotted", color = "steelblue", linewidth = 0.8) +
  
  # Reference line label
annotate("text", x = overall_p_R, y = Inf, 
           label = paste0("Overall (", round(overall_p_R * 100, 1), "%)"), 
           vjust = -0.5, hjust = 0.5, size = 3, color = "steelblue") +
  
  # Confidence intervals
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, linewidth = 0.7) +
  
  # Point estimates
  geom_point(aes(size = n_total), alpha = 0.8) +
  
  # Sample count labels
  geom_text(aes(x = ci_upper + 0.05, label = label), size = 3, hjust = 0, color = "gray30") +
  
  # Scales
  scale_color_manual(
    values = c("Likely_R" = "#D55E00", "Likely_S" = "#0072B2", "Uncertain" = "gray50"),
    name = "Classification"
  ) +
  scale_size_continuous(range = c(2, 5), name = "Sample Size") +
  scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(clip = "off") +
  
  # Labels
  labs(
    title = "Resistance Probability for Low-Confidence mmpR5 Mutants",
    x = "Proportion Resistant (95% CI)",
    y = NULL,
    caption = "Labels: n resistant / n total. Red line: 50% threshold."
  ) +
  
  # Theme
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray40", size = 9),
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10)
  )

print(p_forest)

# Save figure
fig_height <- max(4, nrow(plot_data) * 0.35)
ggsave(
  filename = file.path(OUTPUT_DIR, "fig_forest_low_confidence.png"),
  plot = p_forest,
  width = 8,
  height = fig_height,
  dpi = 300
)
