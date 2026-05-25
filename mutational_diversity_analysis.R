# ============================================================================
# MUTATIONAL DIVERSITY ANALYSIS ACROSS ORFs BY ANTI-HBC STATUS
# ============================================================================
# This script analyzes genetic diversity/variation in HBV mutations
# Measures include: Shannon Diversity, Simpson Diversity, Richness, Evenness
# Stratified by anti-HBc status (Negative, Positive, Unknown)
# ============================================================================

library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)  # For diversity indices
library(gridExtra)

# ============================================================================
# 1. READ DATA
# ============================================================================

df <- read_excel("ws_95.xlsx")

print("Data loaded successfully!")
print(head(df))

# ============================================================================
# 2. DEFINE ORF COLUMNS AND STATUS LABELS
# ============================================================================

orf_columns <- list(
  sp = "mutations_sp",
  core = "mutations_core",
  shb = "mutations_shb",
  mhb = "mutations_mhb",
  lhb = "mutations_lhb",
  prec = "mutations_prec",
  rt = "mutations_rt"
)

status_labels <- c("0" = "Unknown", "1" = "Negative", "2" = "Positive")

# ============================================================================
# 3. FUNCTION TO PARSE MUTATIONS
# ============================================================================

parse_mutations <- function(mutation_string) {
  if (is.na(mutation_string) || mutation_string == "" || trimws(mutation_string) == "") {
    return(character(0))
  }
  
  mutations <- str_trim(unlist(strsplit(as.character(mutation_string), ",")))
  mutations <- mutations[mutations != ""]
  return(mutations)
}

# ============================================================================
# 4. CALCULATE DIVERSITY INDICES FOR EACH ORF BY ANTI-HBC STATUS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("MUTATIONAL DIVERSITY ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

diversity_indices <- data.frame(
  ORF = character(),
  Anti_HBc_Status = character(),
  Anti_HBc_Code = numeric(),
  N_Samples = numeric(),
  N_Mutations = numeric(),
  N_Unique_Mutations = numeric(),
  Shannon_Diversity = numeric(),
  Simpson_Diversity = numeric(),
  Richness = numeric(),
  Pielou_Evenness = numeric(),
  Mean_Mutations_Per_Sample = numeric(),
  SD_Mutations_Per_Sample = numeric(),
  stringsAsFactors = FALSE
)

for (orf_name in names(orf_columns)) {
  col_name <- orf_columns[[orf_name]]
  
  for (status_code in c(0, 1, 2)) {
    status_label <- status_labels[[as.character(status_code)]]
    status_data <- df %>% filter(anti_hbc_status == status_code)
    
    if (nrow(status_data) > 0) {
      # Get all mutations for this ORF and status
      all_mutations <- unlist(lapply(status_data[[col_name]], parse_mutations))
      unique_mutations <- unique(all_mutations)
      
      n_samples <- nrow(status_data)
      n_total_mutations <- length(all_mutations)
      n_unique <- length(unique_mutations)
      
      # Calculate mutation frequencies for diversity indices
      if (n_total_mutations > 0) {
        mutation_counts <- as.numeric(table(all_mutations))
        
        # Shannon Diversity (higher = more diverse)
        shannon <- diversity(mutation_counts, index = "shannon")
        
        # Simpson Diversity (1 - Simpson Index, higher = more diverse)
        simpson <- diversity(mutation_counts, index = "simpson")
        
        # Richness (number of unique mutations)
        richness <- n_unique
        
        # Pielou's Evenness (Shannon / ln(richness))
        # Ranges from 0 to 1, higher = more even distribution
        if (richness > 1) {
          pielou <- shannon / log(richness)
        } else {
          pielou <- NA
        }
        
        # Mean and SD mutations per sample
        mutations_per_sample <- sapply(status_data[[col_name]], function(x) length(parse_mutations(x)))
        mean_per_sample <- mean(mutations_per_sample, na.rm = TRUE)
        sd_per_sample <- sd(mutations_per_sample, na.rm = TRUE)
        
      } else {
        shannon <- 0
        simpson <- 0
        richness <- 0
        pielou <- 0
        mean_per_sample <- 0
        sd_per_sample <- 0
      }
      
      diversity_indices <- rbind(diversity_indices, data.frame(
        ORF = toupper(orf_name),
        Anti_HBc_Status = status_label,
        Anti_HBc_Code = status_code,
        N_Samples = n_samples,
        N_Mutations = n_total_mutations,
        N_Unique_Mutations = n_unique,
        Shannon_Diversity = round(shannon, 4),
        Simpson_Diversity = round(simpson, 4),
        Richness = richness,
        Pielou_Evenness = round(pielou, 4),
        Mean_Mutations_Per_Sample = round(mean_per_sample, 2),
        SD_Mutations_Per_Sample = round(sd_per_sample, 2),
        stringsAsFactors = FALSE
      ))
    }
  }
}

print("Diversity Indices Summary:")
print(diversity_indices)
write.csv(diversity_indices, "div_01_diversity_indices.csv", row.names = FALSE)
cat("\n✓ Saved: div_01_diversity_indices.csv\n")

# ============================================================================
# 5. GRAPH 1: SHANNON DIVERSITY INDEX BY ORF AND ANTI-HBC STATUS
# ============================================================================

cat("\nCreating diversity visualizations...\n")

p_shannon <- ggplot(diversity_indices, aes(x = reorder(ORF, -Shannon_Diversity), y = Shannon_Diversity, fill = Anti_HBc_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  geom_text(aes(label = Shannon_Diversity), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  labs(
    title = "Shannon Diversity Index by ORF and Anti-HBc Status",
    subtitle = "Measures genetic diversity (higher = more diverse mutation patterns)",
    x = "Open Reading Frame",
    y = "Shannon Diversity Index",
    fill = "Anti-HBc Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p_shannon)
ggsave("div_graph_01_shannon_diversity.png", plot = p_shannon, width = 12, height = 7, dpi = 300)
cat("✓ Saved: div_graph_01_shannon_diversity.png\n")

# ============================================================================
# 6. GRAPH 2: SIMPSON DIVERSITY INDEX BY ORF AND ANTI-HBC STATUS
# ============================================================================

p_simpson <- ggplot(diversity_indices, aes(x = reorder(ORF, -Simpson_Diversity), y = Simpson_Diversity, fill = Anti_HBc_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  geom_text(aes(label = Simpson_Diversity), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  labs(
    title = "Simpson Diversity Index by ORF and Anti-HBc Status",
    subtitle = "Alternative diversity measure (higher = more diverse, less dominated by common mutations)",
    x = "Open Reading Frame",
    y = "Simpson Diversity Index",
    fill = "Anti-HBc Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p_simpson)
ggsave("div_graph_02_simpson_diversity.png", plot = p_simpson, width = 12, height = 7, dpi = 300)
cat("✓ Saved: div_graph_02_simpson_diversity.png\n")

# ============================================================================
# 7. GRAPH 3: RICHNESS (NUMBER OF UNIQUE MUTATIONS) BY ORF AND STATUS
# ============================================================================

p_richness <- ggplot(diversity_indices, aes(x = reorder(ORF, -Richness), y = Richness, fill = Anti_HBc_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  geom_text(aes(label = Richness), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  labs(
    title = "Mutation Richness (Number of Unique Mutations) by ORF and Anti-HBc Status",
    subtitle = "Count of distinct mutations (higher = more genetic variants)",
    x = "Open Reading Frame",
    y = "Number of Unique Mutations (Richness)",
    fill = "Anti-HBc Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p_richness)
ggsave("div_graph_03_richness.png", plot = p_richness, width = 12, height = 7, dpi = 300)
cat("✓ Saved: div_graph_03_richness.png\n")

# ============================================================================
# 8. GRAPH 4: PIELOU'S EVENNESS INDEX BY ORF AND STATUS
# ============================================================================

p_evenness <- ggplot(diversity_indices, aes(x = reorder(ORF, -Pielou_Evenness), y = Pielou_Evenness, fill = Anti_HBc_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  geom_text(aes(label = Pielou_Evenness), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  ylim(0, 1) +
  labs(
    title = "Pielou's Evenness Index by ORF and Anti-HBc Status",
    subtitle = "Measures how evenly mutations are distributed (0-1, higher = more even distribution)",
    x = "Open Reading Frame",
    y = "Pielou's Evenness Index",
    fill = "Anti-HBc Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p_evenness)
ggsave("div_graph_04_evenness.png", plot = p_evenness, width = 12, height = 7, dpi = 300)
cat("✓ Saved: div_graph_04_evenness.png\n")

# ============================================================================
# 9. GRAPH 5: SHANNON DIVERSITY - FACETED BY ORF
# ============================================================================

p_shannon_facet <- ggplot(diversity_indices, aes(x = Anti_HBc_Status, y = Shannon_Diversity, fill = Anti_HBc_Status)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = Shannon_Diversity), vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ORF, ncol = 4) +
  scale_fill_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  labs(
    title = "Shannon Diversity Index by Anti-HBc Status (Separate Panel per ORF)",
    subtitle = "Compare genetic diversity across anti-HBc groups for each ORF",
    x = "Anti-HBc Status",
    y = "Shannon Diversity Index",
    fill = "Anti-HBc Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p_shannon_facet)
ggsave("div_graph_05_shannon_faceted.png", plot = p_shannon_facet, width = 14, height = 10, dpi = 300)
cat("✓ Saved: div_graph_05_shannon_faceted.png\n")

# ============================================================================
# 10. GRAPH 6: COMBINED DIVERSITY METRICS - LINE PLOT
# ============================================================================

# Prepare data for line plot
diversity_long <- diversity_indices %>%
  select(ORF, Anti_HBc_Status, Shannon_Diversity, Simpson_Diversity) %>%
  pivot_longer(
    cols = c(Shannon_Diversity, Simpson_Diversity),
    names_to = "Metric",
    values_to = "Value"
  )

p_combined <- ggplot(diversity_long, aes(x = ORF, y = Value, color = Anti_HBc_Status, linetype = Metric, group = interaction(Anti_HBc_Status, Metric))) +
  geom_line(size = 1.2) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Negative" = "#2ecc71", "Positive" = "#e74c3c", "Unknown" = "#95a5a6")) +
  scale_linetype_manual(values = c("Shannon_Diversity" = "solid", "Simpson_Diversity" = "dashed")) +
  labs(
    title = "Diversity Metrics Comparison Across ORFs by Anti-HBc Status",
    subtitle = "Shannon (solid) vs Simpson (dashed) diversity indices",
    x = "Open Reading Frame",
    y = "Diversity Index Value",
    color = "Anti-HBc Status",
    linetype = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_line(color = "gray90")
  )

print(p_combined)
ggsave("div_graph_06_combined_diversity_lines.png", plot = p_combined, width = 13, height = 7, dpi = 300)
cat("✓ Saved: div_graph_06_combined_diversity_lines.png\n")

# ============================================================================
# 11. STATISTICAL COMPARISON: ANOVA FOR DIVERSITY DIFFERENCES
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("STATISTICAL ANALYSIS: Testing for Significant Differences in Diversity\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

statistical_results <- data.frame(
  ORF = character(),
  Metric = character(),
  F_statistic = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

for (orf_name in names(orf_columns)) {
  orf_data <- diversity_indices %>% filter(ORF == toupper(orf_name))
  
  if (nrow(orf_data) > 0 && length(unique(orf_data$Anti_HBc_Status)) > 1) {
    
    # ANOVA for Shannon Diversity
    aov_shannon <- aov(Shannon_Diversity ~ Anti_HBc_Status, data = orf_data)
    summary_shannon <- summary(aov_shannon)
    f_shannon <- summary_shannon[[1]]$`F value`[1]
    p_shannon <- summary_shannon[[1]]$`Pr(>F)`[1]
    sig_shannon <- ifelse(p_shannon < 0.05, "Yes", "No")
    
    statistical_results <- rbind(statistical_results, data.frame(
      ORF = toupper(orf_name),
      Metric = "Shannon Diversity",
      F_statistic = round(f_shannon, 4),
      p_value = round(p_shannon, 4),
      Significant = sig_shannon,
      stringsAsFactors = FALSE
    ))
    
    # ANOVA for Simpson Diversity
    aov_simpson <- aov(Simpson_Diversity ~ Anti_HBc_Status, data = orf_data)
    summary_simpson <- summary(aov_simpson)
    f_simpson <- summary_simpson[[1]]$`F value`[1]
    p_simpson <- summary_simpson[[1]]$`Pr(>F)`[1]
    sig_simpson <- ifelse(p_simpson < 0.05, "Yes", "No")
    
    statistical_results <- rbind(statistical_results, data.frame(
      ORF = toupper(orf_name),
      Metric = "Simpson Diversity",
      F_statistic = round(f_simpson, 4),
      p_value = round(p_simpson, 4),
      Significant = sig_simpson,
      stringsAsFactors = FALSE
    ))
    
    # ANOVA for Richness
    aov_richness <- aov(Richness ~ Anti_HBc_Status, data = orf_data)
    summary_richness <- summary(aov_richness)
    f_richness <- summary_richness[[1]]$`F value`[1]
    p_richness <- summary_richness[[1]]$`Pr(>F)`[1]
    sig_richness <- ifelse(p_richness < 0.05, "Yes", "No")
    
    statistical_results <- rbind(statistical_results, data.frame(
      ORF = toupper(orf_name),
      Metric = "Richness",
      F_statistic = round(f_richness, 4),
      p_value = round(p_richness, 4),
      Significant = sig_richness,
      stringsAsFactors = FALSE
    ))
  }
}

print("Statistical Test Results (ANOVA):")
print(statistical_results)
write.csv(statistical_results, "div_02_statistical_tests.csv", row.names = FALSE)
cat("\n✓ Saved: div_02_statistical_tests.csv\n")

# ============================================================================
# 12. CREATE SUMMARY TABLE FOR INTERPRETATION
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("INTERPRETATION GUIDE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

cat("DIVERSITY METRICS EXPLAINED:\n\n")

cat("1. SHANNON DIVERSITY INDEX (0 to ~4)\n")
cat("   - Measures both richness (number of types) and evenness (how balanced they are)\n")
cat("   - Higher values = more diverse mutation profiles\n")
cat("   - Formula: H = -Σ(pi * ln(pi)) where pi = proportion of each mutation\n")
cat("   - Interpretation: Positive > Negative suggests more diverse mutations in infected\n\n")

cat("2. SIMPSON DIVERSITY INDEX (0 to 1)\n")
cat("   - Probability that two randomly selected mutations are different\n")
cat("   - Higher values = more diverse (less dominated by single mutation)\n")
cat("   - Interpretation: Positive > Negative = more genetic heterogeneity in infected\n\n")

cat("3. RICHNESS (count of unique mutations)\n")
cat("   - Simple count of different mutations present\n")
cat("   - Higher = more different mutation types\n")
cat("   - Interpretation: Raw measure of genetic variants\n\n")

cat("4. PIELOU'S EVENNESS (0 to 1)\n")
cat("   - How uniformly mutations are distributed (range 0-1)\n")
cat("   - Value of 1 = all mutations equally common\n")
cat("   - Value close to 0 = dominated by few mutations\n")
cat("   - Interpretation: High evenness = balanced mutation distribution\n\n")

# ============================================================================
# 13. FINAL SUMMARY
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

cat("SUMMARY:\n")
cat(paste("Total samples analyzed:", nrow(df), "\n"))
cat(paste("  - Negative:", nrow(df %>% filter(anti_hbc_status == 1)), "\n"))
cat(paste("  - Positive:", nrow(df %>% filter(anti_hbc_status == 2)), "\n"))
cat(paste("  - Unknown:", nrow(df %>% filter(anti_hbc_status == 0)), "\n\n"))

cat("DATA FILES CREATED:\n")
cat("  1. div_01_diversity_indices.csv - All diversity metrics by ORF and status\n")
cat("  2. div_02_statistical_tests.csv - ANOVA results testing differences\n\n")

cat("GRAPHS CREATED:\n")
cat("  1. div_graph_01_shannon_diversity.png - Shannon index comparison\n")
cat("  2. div_graph_02_simpson_diversity.png - Simpson index comparison\n")
cat("  3. div_graph_03_richness.png - Richness (unique mutations count)\n")
cat("  4. div_graph_04_evenness.png - Pielou's evenness index\n")
cat("  5. div_graph_05_shannon_faceted.png - Shannon by status per ORF\n")
cat("  6. div_graph_06_combined_diversity_lines.png - Shannon vs Simpson trends\n\n")

cat("✓ Analysis complete! Review the CSV files and graphs for insights on\n")
cat("  mutational diversity patterns across ORFs and anti-HBc groups.\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
