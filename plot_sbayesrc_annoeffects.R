#!/usr/bin/env Rscript

# Script to plot SBayesRC AnnoEffects1 results
# Creates facet plots for ATAC and RNA annotations

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Read the data file
data_file <- "/Users/uqjzeng1/Work/Projects/StagePRS/res/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt_imp.ma.imputed.ma_cepo_dev_long_100kb_peak_0kb_Function_anno_SBayesRC.parSetRes"
cat("Reading data from:", data_file, "\n")

# Read the file (skip header line)
data <- read.table(data_file, header = TRUE, stringsAsFactors = FALSE, sep = "", 
                   comment.char = "", quote = "", fill = TRUE)

# Filter for AnnoEffects1
annoeffects1 <- data %>% 
  filter(Parameter == "AnnoEffects1")

# Extract Intercept value
intercept <- annoeffects1 %>% 
  filter(Annotation == "Intercept") %>% 
  pull(Mean)

if (length(intercept) == 0) {
  stop("Intercept not found in AnnoEffects1")
}
cat("Intercept value:", intercept, "\n")

# Define stages
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

# Filter for Celltype_Stage annotations
celltype_stage <- annoeffects1 %>% 
  filter(Annotation != "Intercept") %>%
  # Check if annotation contains any of the stages
  filter(str_detect(Annotation, paste(stages, collapse = "|")))

# Function to extract cell type and stage from annotation name
parse_annotation <- function(anno) {
  # Remove "_peak" suffix if present
  has_peak <- str_detect(anno, "_peak$")
  anno_clean <- str_replace(anno, "_peak$", "")
  
  # Try to match pattern: Celltype_Stage
  # Stages are: Fetal, Neonatal, Infancy, Childhood, Adolescence, Adult
  # Match from the end to handle cell types with underscores
  for (stage in stages) {
    pattern <- paste0("_", stage, "$")
    if (str_detect(anno_clean, pattern)) {
      # Extract cell type by removing the stage suffix
      celltype <- str_replace(anno_clean, pattern, "")
      return(list(celltype = celltype, stage = stage, has_peak = has_peak))
    }
  }
  
  # If no match, return NA
  return(list(celltype = NA, stage = NA, has_peak = has_peak))
}

# Parse annotations
parsed <- celltype_stage %>%
  rowwise() %>%
  mutate(
    parsed_info = list(parse_annotation(Annotation)),
    celltype = parsed_info$celltype,
    stage = parsed_info$stage,
    has_peak = parsed_info$has_peak,
    annotation_type = if_else(has_peak, "ATAC", "RNA")
  ) %>%
  ungroup() %>%
  filter(!is.na(celltype), !is.na(stage)) %>%
  select(Annotation, Mean, SD, celltype, stage, annotation_type)

# Calculate pnorm(Intercept + Celltype_stage) - pnorm(Intercept)
parsed <- parsed %>%
  mutate(
    intercept_plus_effect = intercept + Mean,
    prob_with = pnorm(intercept_plus_effect),
    prob_intercept = pnorm(intercept),
    prob_diff = prob_with - prob_intercept
  )

# Order stages
parsed$stage <- factor(parsed$stage, levels = stages)

# Get unique cell types
celltypes <- sort(unique(parsed$celltype))
cat("Found", length(celltypes), "cell types\n")
cat("Cell types:", paste(celltypes, collapse = ", "), "\n")

# Create plots for each annotation type (ATAC and RNA)
for (anno_type in c("ATAC", "RNA")) {
  cat("\nCreating plot for", anno_type, "annotations\n")
  
  plot_data <- parsed %>% 
    filter(annotation_type == anno_type)
  
  if (nrow(plot_data) == 0) {
    cat("No data found for", anno_type, "annotations\n")
    next
  }
  
  # Create the plot with rows=stages, columns=celltypes
  p <- ggplot(plot_data, aes(x = celltype, y = stage, fill = prob_diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red",
      midpoint = 0,
      name = expression(paste(Delta, " Probability")),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    labs(
      x = "Cell Type",
      y = "Developmental Stage",
      title = paste("SBayesRC AnnoEffects1:", anno_type, "Annotations"),
      subtitle = expression(paste("Values: ", Phi(Intercept + Effect) - Phi(Intercept)))
    ) +
    coord_fixed(ratio = 1)
  
  # Save the plot
  output_file <- paste0("/Users/uqjzeng1/Work/Projects/StagePRS/res/AnnoEffects1_", anno_type, "_plot.pdf")
  ggsave(output_file, plot = p, width = max(10, length(celltypes) * 0.8), 
         height = 8, units = "in")
  cat("Saved plot to:", output_file, "\n")
  
  # Also save as PNG
  output_file_png <- paste0("/Users/uqjzeng1/Work/Projects/StagePRS/res/AnnoEffects1_", anno_type, "_plot.png")
  ggsave(output_file_png, plot = p, width = max(10, length(celltypes) * 0.8), 
         height = 8, units = "in", dpi = 300)
  cat("Saved plot to:", output_file_png, "\n")
}

# Create a combined plot with both ATAC and RNA as separate facets
# Rows = stages, Columns = cell types, Facets = annotation type
cat("\nCreating combined plot with ATAC and RNA\n")

p_combined <- ggplot(parsed, aes(x = celltype, y = stage, fill = prob_diff)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = expression(paste(Delta, " Probability")),
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  facet_grid(annotation_type ~ ., scales = "free", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(
    x = "Cell Type",
    y = "Developmental Stage",
    title = "SBayesRC AnnoEffects1: ATAC and RNA Annotations",
    subtitle = expression(paste("Values: ", Phi(Intercept + Effect) - Phi(Intercept)))
  ) +
  coord_fixed(ratio = 1)

# Save combined plot
output_file_combined <- "/Users/uqjzeng1/Work/Projects/StagePRS/res/AnnoEffects1_combined_plot.pdf"
ggsave(output_file_combined, plot = p_combined, width = max(12, length(celltypes) * 0.8), 
       height = 14, units = "in")
cat("Saved combined plot to:", output_file_combined, "\n")

output_file_combined_png <- "/Users/uqjzeng1/Work/Projects/StagePRS/res/AnnoEffects1_combined_plot.png"
ggsave(output_file_combined_png, plot = p_combined, width = max(12, length(celltypes) * 0.8), 
       height = 14, units = "in", dpi = 300)
cat("Saved combined plot to:", output_file_combined_png, "\n")

cat("\nDone!\n")

