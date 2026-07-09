# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Read data (header = FALSE because the first row contains actual data values)
dat <- read.csv("dat.csv", header = FALSE)
colnames(dat) <- c("gene", "HH14_1", "HH14_2", "HH14_3", "HH14_4", "HH19_1", "HH19_2", "HH19_3", "HH19_4")

# 2. Reshape from wide to long format
data_long <- dat %>%
  pivot_longer(
    cols = -gene, 
    names_to = "replicate", 
    values_to = "Normalized_Expression"
  ) %>%
  mutate(
    # Isolate "HH14" or "HH19" from the column labels to define the stage factor
    stage = factor(gsub("_.*", "", replicate), levels = c("HH14", "HH19"))
  )

# Define custom colors for the stages
stage_colors <- c("HH14" = "#61CBF4", "HH19" = "#215F9A")

#---------------------------------------------------------
# PANEL 1: HH14-Upregulated Genes
#---------------------------------------------------------

hh14_genes <- data_long %>%
  filter(gene %in% c("ALDH1A2", "RDH10", "LIN28A", "TRIM71", "OSR1", "MEOX1", "LFNG", "HOXB8")) %>%
  mutate(gene = factor(gene, levels = c(
    "ALDH1A2", "RDH10", "LIN28A", "TRIM71", # Row 1
    "OSR1", "MEOX1", "LFNG", "HOXB8"        # Row 2
  )))

plot_hh14 <- ggplot(hh14_genes, aes(x = stage, y = Normalized_Expression, fill = stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black", alpha = 0.8) +
  geom_jitter(aes(color = stage), position = position_jitter(0.15), size = 1.8, alpha = 0.9) +
  # Force x-axes onto all individual panels using axes = "all"
  facet_wrap(~gene, nrow = 2, scales = "free_y", axes = "all") + 
  scale_fill_manual(values = stage_colors) +
  scale_color_manual(values = c("HH14" = "#3fa5cc", "HH19" = "#133f6b")) + 
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    legend.position = "none", 
    strip.text = element_text(size = 11, face = "italic"), 
    strip.background = element_blank(),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(y = "Normalized Expression")

ggsave("HH14_upregulated.png", plot = plot_hh14, width = 20, height = 11, dpi = 300, units = "cm")

#---------------------------------------------------------
# PANEL 2: HH19-Upregulated Genes
#---------------------------------------------------------

hh19_genes <- data_long %>%
  filter(gene %in% c("FGF10", "WNT5A", "GREM1", "DUSP6", "SP8", "GLI3", "RSPO2", "HOXD11")) %>%
  mutate(gene = factor(gene, levels = c(
    "FGF10", "WNT5A", "GREM1", "DUSP6",   # Row 1
    "SP8", "GLI3", "RSPO2", "HOXD11"       # Row 2
  )))

plot_hh19 <- ggplot(hh19_genes, aes(x = stage, y = Normalized_Expression, fill = stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black", alpha = 0.8) +
  geom_jitter(aes(color = stage), position = position_jitter(0.15), size = 1.8, alpha = 0.9) +
  # Force x-axes onto all individual panels using axes = "all"
  facet_wrap(~gene, nrow = 2, scales = "free_y", axes = "all") + 
  scale_fill_manual(values = stage_colors) +
  scale_color_manual(values = c("HH14" = "#3fa5cc", "HH19" = "#133f6b")) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "italic"),
    strip.background = element_blank(),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(y = "Normalized Expression")

ggsave("HH19_upregulated.png", plot = plot_hh19, width = 20, height = 11, dpi = 300, units = "cm")
