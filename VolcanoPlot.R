# Load libraries
library(ggplot2)
library(ggrepel)

# Read data
v <- read.csv("v.csv") 

# Drop the NA rows right away
v <- subset(v, !is.na(padj))

# Transform padj to -log10(padj)
v$log_padj <- -log10(v$padj)

# Define labeled genes and their colors
label_colors <- c(
  ALDH1A2 = "#61CBF4", RDH10   = "#61CBF4",
  LIN28A  = "#61CBF4", TRIM71  = "#61CBF4", MEOX1   = "#61CBF4", 
  HOXB8   = "#61CBF4", LFNG    = "#61CBF4", PRTG    = "#61CBF4", 
  
  FGF10   = "#215F9A", WNT5A   = "#215F9A", GREM1   = "#215F9A", 
  HOXD11  = "#215F9A", GLI3    = "#215F9A", DUSP6   = "#215F9A", 
  RSPO2   = "#215F9A", SPRY2   = "#215F9A", 
  HOXD10  = "#215F9A", HOXD9   = "#215F9A", HOXD12  = "#215F9A", 
  HOXA10  = "#215F9A", HOXA11  = "#215F9A"
)

# Add label and color columns only for labeled genes
v$label <- ifelse(v$gene %in% names(label_colors), v$gene, NA)
v$label_color <- label_colors[v$label]

# Plot
b <- ggplot(v, aes(x = fc, y = log_padj)) +
  
  # Background points (darker gray)
  geom_point(data = subset(v, is.na(label)),
             color = "gray70", fill = "gray70", shape = 21,
             size = 2, stroke = 0, alpha = 0.5) +
  
  # Labeled milestone points
  geom_point(data = subset(v, !is.na(label)),
             aes(fill = label_color),
             color = "black", shape = 21, size = 3.2, stroke = 0.9) +
  
  # Label fonts and spacing
  geom_text_repel(data = subset(v, !is.na(label)),
                  aes(label = paste0("italic('", gene, "')")),
                  parse = TRUE, 
                  size = 7,
                  family = "Arial",
                  max.overlaps = Inf,
                  box.padding = 0.7,     
                  point.padding = 0.4,   
                  segment.size = 0.4, 
                  segment.color = "black",
                  min.segment.length = 0) + 
  
  # Axes titles
  xlab("log2 Fold Change (HH19/HH14)") +
  ylab("-log10 padj") +
  scale_fill_identity() +
  
  
  theme_minimal(base_family = "Arial") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 17)
  )

# Save high-res plot
ggsave("volcano.png", plot = b, width = 22, height = 20, dpi = 300, units = "cm")
