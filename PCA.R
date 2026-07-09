# Load libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Read the raw count data and metadata
dat <- read.csv("raw.csv", header = T, row.names = 1)
info <- read.table("colData.txt", header = T, sep = '\t')

# DESeq2
dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

keep <- rowSums(counts(dds)) >= 10 # Low read filtering
dds <- dds[keep,]

ddsDE <- DESeq(dds)

# For DESeq2 results, used raw & metadata with only 2 conditions (ex: St.36_D vs St.36_V)
res <- results(ddsDE, alpha = 0.05)
res.df <- as.data.frame(res)
res.dfO <- res.df[order(res.df$padj),]
write.csv(res.dfO, "DESeq.csv")

# Generate normalized counts
normCounts <- counts(ddsDE, normalized = T)
normCounts.df <- as.data.frame(normCounts)
write.csv(normCounts.df, "NormCounts.csv")

# PCA
vsd <- vst(dds, blind = FALSE)
#rld <- rlog(dds, blind=FALSE)

a <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

b <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = condition), size = 8, shape = 21, color = "transparent") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 20, family = "Arial"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.1)
  ) +
  guides(fill = "none") +
  scale_fill_manual(values = c(
    "HH14" = "#61CBF4",
    "HH19" = "#215F9A"
  ))

# Calculate the original ranges
x_range <- range(pcaData$PC1)
y_range <- range(pcaData$PC2)

# Calculate the expansion
x_expand <- diff(x_range) * 0.02
y_expand <- diff(y_range) * 0.02

# Define new limits with extra space
x_limits <- c(x_range[1] - x_expand, x_range[2] + x_expand)
y_limits <- c(y_range[1] - y_expand, y_range[2] + y_expand)

# Add Cartesian coordinates with new limits to plot
b <- b + coord_cartesian(xlim = x_limits, ylim = y_limits)
print(b)
ggsave("PCA.png", plot = b, width = 15, height = 15, dpi = 300, units = "cm", limitsize = FALSE)
