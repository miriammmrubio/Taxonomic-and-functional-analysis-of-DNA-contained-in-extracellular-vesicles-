setwd("C:\\Users\\miriam\\TFG\\python tfg")
#Required packages:
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(VennDiagram)
library(utils)
library(writexl)
library(heatmaply)
library(ggdist)
library(stats)
library(ggdendro)
library(readxl)
library(heatmaply)
library(pheatmap)  
library(gplots)
library(dendextend)
library(tibble)
library(gridExtra)
library(dendextend)
library(gridExtra)
library(ggdendroplot)
library(openxlsx)

#HEAT MAP WITH ALL THE PHYLUM (ON X AXIS), THE 20 MOST ABUNDANT KOs OF ALL OF THEM AND COLOURED BY LOG2 OF THE RATIO
#We read the tables generated in Python
#137 m:
tsvheatmapratios137 <- read.table("Tablaheatmapfilosratios.tsv", header = TRUE, sep = "\t", quote = "")

#We convert the 0's into 1's in order to calculate log2.
tsvheatmapratios137[tsvheatmapratios137 == 0] <- 1

# Calculate the log2 of the ratios
numericcolumns <- sapply(tsvheatmapratios137, is.numeric)
tsvheatmapratios137[numericcolumns] <- log2(tsvheatmapratios137[numericcolumns])
tsvheatmapratios137 <- tsvheatmapratios137[, -which(colnames(tsvheatmapratios137) == "Description")]

# Set column 'KO' as an index
rownames(tsvheatmapratios137) <- tsvheatmapratios137$KO
tsvheatmapratios137$KO <- NULL
tsvheatmapratios137 <- as.matrix(tsvheatmapratios137)

heatmap137 <- pheatmap(tsvheatmapratios137, scale = "row", fontsize_row = 15, fontsize_col = 20,
                       cellwidth = 25, cellheight = 15, fontface_row = "bold", fontface_col = "bold")
ggsave("heatmap137.png", heatmap137, width = 10, height = 26)



#The same but for the condition 5-25 metres
tsvheatmapratios255 <- read.table("Tablaheatmapfilosratios255.tsv", header = TRUE, sep = "\t", quote = "")
KO_Des255 <- tsvheatmapratios255[, c("KO", "Description")]
cat(KO_Des255$KO, sep = "\n")
cat(KO_Des255$Description, sep = "\n")
write.csv(KO_Des255, "TablaKODesRatiosFilo255.csv", row.names = FALSE)

tsvheatmapratios255[tsvheatmapratios255 == 0] <- 1

numericcolumns255 <- sapply(tsvheatmapratios255, is.numeric)
tsvheatmapratios255[numericcolumns255] <- log2(tsvheatmapratios255[numericcolumns255])
tsvheatmapratios255 <- tsvheatmapratios255[, -which(colnames(tsvheatmapratios255) == "Description")]


rownames(tsvheatmapratios255) <- tsvheatmapratios255$KO
tsvheatmapratios255$KO <- NULL
tsvheatmapratios255 <- as.matrix(tsvheatmapratios255)

heatmap255 <- pheatmap(tsvheatmapratios255, scale = "row", fontsize_row = 15, fontsize_col = 20,
                       cellwidth = 25, cellheight = 15, fontface_row = "bold", fontface_col = "bold")
ggsave("heatmap255.png", heatmap255, width = 10, height = 24)

