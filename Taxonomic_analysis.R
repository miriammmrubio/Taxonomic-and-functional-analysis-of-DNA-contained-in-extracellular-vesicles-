## This R script have been created following these two tutorials
# 25/1/2023 

## Phyloseq tutorials
# https://mvuko.github.io/meta_phyloseq/#fn1
# https://carpentries-incubator.github.io/metagenomics/aio/index.html

#Required files
#OTU table
# OTUs.txt

#Taxonomy table
# taxonomy.txt

#Metadata table
#metadata.txt
#...................................................................

#Required packages
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(gplots)
library(VennDiagram)
library(ggforce)
library(reshape2)
library(ragg)
library(png)

# setting working directory
setwd("C:\\Users\\miriam\\TFG")

#1. PHYLOSEQ OBJECT


table <- read.table("EV_contig_table_aggregated_sin_n_groups.txt", sep = "\t", header = TRUE)
table$OTU <- paste("OTU", 1:nrow(table), sep = "")


##OTU TABLE
otu <- select(table, -(Tax))
otu <- otu[, c("OTU", setdiff(names(otu), "OTU"))]
rownames(otu) <- otu[,1]
otu[,1] <- NULL
str(otu)
colnames(otu) <- gsub("^X", "",  colnames(otu))

##TAX TABLE
tax_table <- table[c("OTU", "Tax")]
rownames(tax_table) <- tax_table[,1]
tax_table[,1] <- NULL
tax_table$Kingdom <- NA
tax_table$Phylum <- NA
tax_table$Class <- NA
tax_table$Order <- NA
tax_table$Family <- NA
tax_table$Genus <- NA
tax_table$Species <- NA

for (i in 1:nrow(tax_table)) {
  levels <- unlist(strsplit(as.character(tax_table[i, "Tax"]), ";"))
  for (j in 1:length(levels)) {
    level <- trimws(levels[j])
    if (startsWith(level, "k_")) {
      tax_table[i, "Kingdom"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "p_")) {
      tax_table[i, "Phylum"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "c_")) {
      tax_table[i, "Class"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "o_")) {
      tax_table[i, "Order"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "f_")) {
      tax_table[i, "Family"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "g_")) {
      tax_table[i, "Genus"] <- substr(level, 3, nchar(level))
    } else if (startsWith(level, "s_")) {
      tax_table[i, "Species"] <- substr(level, 3, nchar(level))
    } 
  }
}
tax_table[,1] <- NULL


taxmat <- as.matrix(tax_table)

#SAMPLE DATA
samples <- c("137m_V", "137m_M", "25m_V", "25m_M", "5m_V", "5m_M")
biome <- c("vesicles", "cells", "vesicles", "cells", "vesicles", "cells")
depth <- c(137, 137, 25, 25, 5, 5)
metadata <- data.frame(samples,biome,depth)
rownames(metadata) <- metadata[,1]
metadata[,1] <- NULL

#converting tables to phyloseq objects
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
sampledata = sample_data(metadata)

#Generating Physeq object
physeq <- phyloseq(OTU, TAX, sampledata)
physeq


#2. ANALYSIS AT KINGDOM LEVEL 
#KINGDOM PLOT: ABSOLUTE
colors1 <- c("#FFA600", "#008B8B", "#6A5ACD", "#FFFF00")
taxglomk  <- tax_glom(physeq, "Kingdom")
kingdomplot <- plot_bar(taxglomk, fill ="Kingdom") +
  scale_fill_manual(values = colors1) + 
  theme_classic()+
  theme(legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.key.size = unit(1.25, "cm"),
        axis.title.x = element_text(size = 27, face = "bold", color = "black"),
        axis.title.y = element_text(size = 27, face = "bold", color = "black"),
        axis.text.x = element_text(size = 23, face = "bold", color = "black"),
        axis.text.y = element_text(size = 23, face = "bold", color = "black")) + 
  ylab("Abundance") 
ggsave("kingdomplot.png", plot = kingdomplot, width = 15, height = 10)


#KINGDOM PLOT: RELATIVE
physeq_relabund <- transform_sample_counts(physeq, function(x) x/sum(x))
taxglomk_relabund  <- tax_glom(physeq_relabund, "Kingdom")
kingdomplot_relabund <- plot_bar(taxglomk_relabund, fill ="Kingdom") +
  scale_fill_manual(values = colors1) + 
  theme_classic()+
  theme(legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.key.size = unit(1.25, "cm"),
        axis.title.x = element_text(size = 27, face = "bold", color = "black"),
        axis.title.y = element_text(size = 27, face = "bold", color = "black"),
        axis.text.x = element_text(size = 23, face = "bold", color = "black"),
        axis.text.y = element_text(size = 23, face = "bold", color = "black")) + 
  ylab("Relative abundance") 
kingdomplot_relabund
ggsave("KingdomPlotrel.png", plot = kingdomplot_relabund, width = 15, height = 10)


#Filter phyloseq object, selecting the kingdoms bacteria and archaea
tax_table_original <- tax_table(physeq)
tax_table_filtered <- tax_table_original[tax_table_original[, "Kingdom"] %in% c("Bacteria", "Archaea"), ]
selected_otus <- rownames(tax_table_filtered)
bac_arch_physeq <- prune_taxa(selected_otus, physeq)

#Phyloseq objet: cells and vesicles 
cells_physeq <- subset_samples(bac_arch_physeq, biome == "cells")
vesicles_physeq <- subset_samples(bac_arch_physeq, biome == "vesicles")

#Phyloseq (bacteria and archaea) object with data in percent 
percentages <- transform_sample_counts(bac_arch_physeq, function(x) x*100 / sum(x) )

#3. ANALYSIS AT PHYLUM LEVEL

#VENN DIAGRAM: vesicles
otuvenn <- otu_table(vesicles_physeq)
otu3 <- data.frame(t(otuvenn))
#137 m sample
list137 <- colnames(otu3["137m_V", apply(otu3["137m_V",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
#25 m sample
list25 <- colnames(otu3["25m_V", apply(otu3["25m_V",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
#5 m sample
list5 <- colnames(otu3["5m_V", apply(otu3["5m_V",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
taxvenn <- tax_table(bac_arch_physeq)
indices_137 <- match(list137, rownames(taxvenn))
indices_25 <- match(list25, rownames(taxvenn))
indices_5 <- match(list5, rownames(taxvenn))
phylum137 <- unique(as.character(taxvenn[indices_137, "Phylum"]))
phylum25 <- unique(as.character(taxvenn[indices_25, "Phylum"]))
phylum5 <- unique(as.character(taxvenn[indices_5, "Phylum"]))
par(mfrow=c(1,1))
venn(list(phylum137, phylum25, phylum5))
plot.new()
par(mfrow=c(1,1), mar=c(0,0,0,0))
my_colors <- colorRampPalette(c("purple", "orange", "lightblue"))(3)
Venndiagram <- draw.triple.venn(
  area1 = 3+5+2+82, area2 = 5+9+82+9, area3 = 82+9+3+2, n12 = 5+82, n23 = 82+9, n13 = 2+82, n123 = 82, category = c("137m_V", "25m_V", "5m_V"), 
  lty = "blank", fill = my_colors, cex = 3.5, fontface = 2,
  cat.cex = 2.5, cat.fontface = 2
  )
ggsave("VennDiagram.png", plot = Venndiagram, width = 10, height = 10)

#VENN DIAGRAM: cells
otuvennc <- otu_table(cells_physeq)
otu3c <- data.frame(t(otuvennc))
#137
lista137c <- colnames(otu3c["137m_M", apply(otu3c["137m_M",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
#25
lista25c <- colnames(otu3c["25m_M", apply(otu3c["25m_M",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
#5 
lista5c <- colnames(otu3c["5m_M", apply(otu3c["5m_M",], MARGIN=2, function(x) sum(x) > 15 && any(x > 0))])
taxvennc <- tax_table(bac_arch_physeq)
indices_fila137c <- match(lista137c, rownames(taxvennc))
indices_fila25c <- match(lista25c, rownames(taxvennc))
indices_fila5c <- match(lista5c, rownames(taxvennc))
filos137c <- unique(as.character(taxvennc[indices_fila137c, "Phylum"]))
filos25c <- unique(as.character(taxvennc[indices_fila25c, "Phylum"]))
filos5c <- unique(as.character(taxvennc[indices_fila5c, "Phylum"]))
par(mfrow=c(1,1))
venn(list(filos137c, filos25c, filos5c))
plot.new()
par(mfrow=c(1,1), mar=c(0,0,0,0))
my_colors <- colorRampPalette(c("purple", "orange", "lightblue"))(3)
Venndiagramc <- draw.triple.venn(area1 = 12+1+5+86, area2 = 1+2+86+8, area3 = 86+8+5+1, n12 = 1+86, n23 = 86+8, n13 = 5+86, n123 = 86, category = c("137m_M", "25m_M", "5m_M"), 
                                 lty = "blank", fill = my_colors, cex = 3.5, fontface = 2,
                                 cat.cex = 2.5, cat.fontface = 2)
ggsave("VennDiagramcells.png", plot = Venndiagramc, width = 10, height = 10)


#Representation of the 10 most abundant phylum in different samples

colors2 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#A6A6A6", "#FFFF33", "#A65628", "#F781BF", "#66C2A5")
#Sort the Phyla by abundance and pick the top 10
top10P.names = sort(tapply(taxa_sums(bac_arch_physeq), tax_table(bac_arch_physeq)[, "Phylum"], sum), TRUE)[1:10]
top_10P2 <- subset_taxa(bac_arch_physeq, Phylum %in% names(top10P.names))
glom10P2 <- tax_glom(top_10P2,"Phylum")
view(otu_table(glom10P2))
top10phylum <- plot_bar(glom10P2, x="Sample", fill="Phylum", facet_grid = ~Phylum) + 
  scale_fill_manual(values = colors2) +
  theme_classic()+
  geom_bar(aes(fill=Phylum), stat="identity", position="fill")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
        legend.key.size = unit(1.75, "cm"),
        axis.title.x = element_text(size = 45, face = "bold", color = "black"),
        axis.title.y = element_text(size = 45, face = "bold", color = "black"),
        axis.text.x = element_text(size = 40, face = "bold", color = "black"),
        axis.text.y = element_text(size = 40, face = "bold", color = "black"))
ggsave("Top10phylum.png", plot = top10phylum, width = 42, height = 25)


#Representation of the logarithm in base 2 of the ratio "Abundance in cell fraction/abundance in vesicle fraction" of the 10 phyla with the highest abundance in the studied samples. 
ratio_137 <- otu_table(glom10P2)[, "137m_V"] / otu_table(glom10P2)[, "137m_M"]
tax_table(glom10P2)
ratio_25 <- otu_table(glom10P2)[, "25m_V"] / otu_table(glom10P2)[, "25m_M"]
ratio_5 <- otu_table(glom10P2)[, "5m_V"] / otu_table(glom10P2)[, "5m_M"]
ratio_137_log2 <- log2(ratio_137)
ratio_25_log2 <- log2(ratio_25)
ratio_5_log2 <- log2(ratio_5)
ratios_df <- data.frame(
  OTU = rownames(otu_table(glom10P2)[1:10]),
  Ratio_137 = ratio_137_log2[1:10],
  Ratio_25 = ratio_25_log2[1:10],
  Ratio_5 = ratio_5_log2[1:10]
)
txratios <- tax_table(glom10P2)
taxratios <- as.data.frame(txratios)
taxratios$OTU <- rownames(taxratios)
ratios_df_merged <- merge(ratios_df, taxratios, by = "OTU")


ratios_long <- gather(ratios_df_merged, key = "ratio", value = "valor", Ratio_137:Ratio_5)
ratios_long$ratio <- gsub("Ratio_", "", ratios_long$ratio)
ratios <- ggplot(ratios_long, aes(x = ratio, y = valor, fill = Phylum)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic()+
  facet_grid(. ~ Phylum, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Ratio") +
  scale_fill_manual(values = colors2) +
  theme(strip.text.x = element_text(size = 10, face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
        legend.key.size = unit(1.25, "cm"),
        axis.title.x = element_text(size = 38, face = "bold", color = "black"),
        axis.title.y = element_text(size = 38, face = "bold", color = "black"),
        axis.text.x = element_text(size = 38, face = "bold", color = "black"),
        axis.text.y = element_text(size = 38, face = "bold", color = "black"))+
  scale_x_discrete(labels = c("137m", "25m", "5m"))
ggsave("Phylumrratios.png", plot = ratios, width = 25, height = 13)

#Plot the last two graphs with the vesicle data. 
#Representation of the 10 most abundant phylum in different samples
top10P.namesv = sort(tapply(taxa_sums(vesicles_physeq), tax_table(vesicles_physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P.namesv
top10Pv <- subset_taxa(bac_arch_physeq, Phylum %in% names(top10P.namesv))
glom10Pv <- tax_glom(top10Pv,"Phylum")
top10phylumvesicles <- plot_bar(glom10Pv, x="Sample", fill="Phylum", facet_grid = ~Phylum) + 
  scale_fill_manual(values = colors2) +
  geom_bar(aes(fill=Phylum), stat="identity", position="fill")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
        legend.key.size = unit(1.75, "cm"),
        axis.title.x = element_text(size = 45, face = "bold", color = "black"),
        axis.title.y = element_text(size = 45, face = "bold", color = "black"),
        axis.text.x = element_text(size = 40, face = "bold", color = "black"),
        axis.text.y = element_text(size = 40, face = "bold", color = "black"))
ggsave("Top10phylum_vesicles.png", plot = top10phylumvesicles, width = 42, height = 25)

#Representation of the logarithm in base 2 of the ratio "Abundance in cell fraction/abundance in vesicle fraction" of the 10 phyla with the highest abundance in the studied samples. 
ratio_137v  <- otu_table(glom10Pv)[, "137m_V"] / otu_table(glom10Pv)[, "137m_M"]
ratio_25v <- otu_table(glom10Pv)[, "25m_V"] / otu_table(glom10Pv)[, "25m_M"]
ratio_5v <- otu_table(glom10Pv)[, "5m_V"] / otu_table(glom10Pv)[, "5m_M"]
ratio_137_log2v <- log2(ratio_137v)
ratio_25_log2v  <- log2(ratio_25v)
ratio_5_log2v <- log2(ratio_5v)
ratios_dfv <- data.frame(
  OTU = rownames(otu_table(glom10Pv)[1:10]),
  Ratio_137v = ratio_137_log2v[1:10],
  Ratio_25v = ratio_25_log2v[1:10],
  Ratio_5v = ratio_5_log2v[1:10]
)
ratios_dfv
txratiosv <- tax_table(glom10Pv)
taxratiosv <- as.data.frame(txratiosv)
taxratiosv$OTU <- rownames(taxratiosv)
ratios_df_mergedv <- merge(ratios_dfv, taxratiosv, by = "OTU")
ratios_longv <- gather(ratios_df_mergedv, key = "ratio", value = "valor", Ratio_137v:Ratio_5v)
ratios_longv$ratio <- gsub("Ratio_", "", ratios_longv$ratio)

ratiosvesicles <- ggplot(ratios_longv, aes(x = ratio, y = valor, fill = Phylum)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic()+
  facet_grid(. ~ Phylum, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Ratio") +
  scale_fill_manual(values = colors2) +
  theme(strip.text.x = element_text(size = 10, face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
        legend.key.size = unit(1.25, "cm"),
        axis.title.x = element_text(size = 38, face = "bold", color = "black"),
        axis.title.y = element_text(size = 38, face = "bold", color = "black"),
        axis.text.x = element_text(size = 38, face = "bold", color = "black"),
        axis.text.y = element_text(size = 38, face = "bold", color = "black"))+
  scale_x_discrete(labels = c("137m", "25m", "5m"))
ggsave("Phylumratios_vesicles.png", plot = ratiosvesicles, width = 25, height = 13)


#4. ANALYSIS: FAMILY LEVEL
#We explore different family-level phylum
#How to explore a linage
physeq1 <- subset_taxa(physeq, Family != "")

#PLOT RELATIVE ABUNDANCE
# Function
generate_plotrelative <- function(data, y_cutoff, fill_colors, filename) {
  data$Family <- as.factor(data$Family)
  levels(data$Family) <- c(levels(data$Family), paste("Family <", y_cutoff, "abund"))
  
  data$Family[data$Abundance < y_cutoff] <- paste("Family <", y_cutoff, "abund")
  
  plot <- ggplot(data = data, aes(x = Sample, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    scale_fill_manual(values = fill_colors) +
    theme(
      legend.text = element_text(size = 25, face = "bold", color = "black"),
      legend.title = element_text(size = 25, face = "bold", color = "black"),
      legend.key.size = unit(1.25, "cm"),
      axis.title.x = element_text(size = 25, face = "bold", color = "black"),
      axis.title.y = element_text(size = 25, face = "bold", color = "black"),
      axis.text.x = element_text(size = 22, face = "bold", color = "black"),
      axis.text.y = element_text(size = 22, face = "bold", color = "black")
    ) +
    ylab("Relative abundance")
  
  ggsave(filename, plot = plot, width = 17, height = 12)
}

#we apply the relative abundance function to the most abundant phylum
# Cyanobacteria
cyanos <- subset_taxa(physeq1, Phylum == "Cyanobacteria")
unique(cyanos@tax_table@.Data[,2])
cyanos_percentages <- transform_sample_counts(cyanos, function(x) x * 100 / sum(x))
cyanos_glom <- tax_glom(cyanos_percentages, taxrank = "Family")
cyanos_df <- psmelt(cyanos_glom)
colors_cf <- c("#E0FFFF", "#fdc086", "#386cb0", "#f0027f")
generate_plotrelative(cyanos_df, 1, colors_cf, "Cyanosrelative_family.png")

#Proteobacteria
proteo <- subset_taxa(physeq1, Phylum == "Proteobacteria")
unique(proteo@tax_table@.Data[,2])
proteo_percentages <- transform_sample_counts(proteo, function(x) x*100 / sum(x) )
proteo_glom <- tax_glom(proteo_percentages, taxrank = "Family")
proteo_df <- psmelt(proteo_glom)
colorsf  <- c("#7fc97f", "#cab2d6", "#fdc086", "#f0027f", "#386cb0", "#ADD8E6" , "#ffff99", "#d95f02", "#D3D3D3", "#7570b3",  "#bf5b17", "#e6ab02", "#a6761d", "#666666", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#33a02c", "#6a3d9a")
generate_plotrelative(proteo_df, 6, colorsf, "Proteorelative_family.png")

#Bacteroidetes
bacteroidetes <- subset_taxa(physeq1, Phylum == "Bacteroidetes")
unique(bacteroidetes@tax_table@.Data[,2])
bactereo_percentages <- transform_sample_counts(bacteroidetes, function(x) x*100 / sum(x) )
bactereoo_glom <- tax_glom(bactereo_percentages, taxrank = "Family")
bactereo_df <- psmelt(bactereoo_glom)
colorsbf <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#33a02c", "#bf5b17", "#666666", "#E0FFFF", "#D3D3D3", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(bactereo_df, 1, colorsbf, "Bactereorelative_family.png")

#Actinobacteria
actino <- subset_taxa(physeq1, Phylum == "Actinobacteria")
unique(actino@tax_table@.Data[,2])
actino_percentages <- transform_sample_counts(actino, function(x) x*100 / sum(x) )
actino_glom <- tax_glom(actino_percentages, taxrank = "Family")
actino_df <- psmelt(actino_glom)
colorsaf <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#33a02c", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(actino_df, 1, colorsaf, "Actinorelative_family.png")

#Firmicutes
firmic <- subset_taxa(physeq1, Phylum == "Firmicutes")
unique(firmic@tax_table@.Data[,2])
firmi_percentages <- transform_sample_counts(firmic, function(x) x*100 / sum(x) )
firmi_glomf <- tax_glom(firmi_percentages, taxrank = "Family")
firmi_dff <- psmelt(firmi_glomf)
colorsff <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#33a02c", "#D3D3D3", "#666666", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(firmi_dff, 3, colorsff, "Firmirelative_family.png")

#Euryarchaeota
eury <- subset_taxa(physeq1, Phylum == "Euryarchaeota")
unique(eury@tax_table@.Data[,2])
eury_percentages <- transform_sample_counts(eury, function(x) x*100 / sum(x) )
eury_glomf <- tax_glom(eury_percentages, taxrank = "Family")
eury_dff <- psmelt(eury_glomf)
colorsfg <- c("#7fc97f", "#fdc086", "#ffff99", "#33a02c", "#386cb0", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(eury_dff, 3, colorsfg, "Euryrelative_family.png")

#Balneolaeota
bal <- subset_taxa(physeq1, Phylum == "Balneolaeota")
unique(bal@tax_table@.Data[,2])
bal_percentages <- transform_sample_counts(bal, function(x) x*100 / sum(x) )
bal_glomf <- tax_glom(bal_percentages, taxrank = "Family")
bal_dff <- psmelt(bal_glomf)
colorsbr <- c("#7fc97f", "#fdc086", "#fb9a99","#ffff99", "#33a02c", "#386cb0", "#D3D3D3", "#666666", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(bal_dff, 1.5, colorsbr, "Balnerelative_family.png")

#Chloroflexi
physeq2 <- subset_taxa(physeq, Family != "")
chloro <- subset_taxa(physeq2, Phylum == "Chloroflexi")
unique(chloro@tax_table@.Data[,2])
chloro_percentages <- transform_sample_counts(chloro, function(x) x*100 / sum(x) )
chloro_glomf <- tax_glom(chloro_percentages, taxrank = "Family")
chloro_dff <- psmelt(chloro_glomf)
colorsfc <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#33a02c", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(chloro_dff, 2, colorsfc, "Chlororelative_family.png")

#Thaumarchaeota
thau <- subset_taxa(physeq2, Phylum == "Thaumarchaeota")
unique(thau@tax_table@.Data[,2])
thau_percentages <- transform_sample_counts(thau, function(x) x*100 / sum(x) )
thau_glom <- tax_glom(thau_percentages, taxrank = "Family")
thau_df <- psmelt(thau_glom)
colorstf <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f", "#bf5b17", "#666666", "#E0FFFF", "#7570b3", "#d95f02", "#e6ab02", "#a6761d", "#D3D3D3", "#ADD8E6", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(thau_df, 0, colorstf, "Thaurelative_family.png")

#Planctomycetes
plan <- subset_taxa(physeq2, Phylum == "Planctomycetes")
unique(plan@tax_table@.Data[,2])
plan_percentages <- transform_sample_counts(plan, function(x) x*100 / sum(x) )
plan_glom <- tax_glom(plan_percentages, taxrank = "Family")
plan_df <- psmelt(plan_glom)
colorspf <- c("#7fc97f", "#beaed4", "#fdc086", "#b2df8a", "#386cb0","#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_plotrelative(plan_df, 1, colorspf, "Planctorelative_family.png")

#Function absolute abundance 
generate_absolute_plot <- function(data, y_cutoff, fill_colors, filename) {
  data$Family <- as.factor(data$Family)
  data$Family <- ifelse(data$Abundance < y_cutoff, paste("Family <", y_cutoff), as.character(data$Family))
  
  plot <- ggplot(data = data, aes(x = Sample, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    scale_fill_manual(values = fill_colors) +
    theme(
      legend.text = element_text(size = 25, face = "bold", color = "black"),
      legend.title = element_text(size = 25, face = "bold", color = "black"),
      legend.key.size = unit(1.25, "cm"),
      axis.title.x = element_text(size = 25, face = "bold", color = "black"),
      axis.title.y = element_text(size = 25, face = "bold", color = "black"),
      axis.text.x = element_text(size = 22, face = "bold", color = "black", angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 22, face = "bold", color = "black")
    )
  
  ggsave(filename, plot = plot, width = 17, height = 12)
}
#Cyanobacteria
absolutecyanos_glom <- tax_glom(physeq = cyanos, taxrank = "Family")
absolutecyanos_df <- psmelt(absolutecyanos_glom)
absolutecyanos_df$Family <- as.factor(absolutecyanos_df$Family)
absolutecyanos_df$Family <- ifelse(absolutecyanos_df$Abundance < 1000, "Family < 1000", as.character(absolutecyanos_df$Family))
colorscfa <- colores_contrastecf <- c("#7fc97f", "#E0FFFF", "#fdc086", "#386cb0", "#f0027f")
generate_absolute_plot(absolutecyanos_df, 1000, colorscfa, "Cyanosabsolute_family.png")

#Proteobacteria
absoluteproteo_glom <- tax_glom(physeq = proteo, taxrank = "Family")
absoluteproteo_df <- psmelt(absoluteproteo_glom)
absoluteproteo_df$Family <- as.factor(absoluteproteo_df$Family)
absoluteproteo_df$Family <- ifelse(absoluteproteo_df$Abundance < 150000, "Family < 150000", as.character(absoluteproteo_df$Family))
colorsff  <- c("#7fc97f", "#cab2d6", "#fdc086", "#b2df8a","#f0027f", "#386cb0", "#ADD8E6", "#ffff99", "#d95f02", "#D3D3D3", "#7570b3",  "#bf5b17", "#e6ab02", "#a6761d", "#666666", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#33a02c", "#6a3d9a")
generate_absolute_plot(absoluteproteo_df, 150000, colorsff, "Proteoabsolute_family.png")

#Bacteroidetes
absolutebactereo_glom <- tax_glom(physeq = bacteroidetes, taxrank = "Family")
absolutebactereo_df <- psmelt(absolutebactereo_glom)
absolutebactereo_df$Family <- as.factor(absolutebactereo_df$Family)
absolutebactereo_df$Family <- ifelse(absolutebactereo_df$Abundance < 8000, "Family < 8000", as.character(absolutebactereo_df$Family))
colorsbfa <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#33a02c", "#bf5b17", "#666666", "#E0FFFF", "#D3D3D3", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absolutebactereo_df, 8000, colorsbfa, "Bacteroabsolute_family.png")


#Actinobacteria 
absoluteactino_glom <- tax_glom(physeq = actino, taxrank = "Family")
absoluteactino_df <- psmelt(absoluteactino_glom)
absoluteactino_df$Family <- as.factor(absoluteactino_df$Family)
absoluteactino_df$Family <- ifelse(absoluteactino_df$Abundance < 1500, "Family < 1500", as.character(absoluteactino_df$Family))
colorsfa <- c("#7fc97f", "#fb9a99", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absoluteactino_df, 1500, colorsfa, "Actinoabsolute_family.png")

#Firmicutes
absolutefirmi_glom <- tax_glom(physeq = firmic, taxrank = "Family")
absolutefirmi_df <- psmelt(absolutefirmi_glom)
absolutefirmi_df$Family <- as.factor(absolutefirmi_df$Family)
absolutefirmi_df$Family <- ifelse(absolutefirmi_df$Abundance < 1000, "Family < 1000", as.character(absolutefirmi_df$Family))
colorsfaf <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#fb9a99", "#33a02c", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absolutefirmi_df, 1000, colorsfaf, "Firmiabsolute_family.png")

#Euryarchaeota
absoluteeury_glom <- tax_glom(physeq = eury, taxrank = "Family")
absoluteeury_df <- psmelt(absoluteeury_glom)
absoluteeury_df$Family <- as.factor(absoluteeury_df$Family)
absoluteeury_df$Family <- ifelse(absoluteeury_df$Abundance < 10, "Family < 10", as.character(absoluteeury_df$Family))
colorsefa <- c("#7fc97f", "#fdc086", "#fb9a99","#ffff99", "#33a02c", "#386cb0", "#D3D3D3", "#666666", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absoluteeury_df, 10, colorsefa, "Euryabsolute_family.png")

#Balneolaeota
absolutebalne_glom <- tax_glom(physeq = bal, taxrank = "Family")
absolutebalne_df <- psmelt(absolutebalne_glom)
absolutebalne_df$Family <- as.factor(absolutebalne_df$Family)
colorsfab <- c("#7fc97f", "#fdc086", "#ffff99","#386cb0", "#33a02c", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absolutebalne_df, 0, colorsfab, "Balneabsolute_family.png")

#Chloroflexi
absolutechloro_glom <- tax_glom(physeq = chloro, taxrank = "Family")
absolutechloro_df <- psmelt(absolutechloro_glom)
str(absolutechloro_df)
absolutechloro_df$Family <- as.factor(absolutechloro_df$Family)
absolutechloro_df$Family <- ifelse(absolutechloro_df$Abundance < 50, "Family < 50", as.character(absolutechloro_df$Family))
colorsfch <- c("#7fc97f", "#beaed4","#fdc086", "#ffff99","#386cb0", "#33a02c", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absolutechloro_df, 50, colorsfch, "Chloroabsolute_family.png")

#Thaumarchaeota 
absolutethau_glom <- tax_glom(physeq = thau, taxrank = "Family")
absolutethau_df <- psmelt(absolutethau_glom)
str(absolutethau_df)
absolutethau_df$Family <- as.factor(absolutethau_df$Family)
colorsft <- c("#7fc97f", "#beaed4","#fdc086", "#ffff99","#386cb0", "#33a02c", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absolutethau_df, 0, colorsft, "Thauabsolute_family.png")

#Planctomycetes
#Absolute abundance
absoluteplan_glom <- tax_glom(physeq = plan, taxrank = "Family")
absoluteplan_df <- psmelt(absoluteplan_glom)
str(absoluteplan_df)
absoluteplan_df$Family <- as.factor(absoluteplan_df$Family)
absoluteplan_df$Family <- ifelse(absoluteplan_df$Abundance < 500, "Family < 500", as.character(absoluteplan_df$Family))
colorsfap <- c("#7fc97f","#beaed4", "#fdc086", "#b2df8a","#386cb0", "#D3D3D3", "#666666", "#E0FFFF", "#d95f02", "#e6ab02", "#a6761d", "#7570b3", "#ADD8E6", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
generate_absolute_plot(absoluteplan_df, 500, colorsfap, "Planabsolute_family.png")


#5. ANALYSIS: SPECIE LEVEL
#BETA DIVERSITY
# To make this transformation to percentages, 
# we will take advantage of a function of Phyloseq
percentages <- transform_sample_counts(bac_arch_physeq, function(x) x*100 / sum(x) )
head(percentages@otu_table@.Data)

#beta diversity
# display all the possible distance metrics that Phyloseq can use
distanceMethodList

# We will use Bray-curtis since it is one of the most robust and 
# widely used distance metrics to calculate beta diversity.

# We will use the Phyloseq command ordinate() to generate a new object 
# where the distances between our samples will be allocated after calculating them. 
# For this command, we need to specify which method we will use to generate a matrix. 
# In this example, we will use Non-Metric Multidimensional Scaling or NMDS. 
# NMDS attempts to represent the pairwise dissimilarity between objects 
# in a low-dimensional space, in this case, a two-dimensional plot

meta_ord <- ordinate(physeq = percentages, method = "NMDS", 
                     distance = "bray")
meta_ord

beta <- plot_ordination(physeq = percentages, ordination = meta_ord)
ggsave("beta_diversity.png", plot = beta, width = 15, height = 5)
#We cannot do the analysis because there are too few species in the vesicle samples.

## ALPHA DIVERSITY BACTERIA AND ARCHAEA
a_diversitybacarch <- plot_richness(bac_arch_physeq, measures = c("Fisher","Chao1","Shannon"), color = "biome") +
  geom_boxplot(aes(fill = biome), alpha=.7) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) + 
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "gray95", color = NA))+
  theme(legend.text = element_text(size = 22, face = "bold", color = "black"),
        legend.title = element_text(size = 22, face = "bold", color = "black"),
        axis.title.x = element_text(size = 21, face = "bold", color = "black"),
        axis.title.y = element_text(size = 21, face = "bold", color = "black"),
        axis.text.x = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"))
a_diversitybacarch
ggsave("alpha_diversityBacteriaArchaea.png", plot = a_diversitybacarch, width = 15, height = 8)


#PCOA: BRAY PLOT 
color_palette <- c("#1b9e77", "#d95f02", "#7570b3")
bray_ordinate <- ordinate(physeq = percentages, method = "PCoA", distance = "bray")
bray_plot <-  plot_ordination(bac_arch_physeq, bray_ordinate, shape = "biome", title = "Bray Curtis", axes = c(1,2)) +
  geom_point(aes(color = factor(depth)), size = 6) +
  scale_color_manual(values = color_palette) + # especificar la paleta de colores
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "cm"), legend.margin = margin(0, 0, 0, 0))+
  theme(legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.title = element_text(size = 20, face = "bold", color = "black"),
        axis.title.x = element_text(size = 18, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold", color = "black"),
        axis.text.x = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"))+
  labs(color = "depth", x = "PC1 (84.4%)", y = "PC2 (11.5%)")+
  coord_fixed()
bray_plot
ggsave("brayplot.png", plot = bray_plot, width = 10, height = 5)

#BUBBLE PLOT
##Bubble plot of the 30 most abundant species with data from all vesicle samples
#We selected the 30 most abundant OTUs from the samples.
desired_phylum <- c("Actinobacteria", "Bacteroidetes", "Candidatus Marinimicrobia", "Chloroflexi", "Cyanobacteria", "Euryarchaeota", "Planctomycetes", "Proteobacteria", "Thaumarchaeota", "Firmicutes", "Balneolaeota")
topspeciesv <- subset_taxa(vesicles_physeq, Species != "" & !grepl("bacterium", Species))
top30speciesv <- names(sort(taxa_sums(topspeciesv), TRUE)[1:30])
#We create phyloseq object only with data from the 30 most abundant OTUs.
physeqbubblev <- prune_taxa(top30speciesv, bac_arch_physeq)
#Prepare the tables
otububblev <- otu_table(physeqbubblev)
otububblev <- cbind(OTU=rownames(otububblev), otububblev)
txv <- tax_table(physeqbubblev)
taxv <- as.data.frame(txv)
taxv$OTU <- rownames(taxv)
dfv <- merge(otububblev, taxv, by = "OTU")
dfv <- subset(dfv, Species != "")
dfv <- select(dfv, -Kingdom, -Class, -Order, -Genus)
dfv$OTU <- NULL

#generates a long-format df
df_longv <- gather(dfv, key = "Sample", value = "valor", -Species, -Phylum, -Family)
df_longv$valor <- as.numeric(df_longv$valor)
df_longv$log_valor <- log10(df_longv$valor)
df_longv <- df_longv[df_longv$Phylum %in% desired_phylum, ]
df_longv$Filo_Familia <- paste(df_longv$Phylum, df_longv$Family)
df_longv <- df_longv[order(df_longv$Filo_Familia, df_longv$Species),]
ordenspecies_family <- unique(df_longv$Species)

#Bubble plot
colorsbubble <- c("#0074D9", "#FFC300", "#FF5733", "#00BFFF", "#900C3F", "#2ECC40","#FFC0CB", "#39CCCC", "#FFA500")
legend_order <- c("Proteobacteria", "Firmicutes","Euryarchaeota", "Cyanobacteria", "Balneolaeota", "Bacteroidetes", "Actinobacteria")
bubbleplotv <- ggplot(df_longv, aes(x = Sample, y = reorder(Species, match(Species, ordenspecies_family)), size = log_valor, color = Phylum)) +
  geom_point() +
  scale_size(range = c(3, 45)) +
  scale_color_manual(values = colorsbubble, breaks = legend_order) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(vjust = -0.25)) +
  labs(y = "") +
  theme(legend.text = element_text(size = 50, face = "bold", color = "black"),
        legend.title = element_text(size = 50, face = "bold", color = "black"),
        axis.title.x = element_text(size = 50, face = "bold", color = "black"),
        axis.title.y = element_text(size = 50, face = "bold", color = "black"),
        axis.text.x = element_text(size = 40, face = "bold", color = "black"),
        axis.text.y = element_text(size = 40, face = "bold", color = "black"))+
  guides(color = guide_legend(override.aes = list(size=24)))
bubbleplotv

ggsave("bubbleplot.png", plot = bubbleplotv, width = 45, height = 40)

