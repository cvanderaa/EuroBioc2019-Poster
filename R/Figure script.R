
####---- Description ----####

# The script generates the figures of the poster for the EuroBioc2019 
# conference.

# For poster tips, check this website: 
# https://www.animateyour.science/post/how-to-design-an-award-winning-conference-poster

# Abstract: Recent advances in sample preparation, processing and mass 
# spectrometry (MS) have allowed the emergence of MS-based single-cell 
# proteomics (SCP). However, bioinformatics tools to process and analyze these 
# new types of data are still missing. In order to boost the development and the 
# benchmarking of SCP methodologies, we are developing the scpdata experiment 
# package. The package will distribute published and curated SCP data sets in 
# standardized Bioconductor format. The poster will give an overview of the 
# available data sets and show preliminary results regarding data exploration 
# and processing.


####---- Setup environment ----####

# Set working directory
setwd("~/tmp/EuroBioc2019-Poster/") # This is the directory to the git repo.
# Load required packages
library(cowplot)
library(export)
library(ggplot2)
library(kableExtra)
library(knitr)
library(grid)
library(gridExtra)
library(magrittr)
library(MSnbase)
library(scpdata)
library(tidyr)
# Load utility functions
source("./R/utils-0.0.1.R")


####---- Content of the package ----####

desc <- scpdata()$result[, -c(1,2), drop=F]
table2tex(desc, file = "./figs/content.tex", standAlone = FALSE)

####---- Data quality control ----####

sc <- specht2019_peptide2
# Keep single cell runs FP94 and FP97
sel <- !grepl("blank|_QC_|col19|col2[0-4]", pData(sc)$run) &
  grepl("FP9[47]", pData(sc)$run)
sc <- sc[, sel]
run <- "190222S_LCA9_X_FP94BF" # pData(sc)$run[1]
sc <- sc[, pData(sc)$run == run]
# Format the data
channel <- pData(sc)$channel
channel <- gsub("Reporter[.]intensity[.]", "", channel)
channel <- paste0("TMT", as.numeric(channel) + 1)
df <- data.frame(run = run, channel = channel, 
                 cellType = pData(sc)$cell_type, t(exprs(sc)))
df <- pivot_longer(data = df, cols = -(1:3), values_to = "intensity") 
# Get counts per channel
df %>%  group_by(channel) %>% 
  summarise(max = max(intensity, na.rm = TRUE), 
            n = sum(!is.na(intensity)),
            cellType = NA,
            mean = mean(intensity, na.rm = TRUE),
            median = median(intensity, na.rm = TRUE)) -> counts
# Create the plot
qc.pl <- ggplot(data = df, aes(x = channel, y = intensity)) +
  geom_violin(data = df, aes(fill = cellType), na.rm = TRUE) + scale_y_log10() + 
  geom_point(data = counts, aes(x = channel, y = median, shape = "+"), 
             color = "red", size = 5) + 
  geom_text(data = counts, aes(x = channel, y = max*2, label = paste0("n=", n)),
            size = 3, color = "grey50") + 
  scale_fill_manual(name = "Well type",
                    values = c(carrier_mix = "grey80", unused = "grey90",
                               norm = "skyblue3", sc_0 = "wheat", 
                               sc_m0 = "#b3ba82", sc_u = "coral"), 
                    limits = c("carrier_mix", "unused", 
                               "norm", "sc_0", "sc_m0"),
                    labels = c(carrier_mix = "carrier (100c)", 
                               norm = "reference (5c)",
                               sc_0 = "empty", unused = "unused", 
                               sc_m0 = "macrophage (1c)")) +
  scale_x_discrete(limits = paste0("TMT", 1:11)) +
  scale_shape_manual(values = "+",
                     labels = c(`+` = "Median"),
                     name = "")
# Save plot
graph2png(x = qc.pl, file = "./figs/QC.png", height = 2.6, width = 7)


####---- Data validation ----####

scpd <- specht2019_peptide
pca <- nipals(exprs(scpd), ncomp = 3, center = TRUE, scale = TRUE)
pca_pl1 <- customPCA(scpd, pca, x = "PC2", y = "PC1", color = "celltype", 
                     shape = "batch_chromatography", size = 1) + 
  theme(legend.position = "none")
pca_pl2 <- customPCA(scpd, pca, x = "PC3", y = "PC1", color = "celltype", 
                     shape = "batch_chromatography", size = 1) + 
  scale_color_manual(name = "Cell type", 
                     values = c(sc_m0 = "skyblue3",
                                sc_u = "coral"),
                     labels = c(sc_m0 = "macrophages",
                                sc_u = "monocytes")) + 
  scale_shape(name = "Batch")
# Save plot
graph2png(x = plot(grid.arrange(pca_pl1, pca_pl2, ncol = 2, widths=c(2.65, 4))), 
          file = "./figs/PCA.png", height = 2.65, width = 7)

####---- Missingness ----####

sc <- specht2019_peptide
# Format the data
df <- do.call(cbind, lapply(unique(pData(sc)$celltype), function(x){
  .sub <- exprs(sc)[, pData(sc)$celltype == x]
  mis <- rowSums(is.na(.sub))/ncol(.sub)*100
  logFC <- apply(.sub, 1, median, na.rm = TRUE)
  out <- data.frame(mis, logFC)
  colnames(out) <- paste0(c("mis", "logFc"), "_", x)
  return(out)
}))
df$relFC <- df$logFc_sc_m0 - df$logFc_sc_u
df$relFC[df$relFC > 2] <- 2
df$relFC[df$relFC < -2] <- -2
# Scatter plot
sp <- ggplot(data = df, aes(x = mis_sc_m0, y = mis_sc_u, col = relFC)) +
  geom_point(size = 0.5) + 
  scale_color_gradient2(name = "log2(rFC)", low = "darkgreen", breaks = c(-2,2), 
                        labels = c("monocyte", "macrophage"),
                        high = "red3", midpoint = 0, mid = "wheat") +
  theme(legend.direction = "horizontal",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.background = element_blank()) +
  ylab("Missingness (%) in monocytes") + xlab("Missingness (%) in macrophages")
# Macrophage density plot
dp1 <- ggplot(data = df, aes(mis_sc_m0)) + 
  geom_density(fill = "grey") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 1), "cm"),
        panel.background = element_blank())
# monocyte density plot
dp2 <- ggplot(data = df, aes(mis_sc_u)) + 
  geom_density(fill = "grey") +
  coord_flip() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 1, 0), "cm"),
        panel.background = element_blank())
# Get legend
leg <- get_legend(sp)
sp <- sp + theme(legend.position = "none")
# Empty plot
blank <- ggplot() + geom_blank() + theme(panel.background = element_blank())
# Combine and save plots
pl.mis <- plot_grid(dp1, blank, sp, dp2, align = "none",
                    ncol = 2, nrow = 2,
                    rel_widths=c(10, 1), rel_heights=c(1, 10))
graph2png(x = pl.mis, file = "./figs/missing.png", aspectr = 1,
          width = 4, height = 4)
graph2png(x = grid.draw(leg), file = "./figs/missing-leg.png", 
          width = 4, height = 1)


