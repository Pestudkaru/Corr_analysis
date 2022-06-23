# Co-expression analysis using Bio GPS data
# Author G.K.

# Clean the environment
rm(list = ls())
# Set global options
options(stringsAsFactors = F) # cos factors by default are stupid! 

# Load libraries
library("pheatmap")
library("dendsort")
library("ggplot2")
library("ggpubr")
library("gplots")
library("plotly")
library("cowplot")
library("tidyverse") # For all your data wrangling needs
library("Hmisc")

# Average data values from bio GPS
data <- read.csv(file = "/Home Office/Data/2020.12.10.Che project/Data/BioGPS data/U133AGNF1B.gcrma.avg.csv")

# Use the gene names as row names
rownames(data) <- data$X
data <- data[-1]
data[1:5,1:5]

############################################################
###########Data Wrangling & Filtering ######################
############################################################

# Remove all rows with less than n counts across all samples, where n = # samples
low_count_filter <- rowSums(data) > ncol(data)
sum(low_count_filter) # genes that 'passed the check
sprintf("Removing %d low-count genes (%d remaining out of %d).", sum(!low_count_filter), sum(low_count_filter), nrow(data))
data_high <- data[low_count_filter,]

# Log2 transform the data
log_counts <- log2(data_high + 1)

# Only keep genes with variation greater than 1
initial_genes <- nrow(log_counts)
log_counts <- log_counts[apply(log_counts, 1, var) > 1,]
log_counts <- as.data.frame(log_counts)
kept_genes <- nrow(log_counts)
sprintf("Removing genes with variation less than 1 (%d kept out of %d total)", kept_genes, initial_genes) 

# Looking for correlations to NAGK
NAGK <- "218231_at" # NAGK probe name
t(log_counts[NAGK,]) # Average NAGK expression in all the available cell types

# Plot all data as a heatmap
heatmap.2(cor(log_counts), trace='none', main='Sample correlations (log2-transformed)')

# Calculate a pairwise pearson correlation matrix
cor_matrix <- cor(t(log_counts),method = "pearson")
cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$names <- rownames(cor_matrix)
rownames(cor_matrix) <- NULL

# Look at it
cor_matrix[1:5,1:5]
dim(cor_matrix)

# Order the data based on the similarity to NAGK
corr_top <- data.frame(Names = cor_matrix$names, Similarity_to_NAGK = cor_matrix[,NAGK])
# Order from highest to lowest
corr_top <- corr_top[order(corr_top$Similarity_to_NAGK, decreasing = T),] 
# Look at the data
head(corr_top)
tail(corr_top)

# Only keep the 0.7 and greater and remove nagk itself (r=1)
corr_top <- corr_top %>% filter(Similarity_to_NAGK > 0.7 & Similarity_to_NAGK < 1) 
rownames(corr_top) <- NULL
dim(corr_top)
sprintf("Transcripts with r > 0.7 to NAGK %d", nrow(corr_top))

# We end up with a relatively small list of genes, but to see what they are we need to translate the probes to gene names
# I used an online tool to do the conversion of top 300 genes, let's read it in now
gene_names <- read.delim("/Home Office/Data/2020.12.10.Che project/Data/BioGPS data/Translated Gene names.tsv")
# add a symbol column to the data
corr_top$Symbol <- gene_names$Symbol[match(corr_top$Names, gene_names$Affy)]
head(corr_top)

# Let's save a csv of them
write.csv(corr_top,"C:/Home Office/Data/2020.12.10.Che project/For_Paper/top_113_similar_genes.csv")

# we can run gene ontology on those genes using the panther online tool
# Lets read in the GO results table
go_data <- read.delim("/Home Office/Data/2020.12.10.Che project/Data/BioGPS data/2020.11.08_full_go_top_113_genes.tsv")
head(go_data)

############################################################
######################## PLOTTING ##########################
############################################################

# Keep all hits with FDR less than 10^-7
go_data <- go_data %>% filter(FDR <10^-6)
# Add an ordering column
go_data$Order <- nrow(go_data):1
# Add -log 10 of FDR 
go_data$log10fdr <- -log10(go_data$FDR)
# Make description look nicer
go_data$Description <- capitalize(gsub(" \\(.*", "", go_data$GO))  
# Make Enrichment numeric
go_data$Enrichment <- as.numeric(go_data$Enrichment)  

ggplot(go_data, aes(x = Order, y = Enrichment, color = log10fdr, size = Genes_found)) +
  geom_point(stat = "identity") +
  geom_hline(yintercept = 0, linetype='dashed') +
  coord_flip(ylim = c(0, 8)) +
  scale_color_continuous(low="#5D00A9", high="red") + 
  ylab("Enrichment") +
  scale_x_continuous(name = "", label = go_data$Description, breaks = go_data$Order) +
  # scale_y_continuous(breaks = 0:8) +
  theme(axis.text.y = element_text(color="#000000")) +
  ggtitle("TOP 20 GO Hits") + labs(size="Number of Genes", color="-log10(FDR)") +
  theme_classic()

# To save the plot in a more professional and consistent manner
ggsave("/Home Office/Data/2020.12.10.Che project/For_Paper/GO Plot.pdf", height = 10, width = 7)

# We can also make a heatmap of all genes in the "granulocyte activation" group
# Read in the gene names
Gran_genes <- read.delim("/Home Office/Data/2020.12.10.Che project/Data/BioGPS data/granulocyte activation (GO0036230).tsv")
Gran_genes <- Gran_genes$Symbol
# Translate to probeID
Gran_genes_affy <- gene_names$Affy[gene_names$Symbol %in% Gran_genes]

# keep the subset of data
heat_data=log_counts[Gran_genes_affy,]
heat_data=as.data.frame(t(heat_data))
# Rename back to symbols
colnames(heat_data)=gene_names$Symbol[match(colnames(heat_data), gene_names$Affy)]
# pheatmap(heat_data,clustering_distance_cols = "correlation",main = "Genes(columns) clustered by correlation",cluster_rows = F)
# Try re-clustering the columns
mat_cluster_cols <- hclust(dist(cor(heat_data)))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
# use better names
rownames(heat_data)=capitalize(gsub("^X", "", gsub("\\.", "", gsub("_", " ", rownames(heat_data)))))
# Plot
pheatmap(heat_data,cluster_cols= mat_cluster_cols,main = "Genes(columns) clustered by correlation",cluster_rows = F)

# To save the plot in a more professional and consistent manner
graphics.off()
pdf("/Home Office/Data/2020.12.10.Che project/For_Paper/Heatmap_new.pdf", height = 20, width = 8) 
pheatmap(heat_data,cluster_cols= mat_cluster_cols,main = "Genes(columns) clustered by correlation",cluster_rows = F)
dev.off()
graphics.off()

# Also do pairwise correlation plots between top 4 genes from granulocyte act. group
GoIs <- c("CECR1", "ARPC5", "MAN2B1", "COTL1")
AoIs <- c(gene_names$Affy[gene_names$Symbol%in%GoIs], "218231_at")

# Subset data again
corr_data <- log_counts[AoIs,]
corr_data <- as.data.frame(t(corr_data))
colnames(corr_data) <- gene_names$Symbol[match(colnames(corr_data), gene_names$Affy)]

corr_data <- as.data.frame(corr_data)

p1 <- ggscatter(corr_data, x = "NAGK", y = "COTL1", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "NAGK", ylab = "COTL1",title = "Pearson correlation")

p2 <- ggscatter(corr_data, x = "NAGK", y = "MAN2B1", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "NAGK", ylab = "MAN2B1",title = "Pearson correlation")

p3 <- ggscatter(corr_data, x = "NAGK", y = "CECR1", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "NAGK", ylab = "CECR1",title = "Pearson correlation")

p4 <- ggscatter(corr_data, x = "NAGK", y = "ARPC5", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "NAGK", ylab = "ARPC5",title = "Pearson correlation")

plot_grid(p1, p2, p3, p4, align = "hv")

# To save the plot in a more professional and consistent manner
ggsave("/Home Office/Data/2020.12.10.Che project/For_Paper/Correlation_plots_new.pdf", 
       height = 6, width = 6)
