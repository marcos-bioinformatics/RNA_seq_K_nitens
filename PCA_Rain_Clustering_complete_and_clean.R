## Authors: Francisco J. Romero-Campero - fran@us.es
##          Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: May 2020

## Example script to extract gene expression from RNA-seq processed samples and
## perform a Principal Component Analysys (PCA), Rain and Hierarchical Clustering
## analysis over Klebsormidium nitens gene sets.

## Sample folders with .ctab files should be in the same directory as the script.
## In R studio Session -> Set Working Directory -> To Source File Location


# Load required libraries
library(ballgown)
library(genefilter)
library(VennDiagram)
library(robustbase)
library(FactoMineR)
library("factoextra")
library(clusterProfiler)
library(ggplot2)
library(dplyr)

# Experimental design
experimental.design <- read.csv("experimental_design_klebsormidium.csv")
experimental.design

# Read data from each sample, stored in .ctab files inside each sample folder
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)
bg.data
sampleNames(bg.data)

# Extract expression level for each gene
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)

colnames(gene.expression) <- c("zt1_1", "zt1_2", "zt1_3", "zt3_1", "zt3_2", "zt3_3", "zt5_1",
                               "zt5_2", "zt5_3", "zt7_1", "zt7_2", "zt7_3", "zt9_1", "zt9_2",
                               "zt9_3", "zt11_1", "zt11_2", "zt11_3", "zt13_1", "zt13_2",
                               "zt13_3", "zt15_1", "zt15_2", "zt15_3", "zt17_1", "zt17_2",
                               "zt17_3", "zt19_1", "zt19_2", "zt19_3", "zt21_1", "zt21_2",
                               "zt21_3", "zt23_1", "zt23_2", "zt23_3", "zt25_1", "zt25_2",
                               "zt25_3")

# Save gene expression matrix
write.table(x=gene.expression,file = "gene_expression_complete.tsv",sep = "\t",quote = F)

# Calculate mean for triplicates to construct a mean expression matrix
zt1.mean <- (gene.expression[,"zt1_1"] + gene.expression[,"zt1_2"] + gene.expression[,"zt1_3"])/3
zt3.mean <- (gene.expression[,"zt3_1"] + gene.expression[,"zt3_2"] + gene.expression[,"zt3_3"])/3
zt5.mean <- (gene.expression[,"zt5_1"] + gene.expression[,"zt5_2"] + gene.expression[,"zt5_3"])/3
zt7.mean <- (gene.expression[,"zt7_1"] + gene.expression[,"zt7_2"] + gene.expression[,"zt7_3"])/3
zt9.mean <- (gene.expression[,"zt9_1"] + gene.expression[,"zt9_2"] + gene.expression[,"zt9_3"])/3
zt11.mean <- (gene.expression[,"zt11_1"] + gene.expression[,"zt11_2"] + gene.expression[,"zt11_3"])/3
zt13.mean <- (gene.expression[,"zt13_1"] + gene.expression[,"zt13_2"] + gene.expression[,"zt13_3"])/3
zt15.mean <- (gene.expression[,"zt15_1"] + gene.expression[,"zt15_2"] + gene.expression[,"zt15_3"])/3
zt17.mean <- (gene.expression[,"zt17_1"] + gene.expression[,"zt17_2"] + gene.expression[,"zt17_3"])/3
zt19.mean <- (gene.expression[,"zt19_1"] + gene.expression[,"zt19_2"] + gene.expression[,"zt19_3"])/3
zt21.mean <- (gene.expression[,"zt21_1"] + gene.expression[,"zt21_2"] + gene.expression[,"zt21_3"])/3
zt23.mean <- (gene.expression[,"zt23_1"] + gene.expression[,"zt23_2"] + gene.expression[,"zt23_3"])/3
zt25.mean <- (gene.expression[,"zt25_1"] + gene.expression[,"zt25_2"] + gene.expression[,"zt25_3"])/3 

# Create mean expression matrix
mean.expression <- matrix(data = c(zt1.mean, zt3.mean, zt5.mean, zt7.mean, zt9.mean, zt11.mean,
                                   zt13.mean, zt15.mean, zt17.mean, zt19.mean, zt21.mean, zt23.mean,
                                   zt25.mean),ncol=13)

colnames(mean.expression) <- c("ZT1","ZT3","ZT5","ZT7","ZT9","ZT11","ZT13","ZT15","ZT17",
                               "ZT19","ZT21","ZT23","ZT25")

rownames(mean.expression) <- rownames(gene.expression)

################################ Principal Components Analysis ################################ 

# Create transpose of gene expression matrix as a dataframe. This will be the PCA input
pca.gene.expression <- data.frame(colnames(gene.expression),t(gene.expression))
colnames(pca.gene.expression)[1] <- "Time point"

# Apply PCA to this dataframe
res.pca <- PCA(pca.gene.expression, graph = FALSE, scale.unit = TRUE, quali.sup = 1 )

# Show the eigenvalues and the amount of variation retained by each principal component
get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70),main = "")

# Plot samples over PC1 and PC2
pdf("eigenvalues.pdf")
fviz_pca_ind(res.pca, col.ind = experimental.design$sample, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="Principal Components Analysis",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# Contributions of variables to PC1
print(fviz_contrib(res.pca, choice = "var", axes = 1, top = 20))
# Contributions of variables to PC2
print(fviz_contrib(res.pca, choice = "var", axes = 2, top = 20))
# Contributions of samples to PC1
print(fviz_contrib(res.pca, choice = "ind", axes = 1))
# Contributions of samples to PC2
print(fviz_contrib(res.pca, choice = "ind", axes = 2))


################################# Hierarchical Clustering pre-rain #################################

# Compute Hierarchical Clustering on Principal Components
res.hcpc <- HCPC(res.pca, graph=FALSE)    

# Show cluster dendrogram
pdf("cluster.pdf")
fviz_dend(res.hcpc,k=3,
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle"#,
          # labels_track_height = 0.5    # Augment the room for labels
)
dev.off()

# Show sample clustering map
pdf("cluster_ind.pdf")
fviz_cluster(res.hcpc,ellipse = TRUE,
             repel = TRUE,               # Avoid label overlapping
             show.clust.cent = TRUE,     # Show cluster centers
             palette = "jco",            # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Sample Clustering"
)
dev.off()

# Show cluster dendrogram with scale
pdf("hierarchical_clustering_all_genotypes_three_clusters.pdf")
plot(res.hcpc, choice ="tree")
dev.off()



######################################### RAIN #########################################

# Rain provides a non-parametric method for detection of rhythms of pre-specified periods in 
# biological data. It uses hipothesis testing rather than measured values regarding waves.

library(rain)

# We have one series per gene, and for each gene 39 measurements at 13 different time points (nr.series=3)

rain.gene.expression <- gene.expression

# Determine cycling genes in LD
results.rain <- rain(t(rain.gene.expression), #one column per series, one row per time point
                     deltat=2,      #time difference between two data points
                     period=24,     #period to test for
                     verbose=TRUE, 
                     nr.series=3)   #different possibilities to specify multiple experiments

head(results.rain)
sum(results.rain$pVal < 0.05)/nrow(results.rain) ## 64,69%
sum(results.rain$pVal < 0.01)/nrow(results.rain) ## 62,21%
nrow(results.rain) # 17290

##AQUI FALTA UN PLOT DE UNO O DOS GENES PARA COMPROBAR QUE HAY RITMO CIRCADIANO

# Save a subset of genes considered circadian by rain for a p-value<0.01
circadian.genes.rain.p <- subset(results.rain, results.rain$pVal < 0.01)
circadian.genes.rain.p <- rownames(circadian.genes.rain.p)
length(circadian.genes.rain.p)

write(x = circadian.genes.rain.p, file = "circadian_genes_rain_p.txt")

######################### Principal components analysis post-RAIN #########################

# Subset the gene expression matrix with circadian genes by rain's vector
gene.expression.rain <- gene.expression[circadian.genes.rain.p, ]
dim(gene.expression.rain)

# Create transpose of rain expression matrix as a dataframe. This will be the PCA input
pca.gene.expression.circadian <- data.frame(colnames(gene.expression.rain),t(gene.expression.rain))
colnames(pca.gene.expression.circadian)[1] <- "Time point"

# Apply PCA to this dataframe
res.pca.circadian <- PCA(pca.gene.expression.circadian, graph = FALSE, scale.unit = TRUE, quali.sup = 1 )

# Show the eigenvalues and the amount of variation retained by each principal component
get_eigenvalue(res.pca.circadian)
fviz_eig(res.pca.circadian, addlabels = TRUE, ylim = c(0, 70),main = "")

# Plot samples over PC1 and PC2
pdf("eigenvalues_circadian.pdf")
fviz_pca_ind(res.pca.circadian, col.ind = experimental.design$sample, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="Principal Components Analysis",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# Contributions of variables to PC1
print(fviz_contrib(res.pca.circadian, choice = "var", axes = 1, top = 20))
# Contributions of variables to PC2
print(fviz_contrib(res.pca.circadian, choice = "var", axes = 2, top = 20))
# Contributions of samples to PC1
print(fviz_contrib(res.pca.circadian, choice = "ind", axes = 1))
# Contributions of samples to PC2
print(fviz_contrib(res.pca.circadian, choice = "ind", axes = 2))

################################# Hierarchical Clustering post-rain #################################

# Compute Hierarchical Clustering on Principal Components
res.hcpc.circadian <- HCPC(res.pca.circadian, graph=FALSE)    

# Show cluster dendrogram
pdf("cluster_circadian.pdf")
fviz_dend(res.hcpc.circadian,k=3,
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle"#,
          # labels_track_height = 0.5    # Augment the room for labels
)
dev.off()

# Show sample clustering map
pdf("cluster_ind_circadian.pdf")
fviz_cluster(res.hcpc.circadian,ellipse = TRUE,
             repel = TRUE,               # Avoid label overlapping
             show.clust.cent = TRUE,     # Show cluster centers
             palette = "jco",            # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Sample Clustering"
)
dev.off()

# Show cluster dendrogram with scale
pdf("hierarchical_clustering_all_genotypes_three_clusters_circadian.pdf")
plot(res.hcpc.circadian, choice ="tree")
dev.off()

#################################### RAIN Clustering ####################################

# Read circadian genes table and store it in a variable
circadian.genes.rain.p <- read.table(file="circadian_genes_rain_p.txt",as.is=T,header = F)[[1]]
head(circadian.genes.rain.p)

# Create and empty vector with 10754 entries which will be filled in the clustering
cluster.circadian.genes <- vector(mode="character", length=length(circadian.genes.rain.p))
cluster.circadian.genes
length(cluster.circadian.genes)

# Create mean expression matrix for circadian genes
mean.expression.rain <- mean.expression[circadian.genes.rain.p, ]
dim(mean.expression.rain)

# Create a loop to store peak and trough for each circadian gene
for (i in 1:length(circadian.genes.rain.p))
{
  print(i)
  expr.gene <- mean.expression.rain[circadian.genes.rain.p[i], ]
  trough <- substring(names(which.min(expr.gene)),3)
  trough <- as.numeric(strsplit(trough, "_") [[1]] [1])
  peak <- substring(names(which.max(expr.gene)),3)
  peak <- as.numeric(strsplit(peak, "_") [[1]] [1])
  cluster.circadian.genes[i] <- paste(c("peak", peak, "trough", trough), collapse="_")
}

# Create dataframe with each circadian gene + peak and trough points
clustering.rain.genes <- data.frame(circadian.genes.rain.p, cluster.circadian.genes)
colnames(clustering.rain.genes) <- c("Genes", "Cluster")
head(clustering.rain.genes)

write.table(clustering.rain.genes, file = "clustering_rain_genes_p.tsv", quote = FALSE, 
            row.names = FALSE, sep = "\t")

# Do the same clustering but just for the peak, ignoring trough
for (i in 1:length(circadian.genes.rain.p))
{
  print(i)
  expr.gene <- mean.expression.rain[circadian.genes.rain.p[i], ]
  peak <- substring(names(which.max(expr.gene)),3)
  peak <- as.numeric(strsplit(peak, "_") [[1]] [1])
  cluster.circadian.genes[i] <- paste("peak", peak, collapse="_")
}

# Create dataframe with each circadian gene + peak point
clustering.rain.genes.peak <- data.frame(circadian.genes.rain.p, cluster.circadian.genes)
colnames(clustering.rain.genes) <- c("Genes", "Cluster")
head(clustering.rain.genes.peak)

write.table(clustering.rain.genes.peak, file = "clustering_rain_genes_p_peak.tsv", quote = FALSE, 
            row.names = FALSE, sep = "\t")


