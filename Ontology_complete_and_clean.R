## Authors: Francisco J. Romero-Campero - fran@us.es
##          Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: May 2020

## Example script to perform GO term and KEGG pathway enrichment
## analysis over Klebsormidium nitens gene sets

## Copy a text file with a single column with the ids of the genes/proteins
## you want to analyse. The ids should follow the lastest version genome
## annotation.

## In R studio Session -> Set Working Directory -> To Source File Location

# Install the package clusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# Install the annotation package for Klebsormidium nitens
install.packages("./org.Knitens.eg.db/", repos = NULL, type = "source")

# Load both libraries
library(clusterProfiler)
library(org.Knitens.eg.db)

# Load gene lists for PCA/Rain ZT
knitens.1 <- read.table(file = "klebs_go_1", header = FALSE,as.is = TRUE)[[1]]
length(knitens.1)
knitens.3 <- read.table(file = "klebs_go_3", header = FALSE,as.is = TRUE)[[1]]
length(knitens.3)
knitens.5 <- read.table(file = "klebs_go_5", header = FALSE,as.is = TRUE)[[1]]
length(knitens.5)
knitens.7 <- read.table(file = "klebs_go_7", header = FALSE,as.is = TRUE)[[1]]
length(knitens.7)
knitens.9 <- read.table(file = "klebs_go_9", header = FALSE,as.is = TRUE)[[1]]
length(knitens.9)
knitens.11 <- read.table(file = "klebs_go_11", header = FALSE,as.is = TRUE)[[1]]
length(knitens.11)
knitens.13 <- read.table(file = "klebs_go_13", header = FALSE,as.is = TRUE)[[1]]
length(knitens.13)
knitens.15 <- read.table(file = "klebs_go_15", header = FALSE,as.is = TRUE)[[1]]
length(knitens.15)
knitens.17 <- read.table(file = "klebs_go_17", header = FALSE,as.is = TRUE)[[1]]
length(knitens.17)
knitens.19 <- read.table(file = "klebs_go_19", header = FALSE,as.is = TRUE)[[1]]
length(knitens.19)
knitens.21 <- read.table(file = "klebs_go_21", header = FALSE,as.is = TRUE)[[1]]
length(knitens.21)
knitens.23 <- read.table(file = "klebs_go_23", header = FALSE,as.is = TRUE)[[1]]
length(knitens.23)
knitens.25 <- read.table(file = "klebs_go_25", header = FALSE,as.is = TRUE)[[1]]
length(knitens.25)
knitens.dawn <- read.table(file = "klebs_go_dawn", header = FALSE,as.is = TRUE)[[1]]
length(knitens.dawn)
knitens.day <- read.table(file = "klebs_go_day", header = FALSE,as.is = TRUE)[[1]]
length(knitens.day)
knitens.night <- read.table(file = "klebs_go_night", header = FALSE,as.is = TRUE)[[1]]
length(knitens.night)

# Load gene lists for network attributes
knitens.hubs <- read.table(file = "network_hubs", header = FALSE, as.is = TRUE)[[1]]
length(knitens.hubs)

knitens.transitivity <- read.table(file = "network_transitivity", header = FALSE, as.is = TRUE)[[1]]
length(knitens.transitivity)

knitens.degrees <- read.table(file="network_degrees", header = FALSE, as.is = TRUE)[[1]]
length(knitens.degrees)


# Load gene lists for clustering
knitens.hclust1 <- read.table(file="hclust_2_cluster1.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hclust1)
knitens.hclust2 <- read.table(file="hclust_2_cluster2.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hclust2)
knitens.pam1 <- read.table(file="pam_2_cluster1.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.pam1)
knitens.pam2 <- read.table(file="pam_2_cluster2.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.pam2)


# Load gene lists for hubs
knitens.hub1 <- read.table(file="hub_kfl00419_0110_NUMBER1.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub1)
knitens.hub2 <- read.table(file="hub_kfl00755_0040_NUMBER2.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub2)
knitens.hub3 <- read.table(file="hub_kfl00356_0110_NUMBER3.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub3)
knitens.hub4 <- read.table(file="hub_kfl00020_0250_NUMBER4.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub4)
knitens.hub5 <- read.table(file="hub_kfl00126_0270_NUMBER5.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub5)
knitens.hub6 <- read.table(file="hub_kfl00494_0135_NUMBER6.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub6)
knitens.hub7 <- read.table(file="hub_kfl00006_0430_NUMBER7.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub7)
knitens.hub8 <- read.table(file="hub_kfl00173_0150_NUMBER8.txt", header=FALSE, as.is=TRUE)[[1]]
length(knitens.hub8)

# Load gene lists for TFs of interest
knitens.B3 <- read.table(file="1B3_kfl00393_0120.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.Bzip1 <- read.table(file="1BZIP_kfl00031_0340.txt", header=FALSE, as.is=TRUE)[[1]] #### FOUND
knitens.Bzip2 <- read.table(file="1BZIP_kfl00245_0150.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.C2H2 <- read.table(file="1C2H2_kfl00433_0010.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.C3H1 <- read.table(file="1C3H_kfl00015_0110.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.C3H2 <- read.table(file="1C3H_kfl00226_0130.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.CAMTA <- read.table(file="1CAMTA_kfl00375_0090.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.DOF <- read.table(file="1DOF_kfl00762_0030.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.E2F <- read.table(file="1E2F_kfl00334_0090.txt", header=FALSE, as.is=TRUE)[[1]] #### FOUND
knitens.E2FDP <- read.table(file="1E2FDP_kfl00089_0200.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.GRF <- read.table(file="1GRF_kfl00186_0090.txt", header=FALSE, as.is=TRUE)[[1]] #### FOUND
knitens.HD_ZIP <- read.table(file="1HD-ZIP_kfl00396_0130.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.NFYB <- read.table(file="1NF_YB_kfl00013_0370.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.NFYC1 <- read.table(file="1NFYC_kfl00123_0030.txt", header=FALSE, as.is=TRUE)[[1]] #### FOUND
knitens.NFYC2 <- read.table(file="1NFYC_kfl00150_0190.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.WRKY <- read.table(file="1WRKY_kfl00096_0150.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.ZFHD <- read.table(file="1ZF-HD_kfl01106_0010.txt", header=FALSE, as.is=TRUE)[[1]] #### FOUND
knitens.CH33 <- read.table(file="1_CH3.txt", header=FALSE, as.is=TRUE)[[1]]
knitens.MYB <- read.table(file="1_MYB_kfl00026_0540.txt", header=FALSE, as.is=TRUE)[[1]]


################################## GO Analysis ##################################

# Set the background or universe. Here we use the entire annotated genome
knitens.universe <- unique(select(org.Knitens.eg.db, columns = c("GO"),keys=keys(org.Knitens.eg.db,keytype = "GID"))[["GID"]])
length(knitens.universe)

# Perform GO enrichment analysis with the function enrichGO and our annotation package
# over the previously loaded gene set example and the selected universe
ego <- enrichGO(gene          = knitens.MYB, #List of genes to analyze
                universe      = knitens.universe,
                OrgDb         = org.Knitens.eg.db,
                ont           = "BP", #Biological process
                pAdjustMethod = "BH", #Benjamini-Hochberg
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")

# Illustrate the results with different plots
barplot(ego,drop=TRUE,showCategory = 10)
goplot(ego)
dotplot(ego)
emapplot(ego)
cnetplot(ego)

# Generate a data frame with the results and save it into a text file
res.ego <- as.data.frame(ego)
head(res.ego)
write.table(x = res.ego,file = "go_enrichment_result_CH3_kfl00242_0140.tsv",sep = "\t")

# If you want to perform KEGG pathway enrichment analysis install and load the package pathview
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
library("pathview")

################################## KEGG Analysis ##################################

# Select manually the KO terms from knitens.universe
knitens.ko <- select(org.Knitens.eg.db,columns = c("KO"),keys=keys(org.Knitens.eg.db,keytype = "GID"))
ko.universe <- knitens.ko$KO
ko.universe <- ko.universe[!is.na(ko.universe)]

target.ko <- subset(knitens.ko,GID %in% knitens.pam2)$KO
target.ko <- target.ko[!is.na(target.ko)]

# Perform KEGG pathway enrichment analysis with the enrichKEGG function
pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, 
                                 organism = "ko", 
                                 universe = ko.universe,
                                 qvalueCutoff = 0.05))
head(pathway.enrichment)

# Generate a data frame with the results and save it into a text file
write.table(pathway.enrichment, "pathway_enrichment_pam2", quote = FALSE, row.names = FALSE, sep = "\t")


# Generate files with graphical representation of the enriched pathways
# Associate the value 1 with your gene set and 0 to with rest of genes in 
# universe/background
genes.pathway <- rep(0,length(ko.universe))
names(genes.pathway) <- ko.universe
genes.pathway[target.ko] <- 1

# Generate graphical representations 
for(i in 1:nrow(pathway.enrichment))
{
  pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id = pathway.enrichment$ID[i], 
           species = "ko",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")
}

# Similarly you can perform KEGG module enrichment analysis
mkk <- enrichMKEGG(gene = knitens.3, organism = "ko", keyType = "kegg")
head(mkk)

res.mkk <- as.data.frame(mkk)
write.table(x = res.mkk,file = "KEGG_module_enrichment_result_3.tsv",sep = "\t")

