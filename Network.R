## Authors: Francisco J. Romero-Campero - fran@us.es
##          Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: May 2020

## Example script to perform Network Construction and Network analysis over 
## Klebsormidium nitens gene expression matrix.

############################## Mean expression Rain Network construction ##############################

# Load the igraph library
library(igraph)

# Create the correlation matrix
rain.correlation <- cor(t(mean.expression.rain))
dim(rain.correlation)

# Create a threshold list to iterate and two empty vectors to store connectivity and R squared values
thresholds <- seq(from=0.85,to=0.99,by=0.005) 
mean.connectivities <- vector(length=length(thresholds))
scale.free.R2 <- vector(length=length(thresholds))

# Go through all the correlation thresholds
for(i in 1:length(thresholds))
{
  print(thresholds[i])
  # Construct the network corresponding to the current threshold being checked
  
  # Adjacency matrix
  current.adjacency <- (rain.correlation > thresholds[i] & rain.correlation < 1) #Conditions
  # Network
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  # Compute node degrees
  node.degrees <- degree(threshold.network)
  
  # Keep track of the mean connectivity or mean node degree of the current network
  mean.connectivities[i] <- mean(node.degrees)
  
  # Check scale free property graphically
  h <- hist(node.degrees)
  
  # Compute degree frequencies
  degree.frequencies <- table(node.degrees)
  
  # Determine linear regression for logarithmic transformed degree frequencies
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  
  # Extract R squared as a measure of adjustment to the scale free property
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Correlation Threshold",ylab="Mean connectivity")
plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.70,0.99),xlab="Correlation Threshold",ylab="Scale Free Model Fit (R²)")

# Name and show mean connectivities table and R squared adjust to free-scale property table
names(mean.connectivities) <- thresholds
names(scale.free.R2) <- thresholds

mean.connectivities
scale.free.R2 

# Choose the appropiate threshold for network construction
adjacency.0975 <- (rain.correlation > 0.975) & (rain.correlation < 1)
gene.coexpression.network <- graph.adjacency(adjacency.0975, mode="undirected")
write.graph(gene.coexpression.network,file="klebso_gene_coexpression_mean_expression_rain_0975.gml",format="gml")

gene.coexpression.network <- read.graph(file="klebso_gene_coexpression_mean_expression_rain_0975.gml",format="gml")


################################ Network Analysis ################################

# Apply power law to network degree distribution to check free-scale property with KS test
network.degree.distribution <- degree.distribution(gene.coexpression.network)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]] #p-value of Kolmogorov-Smirnov test.

# Calculate degree for each node of the network
network.degrees <- degree(gene.coexpression.network)
network.degrees
degree.histogram <- hist(network.degrees,freq=FALSE,col="blue",xlab="Node degree", ylab="Probability",main="Degree distribution")

# Calculate absolute frequency for degree
degree.frequencies <- table(network.degrees)

# Eliminate nodes with a degree value equal to zero
degree.frequencies.no.0 <- degree.frequencies[-1]

# Subset those nodes over the 95th percentile
quantile.095 <- quantile(network.degrees,prob=0.95)
degree.values <- network.degrees[network.degrees > quantile.095]
degree.names <- names(degree.values)
length(degree.values) #534 Genes
write.table(degree.values,file="network_degrees", quote = FALSE, col.names = FALSE)

# Logaritmic transformation
log10.degrees.frequencies <- log10(degree.frequencies.no.0)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies.no.0)))

# Run a linear regression
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r)
plot(lm.r)

# Hub identifying 
network.hub.scores <- hub.score(gene.coexpression.network)
hub.score.attributes <- network.hub.scores[["vector"]]

write.table(hub.score.attributes,file="hub_score_attributes.txt", 
            quote = FALSE, #No quotes if you want to load it on Cytoscape
            col.names = FALSE) 

# Subset those nodes over the 95th percentile
quantile.095 <- quantile(hub.score.attributes,prob=0.95)
hubs.values <- hub.score.attributes[hub.score.attributes > quantile.095]
hubs.names <- names(hubs.values)
length(hubs.values) #538 Genes
write.table(hubs.names,file="network_hubs", quote = FALSE, col.names = FALSE)

# Calculate each node's transitivity and store it in a text file
node.transitivity <- transitivity(gene.coexpression.network,type = "local",isolates = "zero")
names(node.transitivity) <- names(V(gene.coexpression.network))
write.table(node.transitivity,file="node_transitivity.txt", quote = FALSE, col.names = FALSE)

# Subset those nodes with a transitivity value over 95th percentile
quantile.095 <- quantile(node.transitivity,prob=0.95)
transitivity.values <- node.transitivity[node.transitivity > quantile.095]
transitivity.names <- names(node.transitivity)
length(transitivity.values) #530 Genes
write.table(transitivity.values,file="network_transitivity", quote = FALSE, col.names = FALSE)

# Calculate network's clustering coefficient
network.clustering.coefficient <- transitivity(gene.coexpression.network,type="global") #0.4111 

# Calculate network's average path length
average.path.length <- average.path.length(gene.coexpression.network, directed=FALSE) #11.961

############################## Montecarlo Method for Network simulation ##############################

# Structure of the network
gene.coexpression.network #10754 Nodos y 67135 Aristas

## Vamos a crear una red y posteriormente un conjunto de 1000 redes aleatorias para compararla
## con cualquier característica de nuestra red que queramos.

## En cada paso se añade un nuevo nodo por lo tanto en promedio en cada paso deberíamos
## añadir el siguiente número de aristas numero_de_aristas/numero_de_nodos

# We got to create the following number of edges each step
67135/10754
## Ya que se deben añadir un número entero de aristas redondeamos a la baja, añadimos todos los nodos
## menos el último y en el último paso al introducir el último nodo añadimos le resto de aristas.
67135/6

67135 - 10753*6
sum(c(rep(6,10753),2617))


# Generate random graphs and calculate clustering coefficient
number.of.added.edges <- c(rep(6,10753),2617) 
random.scale.free.graph <- barabasi.game(n=10754, #Número de nodos
                                         out.seq=number.of.added.edges, #Out-degree sequence of the vertices
                                         directed=FALSE) #Type of network
transitivity(random.scale.free.graph)

# Create two vectors to store clustering coefficients and average path length for each network
clustering.coefficients <- vector(length=1000) 
average.path.length.network <- vector(length=1000)

# Simulate 1000 networks and check for small world property
for(i in 1:1000)
{
  print(i)
  random.scale.free.graph <- barabasi.game(n=10754,out.seq=number.of.added.edges,directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
  #average.path.length.network[i] <- average.path.length(random.scale.free.graph) #Takes long
}

# Seleccionar aquellas redes que tengan un coeficiente de clustering mayor que el de nuestra red
# Ver la proporción de ellas que cumplen este criterio

# Check how many networks have a clustering coefficient higher than our network
# Check how many networks have an average path lenght longer than our network
sum(clustering.coefficients > network.clustering.coefficient) / 1000 
sum(average.path.length.network > average.path.length) / 1000 


############################## Network Clustering ##############################

# For clustering co-expressed genes, load WGCNA and cluster libraries
library("WGCNA")
library("cluster")
allowWGCNAThreads()

# For cluster identifying, we select 1 - correlation as a measure of similitude
similarity.matrix <- 1 - rain.correlation

# Use the hclust function on similarity matrix
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

# Cutree function allows us to cut the tree for a k=number of desired clusters
hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)

# PAM function determines clustering basec on centroid method
pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(similarity.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(similarity.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(similarity.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(similarity.matrix),k=10,diss=TRUE)

# Silhouette function calculates the clustering silhouette
sil2 <- silhouette(hclust.2,dist=similarity.matrix)
sil3 <- silhouette(hclust.3,dist=similarity.matrix)
sil4 <- silhouette(hclust.4,dist=similarity.matrix)
sil5 <- silhouette(hclust.5,dist=similarity.matrix)
sil6 <- silhouette(hclust.6,dist=similarity.matrix)
sil7 <- silhouette(hclust.7,dist=similarity.matrix)
sil8 <- silhouette(hclust.8,dist=similarity.matrix)
sil9 <- silhouette(hclust.9,dist=similarity.matrix)
sil10 <- silhouette(hclust.10,dist=similarity.matrix)

plot(sil10,border="blue")

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])
names(hclust.sil.values) <- c("k=2", "k=3", "k=4", "k=5", "k=6", "k=7", "k=8", "k=9", "k=10")

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])
names(pam.sil.values) <- c("k=2", "k=3", "k=4", "k=5", "k=6", "k=7", "k=8", "k=9", "k=10")

# Plot silhouette width vs. number of clusters for PAM and hclust methods
plot(2:10,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.1,1.0),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:10,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3,cex=0.70)

## Cluster visualization
hclust.2
write.table(hclust.2,file="hclust_2.txt",quote=FALSE,col.names=FALSE)

pam.2[["clustering"]]

clustering.pam.2 <- pam.2[["clustering"]]
write.table(clustering.pam.2,file="pam_2.txt",quote=FALSE,col.names=FALSE)

