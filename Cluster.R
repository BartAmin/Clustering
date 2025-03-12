library(clValid)
library(tidyr)
library(opticskxi)
library(fpc)
library(mclust)
library(dbscan)
library(ggplot2)
library(cluster)
library(farff)
library(factoextra)  

# ===================== DATASET SELECTION =====================
#Select dataset by uncommenting -#, dataset 1 is default.

#Dataset 1
url <- "https://raw.githubusercontent.com/milaan9/Clustering-Datasets/master/02.%20Synthetic/2d-20c-no0.arff"

#dataset 2
#url <-https://raw.githubusercontent.com/milaan9/Clustering-Datasets/master/02.%20Synthetic/2d-20c-no0.arff

#dataset 3
# URL of the ARFF file
#url <- "https://raw.githubusercontent.com/milaan9/Clustering-Datasets/master/02.%20Synthetic/3MC.arff"

#dataset 4
#url <- "https://raw.githubusercontent.com/milaan9/Clustering-Datasets/master/02.%20Synthetic/cure-t2-4k.arff"

#Comment lines 29/30 if built-in datasets from cluster are used
temp_file <- tempfile()
download.file(url, temp_file, mode = "wb")  

#All datasets above are 2-dimensional
data <- readARFF(temp_file)[, 1:2] # Read ARFF file

#Dataset 5
#data('moons')
#data <- moons

#Dataset 6
#data('DS3')
#data <- DS3

# ===================== DATA VISUALIZATION =====================
plot(data, pch = 20, main = "Dataset Visualization")


# ===================== PREPARE DISSIMILARITY MATRIX =====================
dissimilarity_matrix <- dist(data)

# Initialize indices dataframe for cluster evaluation
indices <- data.frame(method = c('K-MEANS', 'HAC', 'GMM', 'HDBSCAN', 'OPTICS'))

# ===================== FIT CLUSTERING MODELS =====================

# K-MEANS-CLUSTERING
# Find optimal K (should be at the elbow), in dataset 1 this would be 4
fviz_nbclust(data, FUN = kmeans, method = "wss")  
cluster_kmeans <- kmeans(data, centers = 4, nstart = 50)$cluster


# AGGLOMERATIVE HIERARCHICAL CLUSTERING (AHC)
# Cut the dendrogram where the height increases the most before the 
# clusters merge. For dataset 1 this could be at K = 3 or K = 6
hclust_result <- hclust(dissimilarity_matrix, method = 'average')
plot(hclust_result, main = "Dendrogram of AHC",
     sub = "", hang = -1, cex = 0.8)
abline(h = 8, col = 'red')
cluster_hclust <- cutree(hclust_result, k = 6)

# GAUSSIAN MIXTURE MODELS
cluster_gmm <- Mclust(data)$classification

# HDBSCAN CLUSTERING
# minPts chosen somewhat arbitrarily, can be done more formally with hyperparameter tuning
hdbscan <- hdbscan(data, minPts = 10) 
plot(hdbscan)
cluster_hdbscan = hdbscan$cluster
noise_index_hdbscan <- which(cluster_hdbscan == 0) #store noise indices separately

# OPTICS CLUSTERING
# minPts chosen somewhat arbitrarily, can be done more formally with hyperparameter tuning.
cluster_optics_rp <- optics(data, minPts = 10)

# Extract clusters using DBSCAN approach from OPTICS results, the red line should 
# be just above all the valleys and under the peaks (0.8 for dataset 1)
plot(cluster_optics_rp)
abline(h = 0.8, col = "red")
cluster_optics <- extractDBSCAN(cluster_optics_rp, eps_cl = 0.8)$cluster
noise_index_optics <- which(cluster_optics == 0) #store noise indices separately

# ===================== CLUSTER VALIDATION METRICS =====================

# Function to compute cluster statistics after removing noise points
cluster_stats <- function(dmat, cluster_labels, noise_index, metric) {
  if (length(noise_index) > 0) {
    valid_dmat <- as.dist(as.matrix(dmat)[-noise_index, -noise_index])  
    valid_clusters <- cluster_labels[-noise_index]  
  } else {
    valid_dmat <- dmat
    valid_clusters <- cluster_labels
  }
  
  # Compute the requested metric 
  stats <- cluster.stats(valid_dmat, valid_clusters)
  return(stats[[metric]])  # Extract specific metric
}

# Compute Dunn Index 
# The higher the value the better the cluster quality
indices$dunn2 <- c(
  cluster.stats(dissimilarity_matrix, cluster_kmeans)$dunn2,
  cluster.stats(dissimilarity_matrix, cluster_hclust)$dunn2,
  cluster.stats(dissimilarity_matrix, cluster_gmm)$dunn2,
  cluster_stats(dissimilarity_matrix, cluster_hdbscan, noise_index_hdbscan, "dunn2"),
  cluster_stats(dissimilarity_matrix, cluster_optics, noise_index_optics, "dunn2")
)

# Compute Within-Between Ratios, convert to Between-Within
# The higher the value the better the cluster quality
indices$between_within <- c(
  1 / cluster.stats(dissimilarity_matrix, cluster_kmeans)$wb.ratio,
  1 / cluster.stats(dissimilarity_matrix, cluster_hclust)$wb.ratio,
  1 / cluster.stats(dissimilarity_matrix, cluster_gmm)$wb.ratio,
  1 / cluster_stats(dissimilarity_matrix, cluster_hdbscan, noise_index_hdbscan, "wb.ratio"),
  1 / cluster_stats(dissimilarity_matrix, cluster_optics, noise_index_optics, "wb.ratio")
)

# ===================== PLOT CLUSTER VALIDATION METRICS =====================
pdf("clustering_results.pdf", width = 10, height = 8)  

# Define colors
bar_colors <- c("firebrick1", "dodgerblue")

# Create matrix for bar heights
values_matrix <- rbind(indices$dunn2,indices$between_within)

# Create bar chart
barplot(values_matrix, beside = TRUE, col = bar_colors, names.arg = indices$method, 
        main = "Clustering Performance Metrics", cex.main = 1.5, 
        ylab = "Score", cex.lab = 1.2, cex.axis = 1.1, 
        ylim = c(0, max(values_matrix) * 1.2))  

# Add legend
legend("topleft", legend = c("Dunn Index", "Between-Within Ratio"), 
       fill = bar_colors, bty = "n", cex = 1.1)

# ===================== PLOT CLUSTER RESULTS =====================
# Function to ensure only noise points (cluster 0) are black in plotting
plot_clusters <- function(data, clusters, title) {
  unique_clusters <- sort(unique(clusters))
  
  # Assign colors ensuring noise (0) is black
  colors <- rep("black", length(unique_clusters))  
  valid_clusters <- unique_clusters[unique_clusters != 0]  # Exclude noise
  
  # Assign distinct colors only for valid clusters
  if (length(valid_clusters) > 0) {
    colors[unique_clusters != 0] <- rainbow(length(valid_clusters))
  }
  
  # Map cluster labels to colors
  cluster_colors <- colors[match(clusters, unique_clusters)]  
  
  plot(data, col = cluster_colors, main = title, pch = 20, cex = 1.5, 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")
}

# Apply to each clustering method
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
plot(data, pch = 20, main = "Dataset Visualization")
plot_clusters(data, cluster_kmeans, "K-Means Clustering")
plot_clusters(data, cluster_hclust, "HAC Clustering")
plot_clusters(data, cluster_gmm, "GMM Clustering")
plot_clusters(data, cluster_hdbscan, "HDBSCAN Clustering")
plot_clusters(data, cluster_optics, "OPTICS Clustering")

dev.off() 
