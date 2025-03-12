# Clustering Methods in R

This repository contains an R script for comparing various clustering algorithms to different datasets. The script allows users to select a dataset by uncommenting the desired URL.
It provides visualization tools and computes cluster quality metrics to assess clustering performance. 
If you integrate your own data and it has more than two dimensions, the plotting function will no longer work. To address this, you can either apply Singular Value Decomposition (SVD) (or PCA) to reduce the data to a rank-2 matrix, or simply skip the plotting function.

# The script fits multiple clustering methods, categorized as follows:

- **Centroid-based clustering**: K-Means  
- **Hierarchical-based clustering**: Agglomerative Hierarchical Clustering (AHC)  
- **Distribution-based clustering**: Gaussian Mixture Models (GMM)  
- **Density-based clustering**: OPTICS & HDBSCAN

# Choosing optimal number of clusters

Clustering algorithms whereby it is necessary to determine the number of clusters beforehand have a comment with a suggestion on how to do so. 

# Further reading

For a more comprehensive description on the (dis)-advantages of each algorithm, and the rationale behind the methods of choosing the number of clusters, see: https://bartamin.com/advanced_clustering/
