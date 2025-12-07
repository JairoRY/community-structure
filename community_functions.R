library(igraph)
library(clustAnalytics)

data(karate, package="igraphdata")
wc <- walktrap.community(karate)
clust1 <- unname(membership(wc))
plot(wc, karate)

fc <- fastgreedy.community(karate)
clust2 <- unname(membership(fc))
dendPlot(fc)

# FUNCTION TASK 1
task_1 <- function(memb1, memb2) {
  # Get the unique cluster labels for each clustering
  clusters1 <- unique(memb1)
  clusters2 <- unique(memb2)
  
  # Initialize the result matrix
  js_matrix <- matrix(NA, nrow = length(clusters1), ncol = length(clusters2),
                      dimnames = list(paste0("C1_", clusters1), paste0("C2_", clusters2)))
  
  # Iterate through all pairs of clusters and compute the Jaccard Index
  for (i in 1:length(clusters1)) {
    # Get the set of nodes belonging to the current cluster in memb1
    cluster_A_label <- clusters1[i]
    cluster_A_nodes <- which(memb1 == cluster_A_label)
    
    for (j in 1:length(clusters2)) {
      # Get the set of nodes belonging to the current cluster in memb2
      cluster_B_label <- clusters2[j]
      cluster_B_nodes <- which(memb2 == cluster_B_label)
      
      # A is cluster_A_nodes, B is cluster_B_nodes
      
      # Calculate intersection size
      intersection_size <- length(intersect(cluster_A_nodes, cluster_B_nodes))
      
      # Calculate union size
      union_size <- length(union(cluster_A_nodes, cluster_B_nodes))
      
      # Calculate Jaccard Index
      if (union_size == 0) {
        jaccard_index <- 0
      } else {
        jaccard_index <- intersection_size / union_size
      }
      
      # Store the result
      js_matrix[i, j] <- jaccard_index
    }
  }
  return(js_matrix)
}

task_1(clust1, clust2)
