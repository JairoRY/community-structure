library(igraph)
library(clustAnalytics)

# FUNCTION TASK 1
jaccard_sim <- function(memb1, memb2) {
  # Get the unique cluster labels for each clustering
  clusters1 <- unique(memb1)
  clusters2 <- unique(memb2)
  
  # Initialize the result matrix
  js_matrix <- matrix(NA, nrow = length(clusters1), ncol = length(clusters2),
                      dimnames = list(as.character(clusters1), as.character(clusters2)))
  
  # Iterate through all pairs of clusters and compute the Jaccard Index
  for (i in 1:length(clusters1)) {
    # Get the set of nodes belonging to the current cluster in memb1
    cluster_A_label <- clusters1[i]
    cluster_A_nodes <- which(memb1 == cluster_A_label)
    
    for (j in 1:length(clusters2)) {
      # Get the set of nodes belonging to the current cluster in memb2
      cluster_B_label <- clusters2[j]
      cluster_B_nodes <- which(memb2 == cluster_B_label)
      
      # Calculate intersection and union size
      intersection_size <- length(intersect(cluster_A_nodes, cluster_B_nodes))
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

# FUNCTION TASK 2
match_clusters <- function(js_matrix, name1 = "Clust1", name2 = "Clust2") {
  # Initialize vectors to store results
  best_jaccard_indices <- numeric(nrow(js_matrix))
  match_names <- character(nrow(js_matrix))
  
  # Iterate through each row (cluster from labeling 1)
  for (i in 1:nrow(js_matrix)) {
    # Find the maximum Jaccard index in the current row and its column
    max_index <- max(js_matrix[i, ])
    match_col_index <- which.max(js_matrix[i, ])
    best_jaccard_indices[i] <- max_index
    
    # Get the names for the matched pair
    cluster1_label <- rownames(js_matrix)[i]
    cluster2_label <- colnames(js_matrix)[match_col_index]
    
    # Store the paired cluster names in the format (name1.i, name2.j)
    match_names[i] <- paste0("(", name1, ".", cluster1_label, ",", name2, ".", cluster2_label, ")")
  }
  
  # Name the vector of best Jaccard indices and return it
  names(best_jaccard_indices) <- match_names
  return(best_jaccard_indices)
}

# FUNCTION TASK 3
Wmean <- function(memb1, MC) {
  # Get the cluster sizes and labels for the first clustering
  cluster_sizes <- table(memb1)
  cluster_labels <- names(cluster_sizes) # e.g., "1", "2", "3", "4"
  
  # Extract the cluster labels for the first clustering from the names of the MC vector in the format (Prefix.Label,...)
  mc_names_clean <- sub("^\\(", "", names(MC))
  mc_cluster1_names <- sub(",.*$", "", mc_names_clean)
  mc_labels <- sub(".*\\.", "", mc_cluster1_names)
  
  # Get the raw weights (cluster sizes) corresponding to the order in MC
  size_lookup <- setNames(as.numeric(cluster_sizes), cluster_labels)
  weights_raw <- size_lookup[mc_labels]
  
  # Compute the fractional weights
  total_nodes <- sum(cluster_sizes)
  weights <- weights_raw / total_nodes
  
  # Compute the weighted mean
  weighted_mean_jaccard <- sum(MC * weights)
  return(weighted_mean_jaccard)
}

# -----------------------------------------------------------------------------
# TASK 4: SIGNIFICANCE EVALUATION & JACCARD ANALYSIS

# Define Algorithms and Scoring Functions
alg_list <- list(
  Louvain = cluster_louvain,
  LabelProp = cluster_label_prop,
  Walktrap = cluster_walktrap,
  EdgeBetween = cluster_edge_betweenness
)

# Helper to run full analysis on a graph
analyze_network <- function(g, net_name, gt_memb = NULL, n_iters = 1) {
  # Significance Scoring
  base_scores <- evaluate_significance(g, alg_list = alg_list, gt_clustering = gt_memb)
  accum_scores <- as.matrix(base_scores)
  
  if (n_iters > 1) {
    for (k in 2:n_iters) {
      next_scores <- evaluate_significance(g, alg_list = alg_list, gt_clustering = gt_memb)
      accum_scores <- accum_scores + as.matrix(next_scores)
    }
  }
  # Compute average and select desired metrics
  avg_scores_matrix <- accum_scores / n_iters
  scores <- as.data.frame(avg_scores_matrix)
  row_id   <- grep("Internal density", rownames(scores), ignore.case = TRUE)
  row_exp  <- grep("Expansion", rownames(scores), ignore.case = TRUE)
  row_cond <- grep("Conductance", rownames(scores), ignore.case = TRUE)
  row_mod  <- grep("Modularity", rownames(scores), ignore.case = TRUE)
  
  rows_to_use <- c(row_id, row_exp, row_cond, row_mod)
  if(length(rows_to_use) > 0) {
    print(scores[rows_to_use, , drop=FALSE])
  } else {
    print(scores)
  }
  
  # Jaccard Analysis
  # Determine Reference Clustering
  if (!is.null(gt_memb)) {
    ref_memb <- gt_memb
    ref_name <- "GT"
  } else {
    # CASE No Ground Truth
    if (length(rows_to_use) < 4) {
      stop("Error: Could not find all 4 required metrics (Internal density, Expansion, Conductance, Modularity) in the output rows.")
    }
    
    # Extract numeric vectors for the algorithms
    algs_in_scores <- intersect(names(alg_list), colnames(scores))
    subset_scores <- scores[, algs_in_scores, drop=FALSE]
    
    # Extract values for the 4 metrics
    vals_id   <- as.numeric(subset_scores[row_id, ])
    vals_exp  <- as.numeric(subset_scores[row_exp, ])
    vals_cond <- as.numeric(subset_scores[row_cond, ])
    vals_mod  <- as.numeric(subset_scores[row_mod, ])
    
    # High is Best: Internal Density, Modularity -> Rank descending (-)
    rank_id  <- rank(-vals_id)
    rank_mod <- rank(-vals_mod)
    
    # Low is Best: Expansion, Conductance -> Rank ascending
    rank_exp  <- rank(vals_exp)
    rank_cond <- rank(vals_cond)
    
    # Compute Average Rank for each Algorithm
    rank_matrix <- rbind(rank_id, rank_exp, rank_cond, rank_mod)
    mean_ranks <- colMeans(rank_matrix)
    names(mean_ranks) <- algs_in_scores
    
    # Find Best Algorithm (Lowest Mean Rank)
    best_alg_name <- names(which.min(mean_ranks))
    cat("   Average Ranks (Lower is better):\n")
    print(mean_ranks)
    cat(paste0("-> Selected Reference: ", best_alg_name, "\n"))
    
    # Generate the reference membership
    ref_alg <- alg_list[[best_alg_name]](g)
    ref_memb <- as.numeric(membership(ref_alg))
    ref_name <- best_alg_name
  }
  
  # Compute Jaccard for all algorithms vs Reference
  for (alg_name in names(alg_list)) {
    comm <- alg_list[[alg_name]](g)
    alg_memb <- as.numeric(membership(comm))
    
    JS <- jaccard_sim(ref_memb, alg_memb)
    MC <- match_clusters(JS, name1 = ref_name, name2 = alg_name)
    global_sim <- Wmean(ref_memb, MC)
    
    # Print results
    cat(paste0("\n", ref_name, " vs ", alg_name, "\n"))
    cat(sprintf("Global Jaccard: %.4f\n", global_sim))
    cat("Local Matches:\n")
    print(MC)
  }
}

# RUNNING THE ANALYSIS ON THE 4 NETWORKS

# Network 1: Karate Club
data(karate, package="igraphdata")
# Ground truth is stored in "Faction" attribute
gt_karate <- as.numeric(as.factor(V(karate)$Faction))
analyze_network(karate, "Karate Club", gt_memb = gt_karate)


# Network 2: Synthetic Scale-Free
B_pref <- matrix(0.1, 4, 4)
diag(B_pref) <- 0.9
g_synth <- barabasi_albert_blocks(
  m = 2, 
  p = rep(0.25, 4), 
  B = B_pref, 
  t_max = 200, 
  type = "Hajek", 
  sample_with_replacement = FALSE
)
gt_synth <- V(g_synth)$label
analyze_network(g_synth, "Synthetic Scale-Free", gt_memb = gt_synth)


# Network 3: ENRON
data(enron, package="igraphdata")
g_enron_undirected <- as.undirected(enron, mode = "each")
E(g_enron_undirected)$weight <- 1
g_enron <- simplify(g_enron_undirected, remove.multiple = TRUE, remove.loops = TRUE, 
                    edge.attr.comb = list(weight = "sum", "ignore"))
analyze_network(g_enron, "ENRON Email", n_iters = 10)


# Network 4: User Choice (UKfaculty)
data(UKfaculty, package="igraphdata")
g_faculty_undirected <- as.undirected(UKfaculty, mode = "each")
E(g_faculty_undirected)$weight <- 1
g_faculty <- simplify(g_faculty_undirected, remove.multiple = TRUE, remove.loops = TRUE, 
                      edge.attr.comb = list(weight = "sum", "ignore"))
analyze_network(g_faculty, "UK Faculty", n_iters = 10)
