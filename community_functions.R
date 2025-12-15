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

analyze_network <- function(g, net_name, gt_memb = NULL, n_iters = 1, seed = 123) {
  cat("Network:", net_name, "\n")
  
  pick_row <- function(df, pattern) {
    idx <- grep(pattern, rownames(df), ignore.case = TRUE)
    if (length(idx) == 0) stop("Could not find metric row matching: ", pattern)
    if (length(idx) > 1) {
      warning("Multiple rows matched '", pattern, "'. Using the first: ", rownames(df)[idx[1]])
    }
    idx[1]
  }
  
  base_scores <- evaluate_significance(g, alg_list = alg_list, gt_clustering = gt_memb)
  accum_scores <- as.matrix(base_scores)
  
  if (n_iters > 1) {
    for (k in 2:n_iters) {
      next_scores <- evaluate_significance(g, alg_list = alg_list, gt_clustering = gt_memb)
      accum_scores <- accum_scores + as.matrix(next_scores)
    }
  }
  
  avg_scores_matrix <- accum_scores / n_iters
  scores <- as.data.frame(avg_scores_matrix)
  
  row_id   <- pick_row(scores, "internal density")
  row_exp  <- pick_row(scores, "expansion")
  row_cond <- pick_row(scores, "conductance")
  row_mod  <- pick_row(scores, "modularity")
  
  rows_to_use <- c(row_id, row_exp, row_cond, row_mod)
  print(scores[rows_to_use, , drop = FALSE])
  
  memberships_by_iter <- vector("list", n_iters)
  
  for (k in seq_len(n_iters)) {
    set.seed(seed + k - 1)
    memb_k <- list()
    for (alg_name in names(alg_list)) {
      comm <- suppressWarnings(alg_list[[alg_name]](g))
      memb_k[[alg_name]] <- as.numeric(membership(comm))
    }
    memberships_by_iter[[k]] <- memb_k
  }
  
  if (!is.null(gt_memb)) {
    ref_memb <- as.numeric(gt_memb)
    ref_name <- "GT"
  } else {
    algs_in_scores <- intersect(names(alg_list), colnames(scores))
    subset_scores <- scores[, algs_in_scores, drop = FALSE]
    
    vals_id   <- as.numeric(subset_scores[row_id, ])
    vals_exp  <- as.numeric(subset_scores[row_exp, ])
    vals_cond <- as.numeric(subset_scores[row_cond, ])
    vals_mod  <- as.numeric(subset_scores[row_mod, ])
    
    rank_id   <- rank(-vals_id,   ties.method = "min")
    rank_mod  <- rank(-vals_mod,  ties.method = "min")
    rank_exp  <- rank(vals_exp,   ties.method = "min")
    rank_cond <- rank(vals_cond,  ties.method = "min")
    
    rank_matrix <- rbind(rank_id, rank_exp, rank_cond, rank_mod)
    mean_ranks <- colMeans(rank_matrix)
    names(mean_ranks) <- algs_in_scores
    
    best_alg_name <- names(which.min(mean_ranks))
    cat("   Average Ranks (Lower is better):\n")
    print(mean_ranks)
    cat("-> Selected Reference:", best_alg_name, "\n")
    
    mods <- numeric(n_iters)
    for (k in seq_len(n_iters)) {
      mods[k] <- modularity(g, memberships_by_iter[[k]][[best_alg_name]])
    }
    best_k <- which.max(mods)
    ref_memb <- memberships_by_iter[[best_k]][[best_alg_name]]
    ref_name <- best_alg_name
    cat("-> Using iteration", best_k, "as reference (max modularity =", sprintf("%.4f", mods[best_k]), ")\n")
  }
  
  compare_k <- if (!is.null(gt_memb)) 1 else best_k
  
  for (alg_name in names(alg_list)) {
    alg_memb <- memberships_by_iter[[compare_k]][[alg_name]]
    
    JS <- jaccard_sim(ref_memb, alg_memb)
    MC <- match_clusters(JS, name1 = ref_name, name2 = alg_name)
    global_sim <- Wmean(ref_memb, MC)
    
    cat(paste0("\n", ref_name, " vs ", alg_name, "\n"))
    cat(sprintf("Global Jaccard: %.4f\n", global_sim))
    cat("Local Matches:\n")
    print(MC)
  }
  
  invisible(list(
    scores = scores,
    reference = list(name = ref_name, membership = ref_memb),
    memberships_by_iter = memberships_by_iter
  ))
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
  m = 4, 
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

