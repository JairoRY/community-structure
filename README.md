# Community Structure Analysis

This repository contains R code for detecting and analyzing community structures within complex networks. The project was developed as part of the Complex and Social Networks (CSN) course within the Master in Innovation and Research in Informatics (MIRI) at the Universitat Politècnica de Catalunya (UPC).

The analysis focuses on identifying cohesive groups of nodes that are more densely connected to each other than to the rest of the network, using various algorithmic approaches and evaluation metrics.

## Project Structure

The project currently consists of the following core component:

* **community_functions.R**: A collection of functions designed to facilitate community detection and evaluate the quality of network partitions. This includes implementations or wrappers for calculating modularity, conductance, and other structural properties.

## Getting Started

### Prerequisites

To use these functions, you need to have R installed along with the igraph library, which is the primary tool for network analysis in this project:

```r
install.packages("igraph")

```

### Usage

1. **Clone the repository:**
```bash
git clone https://github.com/JairoRY/MIRI-CSN-Community-Structure.git
cd MIRI-CSN-Community-Structure

```


2. **Source the functions:**
You can load the provided functions into your R environment to use them on your network datasets:
```r
library(igraph)
source("community_functions.R")

```



## Algorithms Evaluated

The analysis typically involves comparing different community detection strategies to understand their strengths and weaknesses in various network topologies:

* **Louvain Method**: A heuristic method based on modularity optimization that is highly efficient for large networks.
* **Walktrap Algorithm**: Detects communities by using random walks, under the premise that short random walks tend to stay within the same community.
* **Fast-Greedy**: A bottom-up hierarchical approach that merges communities to maximize modularity at each step.
* **Edge Betweenness (Girvan-Newman)**: A divisive algorithm that removes edges with the highest betweenness score to reveal the underlying community structure.
