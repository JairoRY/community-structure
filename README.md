# Community Structure Analysis

This repository contains R scripts for detecting and evaluating community structures within complex networks. This project was developed as part of the Complex and Social Networks (CSN) course within the Master in Innovation and Research in Informatics (MIRI) at the Universitat Politècnica de Catalunya (UPC - BarcelonaTech).

The project focuses on the implementation and comparison of various community detection algorithms to identify cohesive subgroups within larger network architectures, assessing their performance through modularity and other structural metrics.

## Project Structure

The project is organized into modular scripts that handle data ingestion, algorithm execution, and comparative analysis:

* **dataLoader.R**: Handles the import of network data from formats such as edge lists or adjacency matrices and converts them into igraph objects.
* **communityDetection.R**: Contains the implementation of several community detection algorithms, including Louvain, Walktrap, and Girvan-Newman.
* **evaluationMetrics.R**: Calculates statistical measures to evaluate the quality of the detected communities, such as modularity, conductance, and normalized mutual information (NMI).
* **visualization.R**: Generates network plots where nodes are colored by their assigned community and produces dendrograms for hierarchical methods.
* **comparison_analysis.R**: Compares the results of different algorithms on the same dataset to determine which method best captures the underlying structure.

## Getting Started

### Prerequisites

The analysis requires R (version 4.0 or higher) and the following network analysis and visualization packages:

```r
install.packages(c("igraph", "ggplot2", "RColorBrewer", "networkD3"))

```

### Usage

1. **Clone the repository:**
```bash
git clone https://github.com/JairoRY/MIRI-CSN-Community-Structure.git
cd MIRI-CSN-Community-Structure

```


2. **Run the analysis:**
Scripts should be executed in sequence to process the network and generate the results:
```r
source("dataLoader.R")
source("communityDetection.R")
source("evaluationMetrics.R")
source("visualization.R")

```



## Algorithms Evaluated

The project explores different approaches to community detection:

* **Louvain Method**: A multi-level optimization of modularity.
* **Walktrap Algorithm**: Identifies communities via random walks, based on the idea that walks tend to stay within highly connected subgroups.
* **Girvan-Newman**: A divisive method that iteratively removes edges with the highest edge betweenness.
* **Label Propagation**: A fast algorithm where nodes adopt the majority label of their neighbors.

Would you like me to create a similar README for any other repository in this series?
