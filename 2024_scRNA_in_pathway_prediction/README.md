Scripts for our manuscript: The utility of single-cell RNA sequencing data in predicting plant metabolic pathway genes

# 1. Clustering coefficient
The Clustering_coefficient folder contains the codes to compute the clustering coefficient and the background values by the dataset and generate visualizations using R scripts. **No manual execution of R scripts is required**.
## 1. WGCNA analysis (R)
The `WGCNA.R` script performs WGCNA (Weighted Gene Co-expression Network Analysis) on the input data and generates an edge list.
## 2. Clustering Coefficient Calculation (Python):
The Python script `calculate_clustering_coefficient.py` uses the edge list to calculate clustering coefficients for both real and random simulation networks.
## 3. Density Plot (R):
The Python script then calls `Simulation_C_density_plot_230328.R`, which generates density plots comparing the clustering coefficient distributions of the real and background values.

# 2.Model buliding