Scripts for our manuscript: The utility of single-cell RNA sequencing data in predicting plant metabolic pathway genes

# 1. Clustering coefficient
The Clustering_coefficient folder contains the workflow code for calculating the clustering coefficient and background values using a Python script named `calculate_clustering_coefficient.py`. The script calls `WGCNA.R` to generate the input data necessary for calculating the clustering coefficient, and finally calls `Simulation_C_density_plot_230328.R` to visualize the results. **No manual execution of R scripts is required**.
* Ensure the required Python and R packages are installed.

## input data
The input data is provided in **CSV format** as a matrix, where each column represents a sample and each row represents a gene. It is essential that the gene column is labeled as **'Gene'**—any other name will not be recognized.


To execute this code, please run the following command:
```bash
python calculate_clustering_coefficient.py input_data.csv
```

# 2.Model buliding
The Model_buliding folder contains the workflow code for 

The results (F1 scores, feature importances, density plots, unknown gene predictions) will be stored in the `Result_All/` folder.

## input data
* train.csv: Training dataset used for model building.

## output files




