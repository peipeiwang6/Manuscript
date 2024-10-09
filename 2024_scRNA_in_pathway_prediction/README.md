**Scripts for our manuscript: The utility of single-cell RNA sequencing data in predicting plant metabolic pathway genes**

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
The Model_buliding folder contains the workflow code for buliding machine learning models using K-Nearest Neighbors (KNN), Random Forest (RF), eXtreme Gradient Boosting (XGBoost), and AutoGluon-Tablular (AutoGluon) on an expression matrix. The data is split into 80% for training and 20% for testing. It also includes an R script for generating visualizations of model performance. The models will be saved after training, allowing for future use without retraining. 

## input files
* Expression matrix: The expression matrix is in CSV format, where rows represent genes and columns represent samples. The first column must be labeled **"Gene"** and contain the gene names. The file name must **start with** "KNN", "RandomForest", "XGBoost", or "Autogluon" to indicate which model will be trained. The expression matrix will be split into 80% for training and 20% for testing.
Example file names:
- `KNN_expression_matrix.csv`
- `RandomForest_expression_matrix.csv`
- `XGBoost_expression_matrix.csv`
+ `Autogluon_expression_matrix.csv`
* Gene classification: Contains two columns: one for genes and another for their classification information. The first column should match the "Gene" column in the expression matrix.
* Unknown gene expression matrix: The format of the unknown gene expression matrix is the same as that of the expression matrix.

## output files
* The results (F1 scores on the cross-validation/test, F1 scores from random simulations, feature importances, density plots, unknown gene predictions) will be stored in the `Result_All/` folder.
* The trained model is saved in the `Result_All/` folder.

Running the KNN, Random Forest, or XGBoost models:
```bash
python KNN_RF_XGBoost_models.py expression_matrix.csv gene_classification.csv unknown_gene_expression.csv
```

Running the AutoGluon models:
```bash
python AutoGluon_models.py expression_matrix.csv gene_classification.csv unknown_gene_expression.csv
```

# 3. Loading Saved Models

