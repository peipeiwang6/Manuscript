### **Scripts for our manuscript: Usefulness of scRNA-seq data in predicting plant metabolic pathway genes**

# 1. Clustering coefficient
The [/Clustering_coefficient](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Clustering_coefficient) folder contains the code for calculating the clustering coefficient and background values using a Python script named `calculate_clustering_coefficient.py`. The script calls `WGCNA.R` to generate the input data necessary for calculating the clustering coefficient, and finally calls `Simulation_C_density_plot_230328.R` to visualize the results. **No manual execution of R scripts is required**.
* Ensure the required Python and R packages are installed.

## input data
* Expression matrix: the input expression data is provided in **CSV format** as a matrix, where each column represents a sample and each row represents a gene. It is essential that the gene column is labeled as **'Gene'** — any different name will not be recognized.
* Gene classification: contains two columns: one for genes and another for their classification information. The first column should match the **"Gene"** column in the expression matrix, and the second column must be named **"Classification"**. 

*Example data can be found in the folder* [/Example_data_for_clustering_coefficient_calculating](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Example_data_for_clustering_coefficient_calculating).

To execute this code, please run the following command:
```bash
python calculate_clustering_coefficient.py input_data.csv
```

# 2.Model buliding
The [/Model_building](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Model_buliding) folder contains the code for building machine learning models using K-Nearest Neighbors (KNN), eXtreme Gradient Boosting (XGBoost), and Random Forest (RF) using an expression matrix. The data is split into 80% for training and 20% for testing. The folder also includes an R script for generating visualizations of model performance. **No manual execution of R scripts is required**. The models will be saved after training, allowing for future use without retraining. 

## input files
* Expression matrix: the expression matrix is in CSV format stored in the `Data/` folder. Rows represent genes, and columns represent samples. The first column must be labeled **"Gene"** and contain the gene names. The file name must **start with** "KNN",  "XGBoost", or "RandomForest" to indicate which model will be trained.
    
    Example file names:
    - `KNN_expression_matrix.csv`
    - `XGBoost_expression_matrix.csv`
    - `RandomForest_expression_matrix.csv`

* Gene classification: contains two columns: one for genes and another for their classification information. The first column should match the **"Gene"** column in the expression matrix, and the second column must be named **"Classification"**. 
* Unknown gene expression matrix: the format of the unknown gene expression matrix is the same as that of the expression matrix.

## output files
* The results will include F1 scores from cross-validation and test sets, F1 scores from random simulations, feature importance rankings, density plots, and predictions for unknown genes.
* The trained model is saved in the `Result_All/` folder.

*Place the files and codes in their respective folders according to the structure provided in the folder* [/Example_files_for_model_building](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Example_data).

Running the KNN, XGBoost, or Random Forest models:
```bash
python KNN_XGBoost_RF_models.py expression_matrix.csv gene_classification.csv unknown_gene_expression.csv
```

# 3. Loading the saved models
The loaded model predicts labels for the new data, and the F1 scores are calculated for each model afterward. The pre-trained models are stored in the specified `model_folder` with a `.pkl` extension. To use the code, please provide the new dataset and its corresponding classifications (labels), along with the folder containing the saved models.

```bash
python load_KNN_RF_XGBoost_models.py new_data.csv new_labels.csv model_folder/
```
