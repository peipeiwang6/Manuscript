### **Scripts for our manuscript: Usefulness of scRNA-seq data in predicting plant metabolic pathway genes**

# 1. WGCNA analysis
The [/WGCNA_analysis](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/WGCNA_analysis) folder contains the code used to calculate gene co-expression within a gene expression matrix based on weighted gene co-expression network analysis (WGCNA). The analysis is implemented in R, and the main script is named `WGCNA.R`.

# 2. Clustering coefficient
The [/Clustering_coefficient](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Clustering_coefficient) folder contains the code for calculating the clustering coefficient and background values using a Python script named `calculate_clustering_coefficient.py`. The script calls `WGCNA.R` to generate the input data necessary for calculating the clustering coefficient, and finally calls `Simulation_C_density_plot_230328.R` to visualize the results. **No manual execution of R scripts is required**.
* Ensure the required Python and R packages are installed.

## input data
* Expression matrix: the input expression data is provided in **CSV format** as a matrix, where each column represents a sample and each row represents a gene. It is essential that the gene column is labeled as **"Gene"** — any different name will not be recognized.
* Gene classification: contains two columns: one for genes and another for their classification information. The first column should match the **"Gene"** column in the expression matrix, and the second column must be named **"Classification"**. 

*Example data can be found in the folder* [/Example_data_for_clustering_coefficient_calculating](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Example_data_for_clustering_coefficient_calculating).

To execute this code, please run the following command:
```bash
python calculate_clustering_coefficient.py input_data.csv
```

# 3.Model buliding
The [/Model_building](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Model_buliding) folder contains the code for building machine learning models based on different algorithms, including FASTAI, neural network (NN), K-Nearest Neighbors (KNN), eXtreme Gradient Boosting (XGBoost), and Random Forest (RF), using a gene expression matrix (`model_building_code.py`). In addition, the models corresponding to these algorithms are also built using AutoGluon based on the code in `AutoGluon_model_building.py`, where the models to be built can be specified and modified using the ‘hyperparameters’ argument. The data is split into 80% for training and 20% for testing. The folder also includes an R script for generating visualizations of model performance. **No manual execution of R scripts is required**. The models will be saved after training, allowing for future use without retraining. 

## input files
* Expression matrix: the expression matrix is in CSV format stored in the `Data/` folder. Rows represent genes, and columns represent samples. The first column must be labeled **"Gene"** and contain the gene names. The file name must **start with** "fastai", "NeuralNetwork", "KNN",  "xgboost", or "RandomForest" to indicate which model will be trained.
    
    Example file names:
    - `fastai_expression_matrix.csv`
    - `NeuralNetwork_expression_matrix.csv`    
    - `KNN_expression_matrix.csv`
    - `xgboost_expression_matrix.csv`
    - `RandomForest_expression_matrix.csv`

* Gene classification: contains two columns: one for genes and another for their classification information. The first column should match the **"Gene"** column in the expression matrix, and the second column must be named **"Classification"**. 
* Unknown gene expression matrix: the format of the unknown gene expression matrix is the same as that of the expression matrix.

## output files
* The results will include F1 scores from cross-validation and test sets, F1 scores from random simulations, feature importance rankings, density plots, and predictions for unknown genes.
* The trained model is saved in the `Result_All/` folder.

*Place the files and codes in their respective folders according to the structure provided in the folder* [/Example_files_for_model_building](https://github.com/peipeiwang6/Manuscript/tree/main/2024_scRNA_in_pathway_prediction/Example_data).

Running the FASTAI, NN, KNN, XGBoost, or Random Forest models:
```bash
python model_building_code.py expression_matrix.csv gene_classification.csv unknown_gene_expression.csv
```

# 4. Loading the saved model
This script loads a pre-trained multi-label classification saved model file (`.pkl`) produced by `model_building_code.py` and predicts the functional classification of genes in a new expression dataset. It outputs a CSV file containing Gene and comma-separated predicted classification names. To run the prediction, provide the prediction script, a new gene expression dataset (CSV) containing a *"Gene"** column and **the same feature columns** used during training, a gene classification file named `class_lables.csv` (a single column of all class types, named **"Classification"**), and a saved model file (.pkl).
    
    Supported algorithms:
    - `KNN`
    - `RandomForest`   
    - `xgboost`
    - `fastai`
    - `NeuralNetwork`   

```bash
python Load_model.py 
        --model_type <MODULE> \
        --model_path <MODEL_PATH> \
        --data_path new_expression_matrix.csv \
        --label_columns class_labels.csv
```