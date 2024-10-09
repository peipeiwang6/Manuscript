import os
import sys
import joblib
import pandas as pd
import numpy as np
from sklearn.metrics import f1_score
from scipy.sparse import csr_matrix

Exp_matrix = sys.argv[1]
gene_classification = sys.argv[2]
Models_folder = sys.argv[3]

def evaluate_models_on_new_data(new_test_data, new_test_labels, model_folder = Models_folder):
    """
    Load saved models, make predictions on a new test set, and calculate the F1 score.
    Parameters:
    - new_test_data: Features of the new test set (DataFrame or numpy array)
    - new_test_labels: Labels of the new test set (numpy array)
    - model_folder: Path to the folder containing the saved models
    Returns:
    - f1_scores: F1 scores of all models (DataFrame)
    """
    f1_scores = []
    if not os.path.exists(model_folder):
        raise FileNotFoundError(f"Model folder {model_folder} does not exist!")
    model_files = [f for f in os.listdir(model_folder) if f.endswith('.pkl')]
    model_files.sort()
    for model_file in model_files:
        model_path = os.path.join(model_folder, model_file)
        try:
            model = joblib.load(model_path)
            print(f"Loaded model {model_file}")
            test_pred = model.predict(new_test_data)
            f1 = f1_score(new_test_labels, test_pred, average=None)
            f1_scores.append(f1)
        except Exception as e:
            print(f"Failed to evaluate model {model_file}: {e}")
            continue
    f1_scores_df = pd.DataFrame(f1_scores)
    return f1_scores_df

gene_class_df = pd.read_csv(gene_classification)
gene_unique = gene_class_df.groupby(['Gene'])['Classification'].apply(lambda x: ','.join(x.values)).reset_index()
gene_test = gene_unique.sample(frac=0.2, random_state=1613).reset_index(drop=True)
gene_which_class_df = gene_class_df.groupby(['Gene'])['Classification'].apply(lambda x: ','.join(x.values))
gene_which_class_df = gene_which_class_df.reset_index()
gene_which_class_df.set_index('Gene', inplace=True)
Gene_classification = {word: index for index, word in enumerate(gene_class_df.Classification.value_counts().index)}
which_classification = {word: 1 for index, word in enumerate(gene_class_df.Classification.value_counts().index)}
rows = []
cols = []
data = []
for i in range(len(gene_which_class_df)):
    for class_list in gene_which_class_df.values[i]:
        for clas in class_list.split(','):
            rows.append(i)
            cols.append(Gene_classification.get(clas))
            data.append(which_classification.get(clas))

class_csr_matrix = csr_matrix((data, (rows, cols)), shape=(len(gene_which_class_df), 14))
class_csr_matrix_df = pd.DataFrame(class_csr_matrix.toarray())
class_csr_matrix_df.columns = [x for x in which_classification.keys()]
class_csr_matrix_df = pd.concat([gene_which_class_df.reset_index().drop('Classification', axis=1), class_csr_matrix_df], axis=1)

total_gene_TPM_df = pd.read_csv(Exp_matrix)
Exp_gene_TPM_df = pd.merge(class_csr_matrix_df.loc[:,'Gene'], total_gene_TPM_df, left_on = "Gene", right_on = "Gene", how = "left").fillna(0)
gene_test_matrix_df = pd.merge(gene_test.loc[:,'Gene'], Exp_gene_TPM_df, on="Gene", how="left")

gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)

new_test_data = gene_test_matrix_df.iloc[:, 1:]
f1_scores_df = evaluate_models_on_new_data(new_test_data, test_multilabel) #calculate the f1
f1_scores_df.columns = class_csr_matrix_df.columns[1:]
f1_scores_df.to_csv("F1_scores_new_test.csv", index=False)
