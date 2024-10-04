import os
import sys
import pandas as pd
import multiprocessing
from scipy.sparse import csr_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import f1_score
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from autogluon.tabular import TabularDataset, TabularPredictor
from autogluon.common.utils.utils import setup_outputdir
from autogluon.core.utils.loaders import load_pkl
from autogluon.core.utils.savers import save_pkl
import os.path
import joblib

class MultilabelPredictor():
    """ Tabular Predictor for predicting multiple columns in table.
        Creates multiple TabularPredictor objects which you can also use individually.
        You can access the TabularPredictor for a particular label via: `multilabel_predictor.get_predictor(label_i)`

        Parameters
        ----------
        labels : List[str]
            The ith element of this list is the column (i.e. `label`) predicted by the ith TabularPredictor stored in this object.
        path : str, default = None
            Path to directory where models and intermediate outputs should be saved.
            If unspecified, a time-stamped folder called "AutogluonModels/ag-[TIMESTAMP]" will be created in the working directory to store all models.
            Note: To call `fit()` twice and save all results of each fit, you must specify different `path` locations or don't specify `path` at all.
            Otherwise files from first `fit()` will be overwritten by second `fit()`.
            Caution: when predicting many labels, this directory may grow large as it needs to store many TabularPredictors.
        problem_types : List[str], default = None
            The ith element is the `problem_type` for the ith TabularPredictor stored in this object.
        eval_metrics : List[str], default = None
            The ith element is the `eval_metric` for the ith TabularPredictor stored in this object.
        consider_labels_correlation : bool, default = True
            Whether the predictions of multiple labels should account for label correlations or predict each label independently of the others.
            If True, the ordering of `labels` may affect resulting accuracy as each label is predicted conditional on the previous labels appearing earlier in this list (i.e. in an auto-regressive fashion).
            Set to False if during inference you may want to individually use just the ith TabularPredictor without predicting all the other labels.
        kwargs :
            Arguments passed into the initialization of each TabularPredictor.

    """

    multi_predictor_file = 'multilabel_predictor.pkl'

    def __init__(self, labels, path=None, problem_types=None, eval_metrics=None, consider_labels_correlation=True, **kwargs):
        if len(labels) < 2:
            raise ValueError("MultilabelPredictor is only intended for predicting MULTIPLE labels (columns), use TabularPredictor for predicting one label (column).")
        if (problem_types is not None) and (len(problem_types) != len(labels)):
            raise ValueError("If provided, `problem_types` must have same length as `labels`")
        if (eval_metrics is not None) and (len(eval_metrics) != len(labels)):
            raise ValueError("If provided, `eval_metrics` must have same length as `labels`")
        self.path = setup_outputdir(path, warn_if_exist=False)
        self.labels = labels
        self.consider_labels_correlation = consider_labels_correlation
        self.predictors = {}  # key = label, value = TabularPredictor or str path to the TabularPredictor for this label
        if eval_metrics is None:
            self.eval_metrics = {}
        else:
            self.eval_metrics = {labels[i] : eval_metrics[i] for i in range(len(labels))}
        problem_type = None
        eval_metric = None
        for i in range(len(labels)):
            label = labels[i]
            path_i = self.path + "Predictor_" + label
            if problem_types is not None:
                problem_type = problem_types[i]
            if eval_metrics is not None:
                eval_metric = eval_metrics[i]
            self.predictors[label] = TabularPredictor(label=label, problem_type=problem_type, eval_metric=eval_metric, path=path_i, **kwargs)

    def fit(self, train_data, tuning_data=None, **kwargs):
        """ Fits a separate TabularPredictor to predict each of the labels.

            Parameters
            ----------
            train_data, tuning_data : str or autogluon.tabular.TabularDataset or pd.DataFrame
                See documentation for `TabularPredictor.fit()`.
            kwargs :
                Arguments passed into the `fit()` call for each TabularPredictor.
        """
        if isinstance(train_data, str):
            train_data = TabularDataset(train_data)
        if tuning_data is not None and isinstance(tuning_data, str):
            tuning_data = TabularDataset(tuning_data)
        train_data_og = train_data.copy()
        if tuning_data is not None:
            tuning_data_og = tuning_data.copy()
        else:
            tuning_data_og = None
        save_metrics = len(self.eval_metrics) == 0
        for i in range(len(self.labels)):
            label = self.labels[i]
            predictor = self.get_predictor(label)
            if not self.consider_labels_correlation:
                labels_to_drop = [l for l in self.labels if l != label]
            else:
                labels_to_drop = [self.labels[j] for j in range(i+1, len(self.labels))]
            train_data = train_data_og.drop(labels_to_drop, axis=1)
            if tuning_data is not None:
                tuning_data = tuning_data_og.drop(labels_to_drop, axis=1)
            print(f"Fitting TabularPredictor for label: {label} ...")
            predictor.fit(train_data=train_data, tuning_data=tuning_data, **kwargs)
            self.predictors[label] = predictor.path
            if save_metrics:
                self.eval_metrics[label] = predictor.eval_metric
        self.save()

    def predict(self, data, **kwargs):
        """ Returns DataFrame with label columns containing predictions for each label.

            Parameters
            ----------
            data : str or autogluon.tabular.TabularDataset or pd.DataFrame
                Data to make predictions for. If label columns are present in this data, they will be ignored. See documentation for `TabularPredictor.predict()`.
            kwargs :
                Arguments passed into the predict() call for each TabularPredictor.
        """
        return self._predict(data, as_proba=False, **kwargs)

    def predict_proba(self, data, **kwargs):
        """ Returns dict where each key is a label and the corresponding value is the `predict_proba()` output for just that label.

            Parameters
            ----------
            data : str or autogluon.tabular.TabularDataset or pd.DataFrame
                Data to make predictions for. See documentation for `TabularPredictor.predict()` and `TabularPredictor.predict_proba()`.
            kwargs :
                Arguments passed into the `predict_proba()` call for each TabularPredictor (also passed into a `predict()` call).
        """
        return self._predict(data, as_proba=True, **kwargs)

    def evaluate(self, data, **kwargs):
        """ Returns dict where each key is a label and the corresponding value is the `evaluate()` output for just that label.

            Parameters
            ----------
            data : str or autogluon.tabular.TabularDataset or pd.DataFrame
                Data to evalate predictions of all labels for, must contain all labels as columns. See documentation for `TabularPredictor.evaluate()`.
            kwargs :
                Arguments passed into the `evaluate()` call for each TabularPredictor (also passed into the `predict()` call).
        """
        data = self._get_data(data)
        eval_dict = {}
        for label in self.labels:
            print(f"Evaluating TabularPredictor for label: {label} ...")
            predictor = self.get_predictor(label)
            eval_dict[label] = predictor.evaluate(data, **kwargs)
            if self.consider_labels_correlation:
                data[label] = predictor.predict(data, **kwargs)
        return eval_dict

    def save(self):
        """ Save MultilabelPredictor to disk. """
        for label in self.labels:
            if not isinstance(self.predictors[label], str):
                self.predictors[label] = self.predictors[label].path
        save_pkl.save(path=self.path+self.multi_predictor_file, object=self)
        print(f"MultilabelPredictor saved to disk. Load with: MultilabelPredictor.load('{self.path}')")

    @classmethod
    def load(cls, path):
        """ Load MultilabelPredictor from disk `path` previously specified when creating this MultilabelPredictor. """
        path = os.path.expanduser(path)
        if path[-1] != os.path.sep:
            path = path + os.path.sep
        return load_pkl.load(path=path+cls.multi_predictor_file)

    def get_predictor(self, label):
        """ Returns TabularPredictor which is used to predict this label. """
        predictor = self.predictors[label]
        if isinstance(predictor, str):
            return TabularPredictor.load(path=predictor)
        return predictor

    def _get_data(self, data):
        if isinstance(data, str):
            return TabularDataset(data)
        return data.copy()

    def _predict(self, data, as_proba=False, **kwargs):
        data = self._get_data(data)
        if as_proba:
            predproba_dict = {}
        for label in self.labels:
            print(f"Predicting with TabularPredictor for label: {label} ...")
            predictor = self.get_predictor(label)
            if as_proba:
                predproba_dict[label] = predictor.predict_proba(data, as_multiclass=True, **kwargs)
            data[label] = predictor.predict(data, **kwargs)
        if not as_proba:
            return data[self.labels]
        else:
            return predproba_dict

class  model_selection():
    def __init__(self, matrix_df, txt):
        self.matrix_df = matrix_df
        self.txt = txt
        self.feature_importance = None
    def save_feature_importance_to_csv(self, file_name):
        if self.feature_importance is not None:
            self.feature_importance.to_csv(file_name, index_label='Feature', header=['Importance'])
        else:
            print("Feature importance is not available.")    
    def KNN_model(self, X_train, y_train):
        param_grid = [{'weights': ["uniform", "distance"], 'n_neighbors': [3, 4, 5, 6, 7]}]
        knn_clf = KNeighborsClassifier()
        grid_search = GridSearchCV(knn_clf, param_grid, cv=5)
        grid_search.fit(X_train, y_train)
        return grid_search.best_estimator_
    def RF_model(self, X_train, y_train):
        param_grid = [{'max_depth':[3, 5, 10], 'max_features': [0.1, 0.25, 0.5, 0.75, 'sqrt', 'log2', None], 'n_estimators': [10, 50, 100, 500]}]
        forest_clf = RandomForestClassifier(random_state=42)
        grid_search = GridSearchCV(forest_clf, param_grid, cv=5)
        grid_search.fit(X_train, y_train)
        return grid_search.best_estimator_
    def XGboost(self, X_train, y_train):
        param_grid = [{'max_depth': [3, 5, 10], 'learning_rate': [0.01, 0.1, 0.5, 1], 'n_estimators': [10, 50, 100],
                       'gamma': [0.1, 5]}]
        bst = XGBClassifier(tree_method='hist', objective='binary:logistic')
        grid_search = GridSearchCV(bst, param_grid, cv=5)
        grid_search.fit(X_train, y_train)
        return grid_search.best_estimator_

def Real_Score(module, gene_train_matrix_df, class_csr_matrix_df, best_model, txt, gene_test_matrix_df, unknown_gene_exp):
    out = []
    test_f1 = []
    global_train = []
    global_test = []
    Real_train_name = module + "_Real_Score_train_" + txt
    Real_test_name = module + "_Real_Score_test_" + txt
    Global_train_name = module + "_Overall_Score_train_" + txt
    Global_test_name = module + "_Overall_Score_test_" + txt

    gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
    test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)

    save_folder = 'save_models'
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for i in range(10):
        gene_train_df = gene_train_matrix_df.sample(frac=1, random_state=i).reset_index(drop=True)
        gene_train_multilabel = pd.merge(gene_train_df.loc[:, 'Gene'], class_csr_matrix_df, on="Gene", how="left")
        input_gene_train_multilabel = np.array(gene_train_multilabel[[x for x in list(gene_train_multilabel)[1:]]] == 1)
        X_train, y_train = gene_train_df.iloc[:, 1:], input_gene_train_multilabel

        model_on = model_selection(matrix_df = gene_train_df, txt = txt)
        if module == 'KNN':
            best_model = model_on.KNN_model(X_train, y_train)
        elif module == 'RandomForest':
            best_model = model_on.RF_model(X_train, y_train)
        elif module == 'xgboost':
            best_model = model_on.XGboost(X_train, y_train)

        y_train_pred = cross_val_predict(best_model, X_train, y_train, cv=5)
        scores = f1_score(y_train, y_train_pred, average=None)
        out.append(scores)
        overall_scores = f1_score(y_train, y_train_pred, average="macro")
        global_train.append(overall_scores)

        X_test, y_test = gene_test_matrix_df.iloc[:, 1:], test_multilabel
        test_pred = best_model.predict(X_test)
        test_scores = f1_score(y_test, test_pred, average=None)
        overall_test_scores = f1_score(y_test, test_pred, average="macro")
        test_pred_df = pd.DataFrame(test_pred)
        test_pred_df.columns = list(gene_train_multilabel)[1:]
        row_name = gene_test['Gene']
        test_pred_label = f"test_label_{i+1}.csv"
        pd.concat([row_name, test_pred_df], axis=1).to_csv(test_pred_label, index=False)
        test_f1.append(test_scores)
        global_test.append(overall_test_scores)

        # save models
        model_filename =os.path.join(save_folder, f'model_{i+1}.pkl')
        joblib.dump(best_model, model_filename)

        if module == 'RandomForest' or module == 'xgboost':
            feature_importances = best_model.feature_importances_    
            feature_names = X_train.columns if hasattr(X_train, 'columns') else range(len(feature_importances))   
            importance_df = pd.DataFrame({  
                'Feature': feature_names,  
                'Importance': feature_importances  
            })  
            importance_name = "feature_importance_" + module +f"_{i+1}.csv"
            importance_df.to_csv(importance_name, index=False)

    df = pd.DataFrame(out)
    label = pd.DataFrame(list(gene_train_multilabel)[1:])
    label.columns = ['Classification']
    mean = pd.DataFrame(df.mean().tolist())
    mean.columns = ['Mean']
    std = pd.DataFrame(df.std().tolist())
    std.columns = ['SD']
    df_t = df.T
    df_t2 = df_t.reset_index(drop = True)
    df_t2.columns = [f'Rep{i+1}' for i in range(10)]
    pd.concat([label, df_t2, mean, std], axis=1).to_csv(Real_train_name, index=False)
    overall_train_df = pd.DataFrame(global_train, columns=["Global_F1_Score"])
    overall_train_df.to_csv(Global_train_name, index = False)

    test_df = pd.DataFrame(test_f1)
    test_label = pd.DataFrame(list(gene_test_multilabel)[1:])
    test_label.columns = ['Classification']
    test_mean = pd.DataFrame(test_df.mean().tolist())
    test_mean.columns = ['Mean']
    test_std = pd.DataFrame(test_df.std().tolist())
    test_std.columns = ['SD']
    test_df_t = test_df.T
    test_df_t2 = test_df_t.reset_index(drop=True)
    test_df_t2.columns = [f'Rep{i+1}' for i in range(10)]
    pd.concat([test_label, test_df_t2, test_mean, test_std], axis=1).to_csv(Real_test_name, index=False)
    overall_test_df = pd.DataFrame(global_test, columns=["Global_F1_Score"])
    overall_test_df.to_csv(Global_test_name, index = False)

    unknown_exp_df = unknown_gene_exp.iloc[:, 1:]
    unknown_pred = best_model.predict(unknown_exp_df)
    unknown_pred_1 = pd.DataFrame(unknown_pred)
    unknown_row_name = unknown_gene_exp['Gene']
    unknown_pred_label = module + "_Unknown_Predict_label_" + txt
    unknown_pred_1.columns = list(gene_train_multilabel)[1:]
    pd.concat([unknown_row_name, unknown_pred_1], axis=1).to_csv(unknown_pred_label, index=False)
    class_csr_matrix_df2 = class_csr_matrix_df.set_index(['Gene'])
    unknown_pred_1.columns = class_csr_matrix_df2.columns
    unknown_pred_result = unknown_pred_1.apply(lambda row: tuple(unknown_pred_1.columns[row == 1]), axis=1)
    unknown_pred_df = pd.DataFrame(unknown_pred_result)
    unknown_pred_df.columns = ['Classification']
    unknown_pred_class = pd.concat([unknown_gene_exp.loc[:,'Gene'], unknown_pred_df], axis=1)
    unknown_pred_name = module + "_Unknown_Predict_Result_" + txt
    unknown_pred_class.to_csv(unknown_pred_name, index = False) 

def Bg_value(module, gene_train_matrix_df, gene_test_matrix_df, class_csr_matrix_df, best_model, txt):
    train_Bg = []
    test_Bg = []
    global_train_bg = []
    global_test_bg = []
    trian_bg_name = module + "_Bg_train_" + txt
    test_bg_name = module + "_Bg_test_" + txt
    Global_Bg_train = module + "_Overall_Bg_train_" + txt
    Global_Bg_test = module + "_Overall_Bg_test_" + txt
    gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
    test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)

    save_folder = 'save_Bg_models'
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for i in range(1000):
        gene_train_df = gene_train_matrix_df
        gene_train_multilabel = pd.merge(gene_train_df.loc[:, 'Gene'], class_csr_matrix_df, on="Gene", how="left")
        input_gene_train_multilabel = np.array(gene_train_multilabel[[x for x in list(gene_train_multilabel)[1:]]] == 1)
        X_train, y_train = gene_train_df.iloc[:, 1:].sample(frac=1, random_state=i).reset_index(drop=True), input_gene_train_multilabel

        model_on = model_selection(matrix_df = gene_train_df, txt = txt)
        if module == 'KNN':
            best_model = model_on.KNN_model(X_train, y_train)
        elif module == 'RandomForest':
            best_model = model_on.RF_model(X_train, y_train)
        elif module == 'xgboost':
            best_model = model_on.XGboost(X_train, y_train)

        y_train_pred = cross_val_predict(best_model, X_train, y_train, cv=5)
        scores = f1_score(y_train, y_train_pred, average=None)
        train_Bg.append(scores)
        global_bg_scores = f1_score(y_train, y_train_pred, average="macro")
        global_train_bg.append(global_bg_scores)

        X_test, y_test = gene_test_matrix_df.iloc[:, 1:], test_multilabel
        test_pred = best_model.predict(X_test)
        test_scores = f1_score(y_test, test_pred, average=None)
        test_Bg.append(test_scores)
        global_test_scores = f1_score(y_test, test_pred, average="macro")
        global_test_bg.append(global_test_scores)
        
        model_filename =os.path.join(save_folder, f'Bg_model_{i+1}.pkl') # save random simulation models
        joblib.dump(best_model, model_filename)

    train_Bg_df = pd.DataFrame(train_Bg)
    train_Bg_df.columns = list(gene_train_multilabel)[1:]
    train_Bg_df.to_csv(trian_bg_name, index=False)
    overall_bg_train_df = pd.DataFrame(global_train_bg, columns=["Global_F1_Score"])
    overall_bg_train_df.to_csv(Global_Bg_train, index = False)

    test_Bg_df = pd.DataFrame(test_Bg)
    test_Bg_df.columns = list(gene_test_multilabel)[1:]
    test_Bg_df.to_csv(test_bg_name, index=False)
    overall_bg_test_df = pd.DataFrame(global_test_bg, columns=["Global_F1_Score"])
    overall_bg_test_df.to_csv(Global_Bg_test, index = False)

def F1score_Calculate(txt, AnyData, raw_name, raw_test, unknown_name):
    gene_class_df = pd.read_csv("all_gene_classification_new.csv")
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

    File_txt = AnyData + '/' + txt
    total_gene_TPM_df = pd.read_csv(File_txt)
    Exp_gene_TPM_df = pd.merge(class_csr_matrix_df.loc[:,'Gene'], total_gene_TPM_df, left_on = "Gene", right_on = "Gene", how = "left").fillna(0)
    gene_train = pd.read_csv(raw_name)
    gene_test = pd.read_csv(raw_test)
    unknown_gene_exp = pd.read_csv(unknown_name)
    gene_train_matrix_df = pd.merge(gene_train.loc[:, 'Gene'], Exp_gene_TPM_df, left_on="Gene", right_on="Gene", how="left")
    gene_train_multilabel = pd.merge(gene_train_matrix_df.loc[:, 'Gene'], class_csr_matrix_df, left_on="Gene", right_on="Gene", how="left")
    input_gene_train_multilabel = np.array(gene_train_multilabel[[x for x in list(gene_train_multilabel)[1:]]] == 1)
    gene_test_matrix_df = pd.merge(gene_test.loc[:,'Gene'], Exp_gene_TPM_df, on="Gene", how="left")

    Path = os.getcwd()
    Out_path = Path + "/Result_All"
    if not os.path.exists(Out_path):
        os.mkdir(Out_path)
    os.chdir(Out_path)
    X_train, y_train = gene_train_matrix_df.iloc[:, 1:], input_gene_train_multilabel
    model_on = model_selection(matrix_df = gene_train_matrix_df, txt = txt)

    if 'KNN' in txt:
        best_model = model_on.KNN_model(X_train, y_train)
        Real_Score('KNN', gene_train_matrix_df, class_csr_matrix_df, best_model, txt, gene_test_matrix_df, unknown_gene_exp)
        Bg_value('KNN', gene_train_matrix_df, gene_test_matrix_df, class_csr_matrix_df, best_model, txt)
        Real_name = 'KNN' + "_Real_Score_train_" + txt
        Test_name = 'KNN' + "_Real_Score_test_" + txt
        Global_real_train = 'KNN' + "_Overall_Score_train_" + txt
        Global_real_test = 'KNN' + "_Overall_Score_test_" + txt
        Bg_name_train = 'KNN' + "_Bg_train_" + txt
        Bg_name_test = 'KNN' + "_Bg_test_" + txt
        Global_Bg_train = 'KNN' + '_Overall_Bg_train_' + txt
        Global_Bg_test = 'KNN' + '_Overall_Bg_test_' + txt
        outPut_train_name = 'KNN_train_' + txt.split('.')[0]
        outPut_test_name = 'KNN_test_' + txt.split('.')[0]
        Global_outPut_train = 'KNN_Global_train_' + txt.split('.')[0]
        Global_outPut_test = 'KNN_Global_test_' + txt.split('.')[0]
        R_density_plot_script_train = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_train + " " + Real_name + " " + outPut_train_name
        os.system(R_density_plot_script_train)
        R_density_plot_script_test = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_test + " " + Test_name + " " + outPut_test_name
        os.system(R_density_plot_script_test)
        R_density_plot_overall_train = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_train + " " + Global_real_train + " " + Global_outPut_train
        os.system(R_density_plot_overall_train)
        R_density_plot_overall_test = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_test + " " + Global_real_test + " " + Global_outPut_test
        os.system(R_density_plot_overall_test)

    if 'RandomForest' in txt:
        best_model = model_on.RF_model(X_train, y_train)
        Real_Score('RandomForest', gene_train_matrix_df, class_csr_matrix_df, best_model, txt, gene_test_matrix_df, unknown_gene_exp)
        Bg_value('RandomForest', gene_train_matrix_df, gene_test_matrix_df, class_csr_matrix_df, best_model, txt)
        Real_name = 'RandomForest' + "_Real_Score_train_" + txt
        Test_name = 'RandomForest' + "_Real_Score_test_" + txt
        Global_real_train = 'RandomForest' + "_Overall_Score_train_" + txt
        Global_real_test = 'RandomForest' + "_Overall_Score_test_" + txt
        Bg_name_train = 'RandomForest' + "_Bg_train_" + txt
        Bg_name_test = 'RandomForest' + "_Bg_test_" + txt
        Global_Bg_train = 'RandomForest' + '_Overall_Bg_train_' + txt
        Global_Bg_test = 'RandomForest' + '_Overall_Bg_test_' + txt
        outPut_train_name = 'RandomForest_train_' + txt.split('.')[0]
        outPut_test_name = 'RandomForest_test_' + txt.split('.')[0]
        Global_outPut_train = 'RandomForest_Global_train_' + txt.split('.')[0]
        Global_outPut_test = 'RandomForest_Global_test_' + txt.split('.')[0]
        R_density_plot_script_train = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_train + " " + Real_name + " " + outPut_train_name
        os.system(R_density_plot_script_train)
        R_density_plot_script_test = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_test + " " + Test_name + " " + outPut_test_name
        os.system(R_density_plot_script_test)
        R_density_plot_overall_train = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_train + " " + Global_real_train + " " + Global_outPut_train
        os.system(R_density_plot_overall_train)
        R_density_plot_overall_test = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_test + " " + Global_real_test + " " + Global_outPut_test
        os.system(R_density_plot_overall_test)

    if 'xgboost' in txt:
        best_model = model_on.XGboost(X_train, y_train)
        Real_Score('xgboost', gene_train_matrix_df, class_csr_matrix_df, best_model, txt, gene_test_matrix_df, unknown_gene_exp)
        Bg_value('xgboost', gene_train_matrix_df, gene_test_matrix_df, class_csr_matrix_df, best_model, txt)
        Real_name = 'xgboost' + "_Real_Score_train_" + txt
        Test_name = 'xgboost' + "_Real_Score_test_" + txt
        Global_real_train = 'xgboost' + "_Overall_Score_train_" + txt
        Global_real_test = 'xgboost' + "_Overall_Score_test_" + txt
        Bg_name_train = 'xgboost' + "_Bg_train_" + txt
        Bg_name_test = 'xgboost' + "_Bg_test_" + txt
        Global_Bg_train = 'xgboost' + '_Overall_Bg_train_' + txt
        Global_Bg_test = 'xgboost' + '_Overall_Bg_test_' + txt
        outPut_train_name = 'xgboost_train_' + txt.split('.')[0]
        outPut_test_name = 'xgboost_test_' + txt.split('.')[0]
        Global_outPut_train = 'xgboost_Global_train_' + txt.split('.')[0]
        Global_outPut_test = 'xgboost_Global_test_' + txt.split('.')[0]
        R_density_plot_script_train = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_train + " " + Real_name + " " + outPut_train_name
        os.system(R_density_plot_script_train)
        R_density_plot_script_test = "Rscript Simulation_F1_density_plot.R" + " " + Bg_name_test + " " + Test_name + " " + outPut_test_name
        os.system(R_density_plot_script_test)
        R_density_plot_overall_train = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_train + " " + Global_real_train + " " + Global_outPut_train
        os.system(R_density_plot_overall_train)
        R_density_plot_overall_test = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_test + " " + Global_real_test + " " + Global_outPut_test
        os.system(R_density_plot_overall_test)

    if 'Autogluon' in txt:
        pass

if __name__ == '__main__':
    Path = os.getcwd()
    os.chdir(Path)
    AnyData = Path + "/Data"
    AnyResult = Path + "/Result_All"
    gene_classification = sys.argv[1]
    gene_class_df = pd.read_csv(gene_classification)
    gene_unique = gene_class_df.groupby(['Gene'])['Classification'].apply(lambda x: ','.join(x.values)).reset_index()
    gene_test = gene_unique.sample(frac=0.2, random_state=1613).reset_index(drop=True)
    gene_train = gene_unique[~gene_unique.Gene.isin(gene_test.Gene.to_list())].reset_index(drop=True)
    gene_train.to_csv('gene_train_raw.csv', index=False)
    gene_test.to_csv('gene_test_raw.csv', index=False)
    raw_name = 'gene_train_raw.csv'
    raw_test = 'gene_test_raw.csv'
    unknown_name = sys.argv[2]
    mission = multiprocessing.Pool(len(os.listdir(AnyData)))
    for txt in os.listdir(AnyData):
        mission.apply_async(F1score_Calculate, args=(txt, AnyData, raw_name, raw_test, unknown_name))
    mission.close()
    mission.join()