import os
import sys
import pandas as pd
import multiprocessing
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import f1_score
from autogluon.tabular import TabularDataset, TabularPredictor
from autogluon.common.utils.utils import setup_outputdir
from autogluon.core.utils.loaders import load_pkl
from autogluon.core.utils.savers import save_pkl
import os.path
import concurrent.futures

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

    def fit(self, train_data, tuning_data=None, hyperparameters=None, **kwargs):
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
            predictor.fit(train_data=train_data, tuning_data=tuning_data, hyperparameters=hyperparameters, **kwargs)
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

def F1score_Calculate(gene_class_df, txt, AnyData, AnyResult, raw_name, unknown_name):
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
    gene_train_matrix_df = pd.merge(gene_train.loc[:,'Gene'], Exp_gene_TPM_df, on="Gene", how="left")
    gene_train_multilabel = pd.merge(gene_train_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
    input_gene_train_multilabel = np.array(gene_train_multilabel[[x for x in list(gene_train_multilabel)[1:]]] == 1)
    gene_test_matrix_df = pd.merge(gene_test.loc[:,'Gene'], Exp_gene_TPM_df, on="Gene", how="left")
    gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
    test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)

    if 'Autogluon' in txt:
        Auto_Input_df = pd.concat([gene_train_matrix_df, gene_train_multilabel.drop(['Gene'], axis=1)], axis=1)
        label = list(gene_train_multilabel)[1:]
        Real_train_name = AnyResult + '/Autogluon' + "_Real_Score_train_" + txt
        Real_test_name = AnyResult + '/Autogluon' + "_Real_Score_test_" + txt
        unknown_pred_name = AnyResult + '/Autogluon' + "_Unknown_Predict_label_" + txt
        Global_train_name = AnyResult + '/Autogluon' + "_Real_Overall_train_" + txt
        Global_test_name = AnyResult + '/Autogluon' + "_Real_Overall_test_" + txt
        out = []
        test_f1 = []
        for run in range(10):
            score_list = []
            Real_gene_df = Auto_Input_df.sample(frac=1, random_state=run).reset_index(drop=True)
            save_path = 'Result_All/Real_model_file_%s/Models-predictEducationClass%s' % (txt.split('.')[0], run)
            problem_types = ['binary'] * 14
            multi_predictor = MultilabelPredictor(labels=label, problem_types=problem_types, eval_metrics=['f1'] * 14, path=save_path)
            multi_predictor.fit(Real_gene_df.iloc[:,1:], num_bag_folds=5, hyperparameters={'KNN': {}},auto_stack=False)
            for ll in label:
                score = pd.DataFrame(multi_predictor.get_predictor(ll).leaderboard()).score_val.to_list()[0]
                score_list.append(score)
            out.append(score_list)

            feature = []
            for lab in label:
                df = multi_predictor.get_predictor(lab).feature_importance(data=Auto_Input_df.iloc[:, 1:], features=list(gene_train_matrix_df)[1:],
                                                                           feature_stage='original', subsample_size=5000, time_limit=None, num_shuffle_sets=10,
                                                                           include_confidence_band=True, confidence_level=0.99, silent=False)
                df['Classification'] = lab
                feature.append(df)
            importance_name = f"feature_importance_Auto_{run+1}.csv"
            pd.concat(feature, axis=0).to_csv(importance_name)

            gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:, 'Gene'], class_csr_matrix_df, on="Gene", how="left")
            test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)
            X_test, y_test = gene_test_matrix_df.iloc[:, 1:], test_multilabel
            test_pred = multi_predictor.predict(X_test)
            test_scores = f1_score(y_test, test_pred, average=None)
            row_name = gene_test['Gene']
            test_pred_df = pd.DataFrame(test_pred)
            test_pred_label = f"test_label_Auto_{run+1}.csv"
            pd.concat([row_name, test_pred_df], axis=1).to_csv(test_pred_label, index=False)
            test_f1.append(test_scores)

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
        overall_F1CV = df_t2.mean(axis=1)
        df_t2['Global_F1_Score'] = overall_F1CV
        overall_train = df_t2['Global_F1_Score']
        overall_train.to_csv(Global_train_name, index=False)

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
        overall_F1test = test_df_t2.mean(axis=1)
        test_df_t2['Global_F1_Score'] = overall_F1test
        overall_test = test_df_t2['Global_F1_Score']
        overall_test.to_csv(Global_test_name, index=False)

        class_csr_matrix_df2 = class_csr_matrix_df.set_index(['Gene'])
        unknown_gene_exp = pd.read_csv(unknown_name)
        unknown_exp_df = unknown_gene_exp.iloc[:, 1:]
        unknown_pred_Auto = multi_predictor.predict(unknown_exp_df)
        unknown_pred_Auto = pd.DataFrame(unknown_pred_Auto)
        unknown_pred_Auto.columns = class_csr_matrix_df2.columns
        unknown_Auto_pred_result = unknown_pred_Auto.apply(lambda row: tuple(unknown_pred_Auto.columns[row == 1]), axis=1)
        unknown_Auto_pred_df = pd.DataFrame(unknown_Auto_pred_result)
        unknown_Auto_pred_df.columns = ['Classification']
        unknown_Auto_pred_class = pd.concat([unknown_gene_exp.loc[:,'Gene'], unknown_Auto_pred_df], axis=1)
        unknown_Auto_pred_class.to_csv(unknown_pred_name, index=False)

# calculate backgroung F1 scores
def Bg_Value(gene_class_df, txt, AnyData, AnyResult, raw_name, raw_test):
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
    Exp_gene_TPM_df = pd.merge(class_csr_matrix_df.loc[:, 'Gene'], total_gene_TPM_df, on="Gene", how="left").fillna(0)
    gene_train = pd.read_csv(raw_name)
    gene_test = pd.read_csv(raw_test)
    gene_train_matrix_df = pd.merge(gene_train.loc[:, 'Gene'], Exp_gene_TPM_df, on="Gene", how="left")
    gene_train_multilabel = pd.merge(gene_train_matrix_df.loc[:, 'Gene'], class_csr_matrix_df, on="Gene", how="left")
    input_gene_train_multilabel = np.array(gene_train_multilabel[[x for x in list(gene_train_multilabel)[1:]]] == 1)
    gene_test_matrix_df = pd.merge(gene_test.loc[:,'Gene'], Exp_gene_TPM_df, on="Gene", how="left")
    gene_test_multilabel = pd.merge(gene_test_matrix_df.loc[:,'Gene'], class_csr_matrix_df, on="Gene", how="left")
    test_multilabel = np.array(gene_test_multilabel[[x for x in list(gene_test_multilabel)[1:]]] == 1)

    if 'Autogluon' in txt:
        label = list(gene_train_multilabel)[1:]
        problem_types = ['binary'] * 14
        Bg_train_name = AnyResult + '/Autogluon' + "_Bg_Score_train_" + txt
        Bg_test_name = AnyResult + '/Autogluon' + "_Bg_Score_test_" + txt
        Global_Bg_train = AnyResult + '/Autogluon' + "_Overall_Bg_train_" + txt
        Global_Bg_test = AnyResult + '/Autogluon' + "_Overall_Bg_test_" + txt
        Bg_train = []
        Bg_test = []
        for num in range(0, 1000):
            save_path = 'Result_All/Overall_BG_train_model_file_%s/agModels-predictEducationClass%s' % (txt.split('.')[0], num)
            bg_score_list = []
            Bg_gene_df = pd.concat([gene_train_matrix_df.sample(frac=1, random_state=num).reset_index(drop=True), gene_train_multilabel.drop(['Gene'], axis=1)], axis=1)
            multi_predictor = MultilabelPredictor(labels=label, problem_types=problem_types, eval_metrics=['f1'] * 14, path=save_path)
            multi_predictor.fit(Bg_gene_df.iloc[:, 1:], num_bag_folds=5, hyperparameters={'KNN': {}},auto_stack=False)
            for ll in label:
                score = pd.DataFrame(multi_predictor.get_predictor(ll).leaderboard()).score_val.to_list()[0]
                bg_score_list.append(score)
            Bg_train.append(bg_score_list)
            X_test, y_test = gene_test_matrix_df.iloc[:, 1:], test_multilabel
            test_pred = multi_predictor.predict(X_test)
            test_scores = f1_score(y_test, test_pred, average=None)
            Bg_test.append(test_scores)
        Bg_trian_df = pd.DataFrame(Bg_train)
        Bg_trian_df.columns = list(gene_train_multilabel)[1:]
        Bg_trian_df.to_csv(Bg_train_name, index = False)
        BG_overall_F1CV = Bg_trian_df.mean(axis=1)
        Bg_trian_df['Global_F1_Score'] = BG_overall_F1CV
        BG_overall_train = Bg_trian_df['Global_F1_Score']
        BG_overall_train.to_csv(Global_Bg_train, index=False)
        Bg_test_df = pd.DataFrame(Bg_test)
        Bg_test_df.columns = list(gene_test_multilabel)[1:]
        Bg_test_df.to_csv(Bg_test_name, index=False)
        BG_overall_F1test = Bg_test_df.mean(axis=1)
        Bg_test_df['Global_F1_Score'] = BG_overall_F1test
        BG_overall_test = Bg_test_df['Global_F1_Score']
        BG_overall_test.to_csv(Global_Bg_test, index=False)

if __name__ == '__main__':
    Path = os.getcwd()
    os.chdir(Path)
    gene_classification = sys.argv[1]
    gene_class_df = pd.read_csv(gene_classification)
    gene_unique = gene_class_df.groupby(['Gene'])['Classification'].apply(lambda x: ','.join(x.values)).reset_index()
    gene_test = gene_unique.sample(frac=0.2, random_state=1613).reset_index(drop=True)
    gene_train = gene_unique[~gene_unique.Gene.isin(gene_test.Gene.to_list())].reset_index(drop=True)
    gene_train.to_csv('gene_train_raw.csv', index=False)
    gene_test.to_csv('gene_test_raw.csv', index=False)
    raw_name = 'gene_train_raw.csv'
    raw_test = 'gene_test_raw.csv'
    txt = sys.argv[2]
    unknown_name = sys.argv[3]
    AnyData = Path + "/Data"
    AnyResult = Path + "/Result_All"
    F1score_Calculate(gene_class_df, txt, AnyData, AnyResult, raw_name, unknown_name)
    Bg_Value(gene_class_df, txt, AnyData, AnyResult, raw_name, raw_test)
    Real_train = AnyResult + '/Autogluon' + "_Real_Score_train_" + txt
    Real_test = AnyResult + '/Autogluon' + "_Real_Score_test_" + txt
    Bg_train = AnyResult + '/Autogluon' + "_Bg_Score_train_" + txt
    Bg_test = AnyResult + '/Autogluon' + "_Bg_Score_test_" + txt
    Global_real_train = AnyResult + '/Autogluon' + "_Real_Overall_train_" + txt
    Global_real_test = AnyResult + '/Autogluon' + "_Real_Overall_test_" + txt
    Global_Bg_train = AnyResult + '/Autogluon' + "_Overall_Bg_train_" + txt
    Global_Bg_test = AnyResult + '/Autogluon' + "_Overall_Bg_test_" + txt
    outPut_train = 'train_' + txt.split('.')[0]
    outPut_test = 'test_' + txt.split('.')[0]
    Global_outPut_train = 'Global_train_' + txt.split('.')[0]
    Global_outPut_test = 'Global_test_' + txt.split('.')[0]
    R_density_plot_script_train = "Rscript Simulation_F1_density_plot.R" + " " + Bg_train + " " + Real_train + " " + outPut_train
    os.system(R_density_plot_script_train)
    R_density_plot_script_test = "Rscript Simulation_F1_density_plot.R" + " " + Bg_test + " " + Real_test + " " + outPut_test
    os.system(R_density_plot_script_test)
    R_density_plot_overall_train = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_train + " " + Global_real_train + " " + Global_outPut_train
    os.system(R_density_plot_overall_train)
    R_density_plot_overall_test = "Rscript Overall_F1_density_plot.R" + " " + Global_Bg_test + " " + Global_real_test + " " + Global_outPut_test
    os.system(R_density_plot_overall_test)