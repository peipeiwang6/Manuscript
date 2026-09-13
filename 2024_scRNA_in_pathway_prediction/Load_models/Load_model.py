#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import numpy as np
import joblib
import torch
from fastai.tabular.all import *
from sklearn.preprocessing import MultiLabelBinarizer

def predict_with_model(model_type, model_path, data_path, label_columns_file=None, threshold=0.2):
    print(f"Loading model from {model_path} ...")
    model = joblib.load(model_path)
    data_df = pd.read_csv(data_path)
    if 'Gene' not in data_df.columns:
        raise ValueError("Data must contain a 'Gene' column.")
    feature_cols = [col for col in data_df.columns if col != 'Gene']
    X_new = data_df[feature_cols]
    if model_type.upper() == 'KNN':
        if label_columns_file is None:
            raise ValueError("For KNN, you must provide label_columns_file.")
        label_df = pd.read_csv(label_columns_file)
        if 'Classification' not in label_df.columns:
            raise ValueError("label_columns_file must contain a 'Classification' column.")
        label_names = label_df['Classification'].tolist()
        pred_bin = model.predict(X_new)
        pred_df = pd.DataFrame(pred_bin, columns=label_names)
        pred_labels = pred_df.apply(lambda row: ','.join(row[row == 1].index), axis=1)
        result = pd.DataFrame({
            'Gene': data_df['Gene'],
            'Classification': pred_labels
        })
    elif model_type.upper() == 'FASTAI':
        if hasattr(model, 'dls') and hasattr(model.dls, 'vocab'):
            label_names = model.dls.vocab
        else:
            if label_columns_file is None:
                raise ValueError("fastai model does not contain vocab info; please provide label_columns_file.")
            label_df = pd.read_csv(label_columns_file)
            if 'Classification' not in label_df.columns:
                raise ValueError("label_columns_file must contain a 'Classification' column.")
            label_names = label_df['Classification'].tolist()
        test_dl = model.dls.test_dl(X_new)
        pred_probs, _ = model.get_preds(dl=test_dl)
        pred_bin = (pred_probs.numpy() > threshold).astype(int)
        pred_df = pd.DataFrame(pred_bin, columns=label_names)
        pred_labels = pred_df.apply(lambda row: ','.join(row[row == 1].index), axis=1)
        result = pd.DataFrame({
            'Gene': data_df['Gene'],
            'Classification': pred_labels
        })
    else:
        raise ValueError(f"Unsupported model_type: {model_type}. Supported: 'KNN', 'fastai'")
    output_path = os.path.splitext(data_path)[0] + '_predicted.csv'
    result.to_csv(output_path, index=False)
    print(f"Predictions saved to {output_path}")
    return result

def main():
    parser = argparse.ArgumentParser(description="Predict using saved KNN or fastai model")
    parser.add_argument('--model_type', type=str, required=True, choices=['KNN', 'fastai'],
                        help="Type of model (KNN or fastai)")
    parser.add_argument('--model_path', type=str, required=True,
                        help="Path to the saved model file")
    parser.add_argument('--data_path', type=str, required=True,
                        help="Path to new data CSV (must contain 'Gene' and feature columns)")
    parser.add_argument('--label_columns', type=str, default=None,
                        help="Path to CSV file with label names (one column named 'Classification'). Required for KNN, optional for fastai.")
    parser.add_argument('--threshold', type=float, default=0.2,
                        help="Threshold for binarization (fastai only)")
    args = parser.parse_args()
    predict_with_model(
        model_type=args.model_type,
        model_path=args.model_path,
        data_path=args.data_path,
        label_columns_file=args.label_columns,
        threshold=args.threshold
    )

if __name__ == '__main__':
    main()