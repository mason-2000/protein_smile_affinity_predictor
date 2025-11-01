#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error
from joblib import dump
import time

def save_output(best_model, best_rmsd, best_w, out_path='model/protein_smile_LR.joblib'):
    
    dump(best_model, out_path)
    print(f'[INFO] model saved in {out_path}')
    
    with open('model/model_summary.txt', 'w') as log:
        
        log.write(f'rmsd: {best_rmsd:.3f}\nchemical weight: {best_w} ')
    
    return

def check_best_model(results):
    
    best_model, best_rmsd, best_w = min(results, key=lambda x: x[1])
    
    return best_model, best_rmsd, best_w

def chemical_weight_tuner(dataset, chemical_col, w):
    """
    Scale chemical (non-embedding) features by a given weight.

    Args:
        dataset (pd.DataFrame): Input dataset.
        chemical_col (list): Columns corresponding to chemical descriptors.
        weight (float): Scaling factor.

    Returns:
        pd.DataFrame: Weighted dataset.
    """

    
    w_dataset = dataset.copy()
    w_dataset[chemical_col] = w_dataset[chemical_col] * w
    
    return w_dataset

def tune_model(dataset, n_threads=6, weights=np.logspace(-1, 2, num=10)):
    """
    Tune ElasticNet hyperparameters by testing different weights
    for chemical vs embedding features.

    Args:
        dataset (pd.DataFrame): Input merged dataset.
        n_threads (int): Number of parallel jobs for training.
        weights (iterable): Range of chemical weights to test.

    Returns:
        tuple: (best_model, best_rmsd, best_weight)
    """
    chemical_col = [column for column in dataset.columns if 'emb_' not in column and column != 'Affinity']
    results = []
    
    for w in weights:
        
        w_dataset = chemical_weight_tuner(dataset, chemical_col, w)
        x_train, x_test, y_train, y_test = split_dataset(w_dataset)
        cv = KFold(n_splits=5)
        
        model = ElasticNetCV(
        l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
        cv=cv,
        max_iter=20000,
        n_jobs=n_threads,
        random_state=42)
        
        print('[INFO] starting the train phase...')
        t_start = time.monotonic()
        model.fit(x_train, y_train)
        t_end = time.monotonic()
        print(f'[INFO] training end in {(t_end-t_start):.2f}s')
        
        y_pred = model.predict(x_test)
        
        rmsd = mean_squared_error(y_test, y_pred)
        print(f'[INFO] rmsd for weight {w} is: {rmsd:.3f}')
        
        results.append((model, rmsd, w))
        
    best_model, best_rmsd, best_w = check_best_model(results)
    
    return best_model, best_rmsd, best_w
    
def split_dataset(dataset):
    
    targets = dataset['Affinity']
    train = dataset.drop('Affinity', axis=1)
    
    x_train, x_test, y_train, y_test = train_test_split(
    train, targets, test_size=0.15, random_state=42)
    
    return x_train, x_test, y_train, y_test

def load_data(input_path='encoded_data/merged_dataset.csv'):
    
    dataset = pd.read_csv(input_path, index_col=0)
    
    return dataset

if __name__ == '__main__':
    
    dataset = load_data()
    best_model, best_rmsd, best_w = tune_model(dataset)
    save_output(best_model, best_rmsd, best_w)
    
    