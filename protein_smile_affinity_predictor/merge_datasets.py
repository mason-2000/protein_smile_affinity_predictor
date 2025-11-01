#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

def save_output(merged_dataset, out_path='encoded_data/'):
    
    merged_dataset.to_csv(out_path + 'merged_dataset.csv')
    print(f'[INFO] dataset saved in {out_path + "merged_dataset.csv"}')

    
def merge_datasets(encoded_proteins, encoded_drugs, res_targets):
    """
    Combine protein and drug embeddings into a single dataset.

    Args:
        encoded_proteins (pd.DataFrame): Protein features.
        encoded_drugs (pd.DataFrame): Drug features.
        affinity_data (pd.DataFrame): Contains 'Protein_Index', 'Drug_Index', and 'Affinity'.

    Returns:
        pd.DataFrame: Merged dataset for model training.
    """
    
    all_features = list(res_encoded_proteins.columns) + list(res_encoded_drugs.columns)
    merged_rows = []
    
    for index, row in targets.iterrows():
        
        p_feats = res_encoded_proteins.loc[row['Protein_Index'], :].values
        d_feats = res_encoded_drugs.loc[row['Drug_Index'], :].values
        
        merged_rows.append(np.concatenate([p_feats, d_feats]))
    
    merged_dataset = pd.DataFrame(merged_rows, columns=all_features, index=targets.index)
    res_targets.drop(['Protein_Index', 'Drug_Index'], axis=1, inplace=True)
    res_targets = res_targets.squeeze()
    
    merged_dataset = pd.concat([merged_dataset, res_targets], axis=1)
    
    print('[INFO] dataset merge completed')
    
    return merged_dataset
        

def rescale_features(encoded_drugs, encoded_proteins, targets, scaler_proteins=None, scaler_drugs=None, scaler_targets=None, chem_only=False,
                     scale_targets=False):
    
    if scaler_proteins == None: scaler_proteins = StandardScaler()
    if scaler_drugs == None: scaler_drugs = StandardScaler()
    if scaler_targets == None: scaler_targets = StandardScaler()

    
    if chem_only:
        
        drugs_chem = [column for column in encoded_drugs.columns if 'emb_d' not in column]
        proteins_chem = [column for column in encoded_proteins.columns if 'emb_p' not in column]
        
        scaled_proteins_chem = scaler_proteins.fit_transform(encoded_proteins[proteins_chem])
        scaled_drugs_chem = scaler_drugs.fit_transform(encoded_drugs[drugs_chem])
        
        encoded_proteins[proteins_chem] = scaled_proteins_chem
        encoded_drugs[drugs_chem] = scaled_drugs_chem
    
    else:
        
        scaled_proteins = scaler_proteins.fit_transform(encoded_proteins)
        scaled_drugs = scaler_drugs.fit_transform(encoded_drugs)
        
        encoded_proteins = pd.DataFrame(data=scaled_proteins, columns=encoded_proteins.columns, index=encoded_proteins.index)
        encoded_drugs = pd.DataFrame(data=scaled_drugs, columns=encoded_drugs.columns, index=encoded_drugs.index)
    
    if scale_targets:
        
        scaled_affinity = scaler_targets.fit_transform(targets['Affinity'].to_frame())
        targets['Affinity'] = scaled_affinity.flatten()
        
    return encoded_proteins, encoded_drugs, targets

def load_data(base_path='encoded_data/'):
    

        encoded_drugs = pd.read_csv(base_path + 'encoded_drugs.csv', index_col=0)
        encoded_proteins = pd.read_csv(base_path + 'encoded_proteins.csv', index_col=0)
        targets = pd.read_csv(base_path + 'drug_protein_affinity.csv')

        return encoded_drugs, encoded_proteins, targets
    
if __name__ == '__main__':
    
    encoded_drugs, encoded_proteins, targets = load_data()
    res_encoded_proteins, res_encoded_drugs, res_targets = rescale_features(encoded_drugs, encoded_proteins, targets)
    merged_dataset = merge_datasets(res_encoded_proteins, res_encoded_drugs, res_targets)
    save_output(merged_dataset)
    
    
    
    
    
    