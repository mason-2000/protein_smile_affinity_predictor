#!/usr/bin/env python3
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import Chem
from transformers import AutoTokenizer, AutoModel
import torch
import numpy as np
import pandas as pd
import multiprocessing as mp

global model, tokenizer
model = None
tokenizer = None

def init_model():
    """
    Initialize the ChemBERTa tokenizer and model.
    This is called once per process in the multiprocessing pool.
    """
    global model, tokenizer
    tokenizer = AutoTokenizer.from_pretrained('DeepChem/ChemBERTa-77M-MLM')
    model = AutoModel.from_pretrained('DeepChem/ChemBERTa-77M-MLM')
    model.eval()


def save_output(encoded_drugs, out_path='encoded_data/'):
    
    print('[INFO] saving output...')
    encoded_drugs.to_csv(out_path + 'encoded_drugs.csv')
    
def drugs_descriptor(smi, desc_names=['MolWt', 'ExactMolWt', 'MolLogP', 'TPSA', 'HeavyAtomMolWt', 'NumHDonors', 'NumHAcceptors', 
             'NumRotatableBonds', 'RingCount', 'NumAromaticRings', 'NumAliphaticRings', 'HeavyAtomCount', 'FractionCSP3', 
             'MolMR', 'LabuteASA', 'MolSurf']):
    """
    Compute a set of standard molecular descriptors from a SMILES string.

    Args:
        smi (str): SMILES representation of a molecule.
        desc_names (list): List of RDKit descriptor names.

    Returns:
        tuple: Computed descriptor values.
    """
    
    mol = Chem.MolFromSmiles(smi)
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
    descriptor_values = calc.CalcDescriptors(mol)
    return descriptor_values
        
def smiles_to_chemberta_embedding(smiles: str) -> np.ndarray:
    """
    Encode a molecule using the ChemBERTa model.

    Args:
        smiles (str): SMILES string to encode.

    Returns:
        np.ndarray: 1D embedding vector (float32).
    """
    
    inputs = tokenizer(smiles, return_tensors='pt', padding=True, truncation=True, max_length=128)
    
    with torch.no_grad():
        outputs = model(**inputs)
        
    embeddings = outputs.last_hidden_state.mean(dim=1)
    print(f'[INFO] SMILES: {smiles[:10]}... encoded')
    
    return embeddings[0].cpu().numpy().astype(np.float32)

def load_data(base_path='raw_data/'):
    
        drugs = pd.read_csv(base_path + 'drugs.csv', index_col=0)
        drugs.drop(['CID', 'Canonical_SMILES'], axis=1, inplace=True)

        return drugs

if __name__ == '__main__':
    
    
    desc_names = ['MolWt', 'ExactMolWt', 'MolLogP', 'TPSA', 'HeavyAtomMolWt', 'NumHDonors', 'NumHAcceptors', 
                 'NumRotatableBonds', 'RingCount', 'NumAromaticRings', 'NumAliphaticRings', 'HeavyAtomCount', 'FractionCSP3', 
                 'MolMR', 'LabuteASA', 'MolSurf']
    
    drugs = load_data()
    
    threads_n = 6
    smiles = [str(smi).strip() for smi in drugs['Isomeric_SMILES']]

    with mp.Pool(processes=threads_n, initializer=init_model) as p:
        
        descs = p.map(drugs_descriptor, smiles)    
        embs = p.map(smiles_to_chemberta_embedding, smiles)
    
    emb_len = len(embs[0])
    emb_col = [f'emb_d{str(c)}' for c in range(1, emb_len+1)]
    
    embedded_df = pd.DataFrame(data=embs, columns=emb_col, index=drugs.index)
    descriptor_df = pd.DataFrame(data=descs, columns=desc_names, index=drugs.index)
    
    encoded_drugs = pd.concat([descriptor_df, embedded_df], axis=1)
    save_output(encoded_drugs)
        
        
        
        
        
        
        
        
        
        
        