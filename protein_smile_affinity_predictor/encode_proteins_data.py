#!/usr/bin/env python3
import pandas as pd 
import torch
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import multiprocessing as mp

global model, alphabet, batch_converter
model = None
alphabet = None
batch_converter = None

def save_output(encoded_proteins, out_path='encoded_data/'):
    
    print('[INFO] saving output...')
    encoded_proteins.to_csv(out_path + 'encoded_proteins.csv')

def proteins_analyzer(seq, low_ph=3.0, high_ph=10):
    """
    Compute basic biochemical properties of a protein sequence.

    Args:
        seq (str): Protein amino acid sequence.
        low_ph (float): Lower bound for charge calculation.
        high_ph (float): Upper bound for charge calculation.

    Returns:
        list: Various physicochemical property values.
    """
    
    X = ProteinAnalysis(seq)
    
    mol_weight = X.molecular_weight()
    aromaticity = X.aromaticity()
    instability = X.instability_index()
    isoelectric_pt = X.isoelectric_point()
    epsilon_reduced, epsilon_oxidized = X.molar_extinction_coefficient()
    gravy = X.gravy()
    low_charge = X.charge_at_pH(low_ph)
    hight_charge = X.charge_at_pH(high_ph)
    helix_f, turn_f, sheet_f = X.secondary_structure_fraction()
    
    chem_properties = [mol_weight, aromaticity, instability, isoelectric_pt, epsilon_reduced, epsilon_oxidized, gravy,
                       low_charge, hight_charge, helix_f, turn_f, sheet_f]
    
    return chem_properties

def embed_single_sequence(seq: str, protein_name) -> np.ndarray:
    """
    Compute ESM2 embedding for a single protein sequence.

    Args:
        seq (str): Protein sequence (AAs).
        protein_name (str): Name or identifier.

    Returns:
        np.ndarray: Averaged embedding vector.
    """
    
    seq = seq.strip().upper()

    data = [(protein_name, seq)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[model.num_layers], return_contacts=False)
        token_reps = results['representations'][model.num_layers]

    emb = token_reps[0, 1:len(seq)+1].mean(0)
    
    print(f'[INFO] protein: {protein_name} encoded ')
    
    return emb.cpu().numpy().astype(np.float32)


def init_model():
    global model, alphabet, batch_converter

    model, alphabet = torch.hub.load('facebookresearch/esm:main', 'esm2_t6_8M_UR50D')
    batch_converter = alphabet.get_batch_converter()
    model.eval()


def load_data(base_path='raw_data/'):
    
        proteins = pd.read_csv(base_path + 'proteins.csv', index_col=0)
        proteins.drop('Accession_Number', axis=1, inplace=True)
        
        return proteins
        
if __name__ == '__main__':
    
    properties_names = ['molecular_weight', 'aromaticity', 'instability_index', 'isoelectric_point', 'molar_extinction_reduced',
                        'molar_extinction_oxidized', 'gravy', 'charge_at_pH3', 'charge_at_pH10', 'helix_fraction',
                        'turn_fraction', 'sheet_fraction']
    
    proteins = load_data()
    thread_n = 6
    seqs = [str(proteins.loc[index, 'Sequence'].strip().upper()) for index in proteins.index]
    p_names = [str(proteins.loc[index, 'Gene_Name'].strip()) for index in proteins.index]
    args = zip(seqs, p_names)
    
    
    with mp.Pool(processes=thread_n,  initializer=init_model) as p:
        
        embs = p.starmap(embed_single_sequence, args)
        chem_properties = p.map(proteins_analyzer, seqs)
        
    emb_len = len(embs[0])
    emb_col = [f'emb_p{str(c)}' for c in range(1, emb_len+1)]
    
    embs_df = pd.DataFrame(data=embs, columns=emb_col, index=proteins.index)
    chem_df = pd.DataFrame(data=chem_properties, columns=properties_names, index=proteins.index)
    encoded_proteins = pd.concat([chem_df, embs_df], axis=1)

    save_output(encoded_proteins)

        
        
    
    
    
    