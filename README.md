protein_smile_affinity_predictor

The project investigates how weighting chemical properties relative to molecular and protein embeddings affects protein–ligand affinity prediction. An ElasticNet regression with cross-validation is used to evaluate the optimal weighting factor and model performance in terms of RMSD.

Drug–Protein Interaction Encoding and Modeling Pipeline

This repository contains a complete data processing and modeling pipeline for encoding, merging, and modeling drug–protein interaction data.
The workflow transforms raw biochemical data (SMILES and amino acid sequences) into high-dimensional representations using ChemBERTa and ESM2, enriches them with traditional chemical and physicochemical descriptors, merges both representations, and finally trains an ElasticNet regression model to predict affinity scores.

Overview

The pipeline performs the following main steps:

1. Encode Drug Molecules (encode_drugs_data.py)

Computes traditional molecular descriptors using RDKit.

Generates molecular embeddings with ChemBERTa (DeepChem/ChemBERTa-77M-MLM).

Saves the combined features to encoded_data/encoded_drugs.csv.

2. Encode Protein Sequences (encode_proteins_data.py)

Calculates physicochemical and structural descriptors with BioPython.

Generates ESM2 embeddings from amino acid sequences.

Saves results to encoded_data/encoded_proteins.csv.

3. Merge and Normalize Data (merge_datasets.py)

Scales the encoded features (optionally only the chemical ones).

Merges drug and protein features based on index mappings defined in drug_protein_affinity.csv.

Produces a single dataset stored as encoded_data/merged_dataset.csv.

4. Model Training and Tuning (tune_model.py)

Splits the merged dataset into training and test sets.

Tunes a weighted ElasticNet regression model over multiple chemical/embedding weight ratios.

Saves the best-performing model to models/protein_smile_LR.joblib.

Folder Structure
project_root/
│
├── raw_data/
│   ├── drugs.csv
│   ├── proteins.csv
│   └── drug_protein_affinity.csv
│
├── encoded_data/
│   ├── encoded_drugs.csv
│   ├── encoded_proteins.csv
│   └── merged_dataset.csv
│
├── models/
│   └── protein_smile_LR.joblib
│
├── encode_drugs_data.py
├── encode_proteins_data.py
├── merge_datasets.py
├── tune_model.py
└── README.md

Installation
1. Create a virtual environment
python3 -m venv venv
source venv/bin/activate   # on Linux/Mac
venv\Scripts\activate      # on Windows

2. Install dependencies
pip install -r requirements.txt


If no requirements.txt is available, manually install:

pip install rdkit-pypi torch transformers biopython scikit-learn pandas numpy joblib

Input Data Format
drugs.csv
Column	Description
CID	PubChem Compound ID (optional, dropped)
Canonical_SMILES	Optional canonical representation (dropped)
Isomeric_SMILES	SMILES string used for encoding
Other columns	Metadata, if available
proteins.csv
Column	Description
Accession_Number	Optional ID (dropped)
Gene_Name	Protein or gene identifier
Sequence	Amino acid sequence (single-letter format)
drug_protein_affinity.csv
Column	Description
Protein_Index	Row index in encoded_proteins.csv
Drug_Index	Row index in encoded_drugs.csv
Affinity	Experimental affinity or binding strength
Step-by-Step Usage
1. Encode drugs
python encode_drugs_data.py


Loads raw_data/drugs.csv

Encodes molecules using RDKit + ChemBERTa

Output: encoded_data/encoded_drugs.csv

2. Encode proteins
python encode_proteins_data.py


Loads raw_data/proteins.csv

Computes physicochemical features + ESM2 embeddings

Output: encoded_data/encoded_proteins.csv

3. Merge datasets
python merge_datasets.py


Loads encoded drugs and proteins

Scales features and merges them with affinity data

Output: encoded_data/merged_dataset.csv

4. Tune and train model
python tune_model.py


Loads merged dataset

Tunes ElasticNet over different chemical/embedding weight ratios

Saves model to models/protein_smile_LR.joblib

Writes summary in model_summary

Model Details
Parameter	Description
Model Type	ElasticNetCV
Cross-validation	5-fold
Regularization ratios	[0.1, 0.3, 0.5, 0.7, 0.9]
Iterations	20,000
Evaluation Metric	RMSE (Root Mean Square Deviation)

The model predicts the binding affinity between each drug–protein pair using:

Chemical descriptors + molecular embeddings (ChemBERTa)

Protein physicochemical descriptors + embeddings (ESM2)

Parallelization

Both encoding scripts (encode_drugs_data.py and encode_proteins_data.py) use multiprocessing to accelerate feature computation.

The number of processes is controlled via the variable num_threads (default: 6).

Each process initializes its own model instance to avoid GPU/CPU conflicts.

Output Files
File	Description
encoded_drugs.csv	Drug descriptors + ChemBERTa embeddings
encoded_proteins.csv	Protein descriptors + ESM2 embeddings
merged_dataset.csv	Combined dataset for modeling
protein_smile_LR.joblib	Trained ElasticNet model
model_summary	Best RMSD and optimal chemical weight
Reproducibility Notes

The random seed (random_state=42) ensures reproducibility.

Model outputs and encodings are deterministic given the same inputs.

Feature scaling is standardized (mean = 0, std = 1) before training.

Citation and References

ChemBERTa: Chithrananda, S., Grand, G., & Ramsundar, B. (2020). ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction.

ESM2: Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model.
