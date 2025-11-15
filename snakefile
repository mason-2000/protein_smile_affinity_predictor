rule all:
    input: 
        "model/protein_smile_LR.joblib",
        "model/model_summary.txt"

rule encode_drugs:
    input: 
        "raw_data/drugs.csv"   
    output: 
        "encoded_data/encoded_drugs.csv"
    shell:
        "python3 encode_drugs_data.py"

rule encode_proteins:
    input: 
        "raw_data/proteins.csv"   
    output: 
        "encoded_data/encoded_proteins.csv"
    shell:
        "python3 encode_proteins_data.py"

rule merge_data:
    input:
        "encoded_data/encoded_proteins.csv",
        "encoded_data/encoded_drugs.csv",
        "encoded_data/drug_protein_affinity.csv"
    output:
        "encoded_data/merged_dataset.csv"
    shell:
        "python3 merge_datasets.py"

rule train_model:
    input:
        "encoded_data/merged_dataset.csv"
    output:
        "model/protein_smile_LR.joblib",
        "model/model_summary.txt"
    shell:
        "python3 train_tune_model.py"
