# utils.py
from rdkit import Chem
from rdkit.Chem import Draw, DataStructs, rdFingerprintGenerator
import pandas as pd

def calculate_tanimoto(smiles1, smiles2, radius=2, nbits=2048):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)

    if mol1 and mol2:
        fp1 = morgan_generator.GetFingerprint(mol1)
        fp2 = morgan_generator.GetFingerprint(mol2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    return None

def calculate_tanimoto_from_list(molecule_str, drug_list, radius=2, nbits=2048):
    similarities = []
    for drug_name, drug_smiles in drug_list:
        similarity = calculate_tanimoto(molecule_str, drug_smiles, radius, nbits)
        if similarity is not None:
            similarities.append((drug_name, drug_smiles, similarity))

    similarities.sort(key=lambda x: x[2], reverse=True)
    df = pd.DataFrame(similarities, columns=["Drug Name", "SMILES", "Tanimoto Similarity"])
    styled_df = df.style.background_gradient(cmap="coolwarm", subset=["Tanimoto Similarity"])
    return styled_df

def compute_top_similarities(query_smiles, df, radius=2, nbits=2048):
    results = []
    for _, row in df.iterrows():
        sim = calculate_tanimoto(query_smiles, row['smiles'], radius, nbits)
        if sim is not None:
            results.append((sim, row['substance'], row['smiles'], row['class1'], row['class3']))
    results.sort(key=lambda x: x[0], reverse=True)
    return pd.DataFrame(results[:], columns=["Similarity", "Drug Name", "SMILES", "Class 1", "Class 3"])

drug_list_10 = [
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
    ("Acetaminophen", "CC(=O)NC1=CC=C(O)C=C1"),
    ("Penicillin V", "CC1(C)S[C@H](N2C=NC3=C2C(=O)NC(=O)N3)C(NC1=O)C(=O)O"),
    ("Amoxicillin", "CC1(C)S[C@H](N2C=NC3=C2C(=O)NC(=O)N3)C(NC1=O)C(=O)OCCN"),
    ("Metformin", "CNC(=N)NC(=N)NCCN"),
    ("Lisinopril", "CC(C)C[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](N)CC2=CC=CC=C2)C(=O)O"),
    ("Simvastatin", "CC(C)C1CCC2C(CCC3C2C1CC4C3(CCC(C4)O)C(=O)OCC=C(C)C)C"),
    ("Warfarin", "CC(=O)CC1=C(O)C2=C(C=CC=C2)C(=O)C1")
]
