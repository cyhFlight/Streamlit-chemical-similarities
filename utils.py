# functions and utils
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFingerprintGenerator
import pandas as pd


def calculate_tanimoto(smiles1, smiles2, radius=2, nbits=2048):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)

    if mol1 and mol2:
        fp1 = morgan_generator.GetFingerprint(mol1)
        fp2 = morgan_generator.GetFingerprint(mol2)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        return similarity
    else:
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
