import streamlit as st
from streamlit_ketcher import st_ketcher
# import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw


# UI parameters
ketcher_height = 500
st.set_page_config(layout="wide")
st.title("Ketcher Molecule Editor")

## Ketcher HTML/JavaScript embedding not working
## https://unpkg.com/ketcher@latest/dist/ketcher.min.js not found

# Streamlit_ketcher module to obtain the molecule of interest
default_mol = "C[N+]1=CC=C(/C2=C3\C=CC(=N3)/C(C3=CC=CC(C(N)=O)=C3)=C3/C=C/C(=C(\C4=CC=[N+](C)C=C4)C4=N/C(=C(/C5=CC=CC(C(N)=O)=C5)C5=CC=C2N5)C=C4)N3)C=C1"
molecule_str = st.text_input("Input molecule SMILES", default_mol)
smiles_code = st_ketcher(molecule_str, height=ketcher_height)
st.markdown("Click :blue-background[Apply] to obtain the SMILES of the molecule")
st.markdown(f"SMILES code: ``{smiles_code}``")

# Use RDkit to draw
if molecule_str:
    try:
        mol = Chem.MolFromSmiles(molecule_str)
        if mol:
            st.image(Draw.MolToImage(mol), caption="Molecule structure", use_container_width=False)
    except Exception as e:
        st.error(f"Error: {e}")


###################################
# Compare two molecules

# Input field for the second molecule's SMILES
second_molecule_str = st.text_input("Input second molecule SMILES")

# Adjustable parameters for fingerprint calculation
default_radius = 2
default_nbits = 2048
radius = st.number_input("Morgan Fingerprint Radius", min_value=1, max_value=5, value=default_radius)
nbits = st.number_input("Fingerprint Bit Vector Size", min_value=512, max_value=4096, value=default_nbits, step=512)

# Function to calculate Tanimoto similarity between two SMILES strings
def calculate_tanimoto(smiles1, smiles2, radius=default_radius, nbits=default_nbits):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 and mol2:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=nbits)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=nbits)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        return similarity
    else:
        return None

# Display the Tanimoto similarity between the two input molecules
if molecule_str and second_molecule_str:
    similarity = calculate_tanimoto(molecule_str, second_molecule_str, radius, nbits)
    if similarity is not None:
        st.markdown(f"Tanimoto similarity between the two molecules: **{similarity:.2f}**")
    else:
        st.error("Invalid SMILES input for one or both molecules.")



###############################
# Compare with ten drug molecules


