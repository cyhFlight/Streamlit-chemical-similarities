# main file to run
import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw
from utils import calculate_tanimoto_from_list, calculate_tanimoto

# UI parameters
ketcher_height = 500
st.set_page_config(layout="wide")
st.title("Ketcher Molecule Editor")

# Input SMILES
default_mol = "C[N+]1=CC=C(/C2=C3\C=CC(=N3)/C(C3=CC=CC(C(N)=O)=C3)=C3/C=C/C(=C(\C4=CC=[N+](C)C=C4)C4=N/C(=C(/C5=CC=CC(C(N)=O)=C5)C5=CC=C2N5)C=C4)N3)C=C1"
molecule_str = st.text_input("Input molecule SMILES", default_mol)
smiles_code = st_ketcher(molecule_str, height=ketcher_height)
st.markdown("Click :blue-background[Apply] to obtain the SMILES of the molecule")
st.markdown(f"SMILES code: ``{smiles_code}``")

# Draw molecule
if molecule_str:
    try:
        mol = Chem.MolFromSmiles(molecule_str)
        if mol:
            st.image(Draw.MolToImage(mol), caption="Molecule structure", use_container_width=False)
    except Exception as e:
        st.error(f"Error: {e}")

st.divider()

# Parameters for Morgan Fingerprints
st.title("Calculate Tanimoto similarities using Morgan Fingerprints")
default_radius = 2
default_nbits = 2048

col1, col2 = st.columns(2)
with col1:
    mfp_radius = st.number_input("Morgan Fingerprint Radius", min_value=1, max_value=5, value=default_radius)
with col2:
    mfp_nbits = st.number_input("Fingerprint Bit Vector Size", min_value=512, max_value=4096, value=default_nbits, step=512)

# Reference drug list
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

# Display reference drug images
st.subheader("Chemical Structures of Reference Drugs")
cols = st.columns(5)
drug_images = []
for drug_name, drug_smiles in drug_list_10:
    mol = Chem.MolFromSmiles(drug_smiles)
    if mol:
        img = Draw.MolToImage(mol)
        drug_images.append((drug_name, img))
for idx, (name, img) in enumerate(drug_images):
    with cols[idx % 5]:
        st.image(img, caption=name, use_container_width=False)

# Calculate similarities
st.subheader("Tanimoto Similarity to Common Drugs")
if st.button('Calculate', type='primary'):
    styled_df = calculate_tanimoto_from_list(molecule_str, drug_list=drug_list_10, radius=mfp_radius, nbits=mfp_nbits)
    st.write(styled_df)
