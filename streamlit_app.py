# streamlit_app.py
import mols2grid
import streamlit as st
import pandas as pd
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw
from utils import calculate_tanimoto_from_list, compute_top_similarities, drug_list_10

# UI parameters
ketcher_height = 500
st.set_page_config(layout="wide")
st.title("Ketcher Molecule Editor")

# Default molecule
default_mol = "C[N+]1=CC=C(/C2=C3\C=CC(=N3)/C(C3=CC=CC(C(N)=O)=C3)=C3/C=C/C(=C(\C4=CC=[N+](C)C=C4)C4=N/C(=C(/C5=CC=CC(C(N)=O)=C5)C5=CC=C2N5)C=C4)N3)C=C1"

# Initialize session state
if 'smiles' not in st.session_state:
    st.session_state.smiles = default_mol
if 'text_input' not in st.session_state:
    st.session_state.text_input = default_mol

def draw_mol():
    st.session_state.smiles = st.session_state.text_input
    st.markdown("Your quest molecule:")
    st.markdown(f"SMILES code: ``{st.session_state.smiles}``")
    # RDKit molecule drawing
    try:
        mol = Chem.MolFromSmiles(st.session_state.smiles)
        if mol:
            st.image(Draw.MolToImage(mol), caption="Molecule structure", use_container_width=False)
    except Exception as e:
        st.error(f"Error: {e}")

# Text input box
st.session_state.text_input = st.text_input("Input molecule SMILES", st.session_state.text_input)

# Ketcher editor
st.session_state.text_input = st_ketcher(st.session_state.text_input, height=ketcher_height)

# Sync mol smiles
st.markdown("Click :blue-background[Apply] to apply the molecule you draw")
if st.button("Confirm Molecule", type="primary"):
    draw_mol()

# Load FDA drugs data from csv
st.divider()
st.title("Load FDA approved drugs")
df_fda = pd.read_csv("data/FDA_drugs_smiles.csv")
st.markdown("The dimension of the loaded data frame is:")
st.write(df_fda.shape)
st.write(df_fda)

st.divider()
# Adjustable parameters for fingerprint calculation
default_radius = 2
default_nbits = 2048
st.title("Calculate Tanimoto similarities using Morgan Fingerprints")

col1, col2 = st.columns(2)
with col1:
    mfp_radius = st.number_input("Morgan Fingerprint Radius", min_value=1, max_value=5, value=default_radius)
with col2:
    mfp_nbits = st.number_input("Fingerprint Bit Vector Size", min_value=512, max_value=4096, value=default_nbits, step=512)

# Show quest molecule again for easy comparison
draw_mol()
# st.subheader("Tanimoto Similarity to FDA Approved Drugs")
if st.button('Search Similar FDA Approved Drugs', type="primary"):
    top_sim_df = compute_top_similarities(st.session_state.smiles, df_fda, radius=mfp_radius, nbits=mfp_nbits)
    st.subheader("FDA Approved Drugs Sorted by Similarity")
    st.dataframe(top_sim_df.style.background_gradient(cmap="coolwarm", subset=["Similarity"])) # with background gradient

    raw_html = mols2grid.display(top_sim_df,
                            subset=["img", "Drug Name"],
                            mapping={"SMILES": "SMILES", "Drug Name": "Name"})._repr_html_()
    st.components.v1.html(raw_html, width=1200, height=650, scrolling=False)


#### The below is just demo, in expander format
# Display images of reference drug molecules
with st.expander("Chemical Structures of 10 Reference Drugs (Deprecated)"):
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

    # Tanimoto similarity table
    st.subheader("Tanimoto Similarity to Common Drugs")
    if st.button('Calculate Similarity', type="primary"):
        styled_df = calculate_tanimoto_from_list(
            st.session_state.smiles, drug_list=drug_list_10,
            radius=mfp_radius, nbits=mfp_nbits
        )
        st.write(styled_df)
