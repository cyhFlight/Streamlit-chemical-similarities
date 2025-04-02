import streamlit as st
from streamlit_ketcher import st_ketcher
# import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw


# UI parameters
ketcher_height = 800
st.set_page_config(layout="wide")




st.title("Ketcher Molecule Editor")

## Ketcher HTML/JavaScript embedding not working
## https://unpkg.com/ketcher@latest/dist/ketcher.min.js not found

# Streamlit_ketcher module
default_mol = "C[N+]1=CC=C(/C2=C3\C=CC(=N3)/C(C3=CC=CC(C(N)=O)=C3)=C3/C=C/C(=C(\C4=CC=[N+](C)C=C4)C4=N/C(=C(/C5=CC=CC(C(N)=O)=C5)C5=CC=C2N5)C=C4)N3)C=C1"
molecule_str = st.text_input("Molecule", default_mol)
smiles_code = st_ketcher(molecule_str, height=ketcher_height)
st.markdown(f"SMILES code: ``{smiles_code}``")

print(smiles_code)

# Use RDkit to draw
if molecule_str:
    try:
        mol = Chem.MolFromSmiles(molecule_str)
        if mol:
            st.image(Draw.MolToImage(mol), caption="Molecule structure", use_container_width=True)
    except Exception as e:
        st.error(f"Error: {e}")
