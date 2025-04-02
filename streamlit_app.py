import streamlit as st
from streamlit_ketcher import st_ketcher
# import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw

# Streamlit 标题
st.title("Ketcher 分子编辑器")

# Ketcher 的 HTML/JavaScript 嵌入
## not working https://unpkg.com/ketcher@latest/dist/ketcher.min.js not found

molecule = st.text_input("Molecule", DEFAULT_MOL)
smile_code = st_ketcher(molecule)
st.markdown(f"Smile code: ``{smile_code}``")

# 从 URL 参数获取 SMILES
smiles = st.text_input("输入在 Ketcher 中创建的分子 SMILES:", "")

# 显示分子结构
if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            st.image(Draw.MolToImage(mol), caption="分子结构", use_column_width=True)
    except Exception as e:
        st.error(f"错误: {e}")
