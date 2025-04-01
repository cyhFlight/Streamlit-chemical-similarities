import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw

# Streamlit 标题
st.title("Ketcher 分子编辑器")

# Ketcher 的 HTML/JavaScript 嵌入
ketcher_html = """
<!DOCTYPE html>
<html>
<head>
    <script src="https://unpkg.com/ketcher@latest/dist/ketcher.min.js"></script>
    <style>
        #ketcher-container {
            width: 100%;
            height: 400px;
        }
    </style>
</head>
<body>
    <div id="ketcher-container"></div>
    <button onclick="getSmiles()">获取 SMILES</button>
    <input type="text" id="smilesOutput" readonly>
    <script>
        let ketcher = Ketcher.create('#ketcher-container');
        function getSmiles() {
            ketcher.getSmiles()
                .then(smiles => {
                    document.getElementById('smilesOutput').value = smiles;
                    fetch('/?smiles=' + encodeURIComponent(smiles), {method: 'POST'});
                })
                .catch(error => console.error(error));
        }
    </script>
</body>
</html>
"""

# 在 Streamlit 中渲染 Ketcher
components.html(ketcher_html, height=500)

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
