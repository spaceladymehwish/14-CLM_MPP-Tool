import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Draw

st.title("Real-Time Molecular Property Predictor")
smiles_input = st.text_input("Enter a SMILES string:")

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    if mol:
        st.image(Draw.MolToImage(mol, size=(300, 300)))

        response = requests.post(BACKEND_URL, json={"smiles": smiles_input})
        if response.status_code == 200:
            result = response.json()
            st.markdown(f"### Prediction: {result['prediction']}")
            st.markdown(f"**Confidence (probability): {result['confidence']:.3f}**")
        else:
            st.error("Invalid SMILES or server error.")
    else:
        st.error("Invalid SMILES string.")
