import os
import dotenv
import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Draw
dotenv.load_dotenv(".env")
url = os.getenv("BACKEND_URL")

st.title("Real-Time Molecular Property Predictor")


smiles_input = st.text_input("Enter a SMILES string:")

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    if mol:
        st.image(Draw.MolToImage(mol, size=(300, 300)))

        try:
            response = requests.post(url, json={"smiles": smiles_input})
            if response.status_code == 200:
                result = response.json()
                st.markdown(f"### Prediction: {result['prediction']}")
                st.markdown(f"**Confidence (probability): {result['confidence']:.3f}**")
            else:
                st.error(f"Error from server: {response.json().get('error', 'Unknown error')}")
        except Exception as e:
            st.error(f"Failed to connect to Flask server: {e}")
    else:
        st.error("Invalid SMILES string.")
