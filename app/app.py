import sys
import os
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import torch
from src.CLM import descriptor_calc, model

st.title("Real-Time Molecular Property Predictor")
smiles_input = st.text_input("Enter a SMILES string:")

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    if mol:
        # Show molecule
        st.image(Draw.MolToImage(mol, size=(300, 300)))

        # Compute descriptors
        desc_values = descriptor_calc.CalcDescriptors(mol)
        x_tensor = torch.tensor([desc_values], dtype=torch.float32)

        with torch.no_grad():
            prob = model(x_tensor).item()
            pred = 1 if prob >= 0.5 else 0

        st.markdown(f"### Prediction: {'Active' if pred == 1 else 'Inactive'}")
        st.markdown(f"**Confidence (probability): {prob:.3f}**")

    else:
        st.error("Invalid SMILES string.")
