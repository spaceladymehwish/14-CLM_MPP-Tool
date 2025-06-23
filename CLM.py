import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import torch
import torch.nn as nn
import numpy as np

# Define your exact 140 descriptors here
your_140_descriptor_list = [
    'BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v',
    'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3',
    'LabuteASA', 'MolLogP', 'MolMR', 'MolWt', 'TPSA', 'HeavyAtomCount', 'NumHAcceptors',
    'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumValenceElectrons',
    'RingCount', 'FractionCSP3', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
    'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
    'NumAromaticHeterocycles', 'NumAromaticRings', 'NumSaturatedCarbocycles',
    'NumSaturatedHeterocycles', 'NumSaturatedRings', 'EState_VSA1', 'EState_VSA2', 'EState_VSA3',
    'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9',
    'EState_VSA10', 'EState_VSA11', 'EState_VSA12', 'EState_VSA13', 'EState_VSA14', 'PEOE_VSA1',
    'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8',
    'PEOE_VSA9', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
    'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7',
    'SMR_VSA8', 'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5',
    'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10', 'SlogP_VSA11',
    'SlogP_VSA12', 'SlogP_VSA13', 'SlogP_VSA14', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert',
    'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_S',
    'fr_ether', 'fr_furan', 'fr_guanido', 'fr_imidazole', 'fr_imide', 'fr_amine', 'fr_lactam',
    'fr_methoxy', 'fr_morpholine', 'fr_nitro', 'fr_oxazole', 'fr_piperidine', 'fr_pyridine',
    'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_urea',
    'fr_benzene', 'fr_ester', 'fr_amide', 'fr_alkyl_carbamate', 'fr_ketone', 'fr_nitrile',
    'fr_thiophene', 'fr_azo', 'fr_azide', 'fr_halogen', 'fr_hdrzone', 'fr_isocyan', 'fr_isothiocyan',
    'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_term_acetylene'
]
descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(your_140_descriptor_list)

# Define the model architecture
class ImprovedMolecularNN(nn.Module):
    def __init__(self, input_dim):
        super(ImprovedMolecularNN, self).__init__()
        self.fc1 = nn.Linear(input_dim, 512)
        self.bn1 = nn.BatchNorm1d(512)
        self.fc2 = nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.fc3 = nn.Linear(256, 128)
        self.bn3 = nn.BatchNorm1d(128)
        self.fc4 = nn.Linear(128, 64)
        self.bn4 = nn.BatchNorm1d(64)
        self.fc5 = nn.Linear(64, 1)

        self.leaky_relu = nn.LeakyReLU(0.1)
        self.dropout = nn.Dropout(0.4)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.leaky_relu(self.bn1(self.fc1(x)))
        x = self.dropout(x)
        x = self.leaky_relu(self.bn2(self.fc2(x)))
        x = self.dropout(x)
        x = self.leaky_relu(self.bn3(self.fc3(x)))
        x = self.dropout(x)
        x = self.leaky_relu(self.bn4(self.fc4(x)))
        x = self.sigmoid(self.fc5(x))
        return x

# Load trained model
input_dim = 140
model = ImprovedMolecularNN(input_dim)
model = torch.load("2D_rdkit_best_nn_model_full.pth", map_location='cpu', weights_only=False)
model.eval()  # Set model to evaluation mode

# Streamlit UI
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
