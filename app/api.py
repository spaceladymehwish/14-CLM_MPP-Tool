import os
from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import torch
from constants import descriptor_list, file_dir, model_file_name
from model import ImprovedMolecularNN

app = Flask(__name__)

model_file_path = os.path.join(file_dir, model_file_name)
input_dim = 140

model = ImprovedMolecularNN(input_dim)
model.load_state_dict(torch.load(model_file_path, map_location='cpu'))
model.eval()

descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_list)


@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    smiles = data.get("smiles")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    desc_values = descriptor_calc.CalcDescriptors(mol)
    x_tensor = torch.tensor([desc_values], dtype=torch.float32)

    with torch.no_grad():
        prob = model(x_tensor).item()
        pred = "Active" if prob >= 0.5 else "Inactive"

    return jsonify({
        "prediction": pred,
        "confidence": round(prob, 3)
    })
