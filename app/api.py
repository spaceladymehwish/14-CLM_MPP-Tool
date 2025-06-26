import os
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import torch
from constants import descriptor_list, file_dir, model_file_name
from model import ImprovedMolecularNN

app = FastAPI()


model_file_path = os.path.join(file_dir, model_file_name)
input_dim = 140

model = ImprovedMolecularNN(input_dim)
model.load_state_dict(torch.load(model_file_path, map_location='cpu'))
model.eval()

descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_list)


class MoleculeRequest(BaseModel):
    smiles: str


class PredictionResponse(BaseModel):
    prediction: str
    confidence: float


@app.post("/predict", response_model=PredictionResponse)
def predict_molecule(req: MoleculeRequest):
    mol = Chem.MolFromSmiles(req.smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")

    desc_values = descriptor_calc.CalcDescriptors(mol)
    x_tensor = torch.tensor([desc_values], dtype=torch.float32)

    with torch.no_grad():
        prob = model(x_tensor).item()
        pred = "Active" if prob >= 0.5 else "Inactive"

    return PredictionResponse(prediction=pred, confidence=round(prob, 3))
