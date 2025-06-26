Molecular Property Predictor: 2D Descriptor-Based Deep Learning Classifier for Antibiotics
The molecular property predictor is a deep learning model built to classify antibiotic molecules as either active or inactive, using features derived from their 2D chemical structure. It plays a critical role in virtual screening by enabling rapid and automated evaluation of candidate molecules.

Feature Computation: 2D RDKit Descriptors
Each molecule is represented as a SMILES string, from which 2D molecular descriptors are computed using RDKit, an open-source cheminformatics library. These descriptors include:

Topological descriptors (e.g., molecular weight, number of rings, rotatable bonds)

Electro-topological state indices

Molecular fragment counts

Van der Waals surface area (VSA) features

Partial charge and polarity-related features

Substructure-based fingerprints

A total of 140 numerical features are extracted per molecule, summarizing its physicochemical, structural, and functional properties.

Deep Learning Model
The model is a feedforward neural network trained to predict binary activity labels based on the 2D descriptors:

Input: A 140-dimensional vector of 2D molecular features.

Architecture:

Multiple fully connected layers

Nonlinear activations (e.g., ReLU or LeakyReLU)

Regularization techniques like Dropout and Batch Normalization

Early stopping to prevent overfitting

Output: A single value between 0 and 1, interpreted as the probability of the molecule being active.

Use Case
This property predictor is especially useful in antibiotic discovery workflows and molecule generation pipelines (e.g., Monte Carlo Tree Search), where evaluating thousands of molecules for potential activity is computationally demanding.

Given a new SMILES string:

Its 2D descriptors are computed.

The descriptor vector is passed through the trained neural network.

The model outputs an activity prediction: active (1) or inactive (0).
