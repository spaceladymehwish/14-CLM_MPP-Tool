FROM continuumio/miniconda3

# Create conda env with RDKit
RUN conda create -y -n rdkit-env -c conda-forge python=3.10 rdkit \
    && conda clean -afy

# Activate env
SHELL ["conda", "run", "-n", "rdkit-env", "/bin/bash", "-c"]

# Set working directory
WORKDIR /app

# Copy and install pip requirements
COPY requirements.txt .
RUN conda run -n rdkit-env pip install --upgrade pip && conda run -n rdkit-env pip install -r requirements.txt

# Copy app files
COPY . .

EXPOSE 8501

# Run Streamlit via conda env
CMD ["conda", "run", "-n", "rdkit-env", "streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]