# Use conda/miniconda base image
FROM continuumio/miniconda3

# Create conda environment with Python 3.10 and RDKit
RUN conda create -n myenv -c conda-forge python=3.10 rdkit=2022.9.5

# Activate conda environment for all subsequent commands
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Copy your app code into the container
WORKDIR /app
COPY . /app

# Install pip dependencies inside conda env
RUN conda run -n myenv pip install --no-cache-dir -r requirements.txt

# Expose port (Streamlit default)
EXPOSE 8501

# Run Streamlit in the activated conda env
CMD ["conda", "run", "-n", "myenv", "streamlit", "run", "CLM.py", "--server.port=8501", "--server.address=0.0.0.0"]
