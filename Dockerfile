FROM python:3.10-slim

# Install OS-level dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    libeigen3-dev \
    libxrender1 \
    libsm6 \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy FastAPI backend code
COPY ./app /app

# Install Python packages (only FastAPI-related)
RUN pip install --no-cache-dir \
    fastapi \
    uvicorn \
    rdkit \
    torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu

# Expose FastAPI port
EXPOSE 8000

# Start FastAPI server
CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8000"]
