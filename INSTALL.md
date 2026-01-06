# Installation Guide

This guide provides detailed instructions for setting up the computational environment needed to run the glacier flow simulations.

---

## System Requirements

### Minimum Requirements
- **OS**: Linux (Ubuntu 20.04+), macOS (10.15+), or Windows (WSL2)
- **RAM**: 4 GB (8 GB recommended)
- **Disk Space**: 2 GB free space
- **Python**: 3.8 or higher

### Recommended for Large-Scale Simulations
- **RAM**: 16 GB+
- **CPU**: Multi-core processor (4+ cores)
- **MPI**: For parallel computing capabilities

---

## Installation Methods

### Method 1: Conda (Recommended)

This is the easiest and most reliable method for installing FEniCSx and all dependencies.

#### Step 1: Install Miniconda/Anaconda

If you don't have conda installed:

```bash
# Download Miniconda (Linux/macOS)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts and restart terminal
```

For macOS with Apple Silicon (M1/M2):
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

#### Step 2: Create FEniCSx Environment

```bash
# Create a new conda environment
conda create -n fenicsx-env python=3.10
conda activate fenicsx-env

# Install FEniCSx and dependencies
conda install -c conda-forge fenics-dolfinx mpich pyvista
```

#### Step 3: Install Additional Python Packages

```bash
pip install matplotlib numpy scipy
```

#### Step 4: Verify Installation

```bash
python -c "import dolfinx; print(f'DOLFINx version: {dolfinx.__version__}')"
python -c "import ufl; print(f'UFL version: {ufl.__version__}')"
```

Expected output:
```
DOLFINx version: 0.8.0 (or later)
UFL version: 2024.1.0 (or similar)
```

---

### Method 2: Docker (Cross-Platform)

For users who want a completely isolated environment or are using Windows.

#### Step 1: Install Docker

Download and install Docker Desktop from: https://www.docker.com/products/docker-desktop

#### Step 2: Pull FEniCSx Image

```bash
docker pull dolfinx/dolfinx:stable
```

#### Step 3: Run Container with Project

```bash
# Navigate to project directory
cd /path/to/picard-fem-glacier-flow

# Run container with volume mount
docker run -ti -v $(pwd):/home/shared dolfinx/dolfinx:stable
```

Inside the container:
```bash
cd /home/shared
python src/main.py
```

---

### Method 3: From Source (Advanced)

For users who need the latest development version or custom builds.

**⚠️ Warning**: This method is complex and may require solving dependency conflicts.

See official documentation: https://docs.fenicsproject.org/dolfinx/main/python/installation.html

---

## Project-Specific Setup

### Clone Repository

```bash
git clone https://github.com/yourusername/picard-fem-glacier-flow.git
cd picard-fem-glacier-flow
```

### Verify Data Files

Ensure the mesh data is present:

```bash
ls -lh data/arolla.*
```

You should see:
```
arolla.h5      (≈77 KB)
arolla.xdmf    (≈1 KB)
```

### Test Installation

Run a quick test to ensure everything works:

```bash
# Activate conda environment (if using conda)
conda activate fenicsx-env

# Run linear solver
python src/main.py
```

If successful, you should see:
```
Mesh loaded: ...
[Solver output]
Plot saved to plots/velocity_field.png
```

---

## Troubleshooting

### Issue: `ModuleNotFoundError: No module named 'dolfinx'`

**Solution**: Ensure you've activated the correct conda environment:
```bash
conda activate fenicsx-env
```

### Issue: `FileNotFoundError: data/arolla.xdmf`

**Solution**: Run scripts from the project root directory:
```bash
cd picard-fem-glacier-flow
python src/main.py  # Not: cd src && python main.py
```

### Issue: MPI Errors

**Solution**: Ensure MPI is properly installed:
```bash
# Test MPI installation
mpirun --version

# If issues persist, reinstall via conda
conda install -c conda-forge mpich --force-reinstall
```

### Issue: Plotting Fails on Headless Systems

**Solution**: Use matplotlib's Agg backend:
```python
# Add to top of script
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```

### Issue: Memory Errors with Large Meshes

**Solution**: 
1. Reduce mesh resolution in preprocessing
2. Use iterative solvers instead of direct LU:
   ```python
   petsc_options={"ksp_type": "cg", "pc_type": "hypre"}
   ```

---

## Optional: ParaView for Advanced Visualization

ParaView is recommended for exploring 3D velocity fields and creating publication-quality visualizations.

### Installation

**Linux (Ubuntu/Debian)**:
```bash
sudo apt update
sudo apt install paraview
```

**macOS** (via Homebrew):
```bash
brew install --cask paraview
```

**Windows or from Official Site**:
Download from: https://www.paraview.org/download/

### Usage

```bash
# Open VTU/VTX files
paraview glacier_solution_model.vtu
```

---

## Environments for Different Use Cases

### Quick Testing (Minimal)
```bash
conda create -n fem-minimal python=3.10
conda activate fem-minimal
conda install -c conda-forge fenics-dolfinx matplotlib
```

### Full Development (Recommended)
```bash
conda create -n fem-dev python=3.10
conda activate fem-dev
conda install -c conda-forge fenics-dolfinx mpich pyvista jupyter
pip install matplotlib scipy pandas black pytest
```

### High-Performance Computing
For cluster/HPC environments, consult your system administrator. Typically:
```bash
module load python/3.10
module load openmpi
pip install --user fenics-dolfinx
```

---

## Next Steps

After successful installation:

1. **Run Examples**: Start with `examples/poisson2D.py` to verify FEM basics
2. **Linear Model**: Execute `python src/main.py`
3. **Non-Linear Model**: Run `python src/nonlinear_solver.py`
4. **Read Documentation**: Review `solutions.pdf` for theoretical background

---

## Getting Help

- **FEniCSx Documentation**: https://docs.fenicsproject.org/dolfinx/
- **FEniCS Discourse**: https://fenicsproject.discourse.group/
- **Project Issues**: [GitHub Issues page for this repo]

---

*Last Updated: January 2026*
