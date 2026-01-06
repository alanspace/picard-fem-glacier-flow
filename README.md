# Finite Element Methods for Non-Newtonian Glacier Flow Dynamics

**A Computational Study of the Haut Glacier d'Arolla using FEniCSx and Picard Iteration**

---

## 📋 Overview

This repository presents a comprehensive computational study of glacier flow dynamics using Finite Element Methods (FEM). The project implements the **First-Order Stokes approximation** to simulate the Haut Glacier d'Arolla in Switzerland, solving both linear and non-linear formulations of the governing equations.

### Key Features

- **Full FEM Implementation**: Complete weak formulation derivation and implementation using FEniCSx
- **Non-Linear Solver**: Picard iteration method for handling strain-rate dependent viscosity  
- **Mathematical Rigor**: Includes proofs of Galerkin orthogonality and best approximation properties
- **Realistic Glacier Simulation**: Applied to real-world glacier geometry data
- **Comprehensive Analysis**: Systematic study of regularization parameters and convergence behavior

---

## 🎯 Project Objectives

This work addresses fundamental questions in computational glaciology:

1. **Mathematical Analysis**: Classification and theoretical properties of glacier flow PDEs
2. **Numerical Methods**: Development of efficient FEM solvers for elliptic systems
3. **Non-Linear Dynamics**: Implementation of iterative schemes for strain-rate dependent rheology
4. **Physical Validation**: Comparison between linear and non-linear ice flow models

---

## 🔬 Physical & Mathematical Background

### Governing Equations

Glacier flow is modeled using the **First-Order (Shallow Ice) Stokes approximation**:

```
2∂ₓ(η∂ₓu) + ½∂ᵤ(η∂ᵤu) = ρg∂ₓh
```

Where:
- `η`: Effective viscosity (Pa·s)
- `u`: Horizontal velocity (m/s)
- `ρ`: Ice density (910 kg/m³)
- `g`: Gravitational acceleration (9.81 m/s²)
- `h`: Surface elevation (m)

### Non-Newtonian Rheology

Ice exhibits **power-law (Glen's law)** viscosity:

```
η(ε̇) = ½A^(-1/n) ε̇ₑ^((1-n)/n)
```

With `n=3` (Glen's flow law exponent) and `A` the flow parameter.

### PDE Classification

- **Type**: Second-order elliptic PDE
- **Discriminant**: Δ = B² - AC < 0
- **Physical Meaning**: Equilibrium state with instantaneous force balance

---

## 🛠️ Implementation Details

### Technology Stack

- **Python 3.x**
- **FEniCSx**: Modern finite element framework
- **NumPy**: Numerical computations
- **Matplotlib**: Visualization and plotting
- **ParaView** (optional): Advanced 3D visualization

### Computational Approach

#### 1. Weak Formulation

The variational form seeks `u ∈ V` such that ∀v ∈ V₀:

```
∫_Ω (2η∂ₓu∂ₓv + ½η∂ᵤu∂ᵤv) dΩ = -∫_Ω ρg∂ₓh·v dΩ
```

#### 2. Discretization

- **Element Type**: Continuous Piecewise Linear (P1) elements
- **Function Space**: H¹(Ω) with Dirichlet boundary conditions
- **Mesh**: Unstructured triangular mesh from real glacier geometry

#### 3. Non-Linear Solver

**Picard Iteration** (Fixed-Point Method):

```python
1. Initialize: u₀ = 0
2. For k = 1, 2, ... until convergence:
   a. Compute η_{k-1} from u_{k-1}
   b. Solve linear system with fixed η_{k-1} → u_k
   c. Check: ||u_k - u_{k-1}|| < tolerance
```

---

## 📁 Repository Structure

```
picard-fem-glacier-flow/
├── README.md                    # This file
├── solutions.pdf                # Complete technical report
├── solutions.tex                # LaTeX source for report
│
├── src/                         # Source code
│   ├── main.py                  # Main simulation driver
│   ├── nonlinear_solver.py      # Picard iteration implementation
│   ├── helpfunctions.py         # Utility functions for FEM
│   ├── compare_profiles.py      # Analysis and comparison tools
│   └── generate_plot_analytical.py  # Plotting utilities
│
├── data/                        # Input data
│   ├── arolla.h5                # Glacier geometry (HDF5 format)
│   └── arolla.xdmf              # XDMF mesh descriptor
│
├── plots/                       # Generated visualizations
│   ├── velocity_field.png       # Linear model results
│   ├── nonlinear_velocity.png   # Non-linear model results
│   └── profile_comparison.png   # Linear vs non-linear comparison
│
├── examples/                    # Tutorial examples
│   ├── poisson2D.py             # 2D Poisson equation example
│   ├── poisson_project.py       # Basic FEM demonstration
│   └── 2D_Hat.py                # Hat function illustration
│
└── docs/                        # Additional documentation
    ├── Project2025.pdf          # Original project description
    └── references/              # Reference papers
        └── A cut finite element method for...pdf
```

---

## 🚀 Getting Started

### Prerequisites

```bash
# Install FEniCSx (via conda - recommended)
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista

# Or using pip (may require additional system dependencies)
pip install fenics-dolfinx
```

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/picard-fem-glacier-flow.git
cd picard-fem-glacier-flow

# Verify installation
python -c "import dolfinx; print(dolfinx.__version__)"
```

### Running the Simulations

#### Linear Model

```bash
cd src
python main.py
```

This will:
- Load the Arolla glacier mesh
- Solve the linear Stokes equation (η = constant)
- Generate velocity field visualizations
- Export results in VTU format for ParaView

#### Non-Linear Model

```bash
cd src
python nonlinear_solver.py --epsilon 1e-10 --max-iter 50 --tolerance 1e-6
```

**Command-line arguments:**
- `--epsilon`: Regularization parameter (default: 1e-10)
- `--max-iter`: Maximum Picard iterations (default: 50)
- `--tolerance`: Convergence tolerance (default: 1e-6)

### Visualization

```bash
# Generate comparison plots
python src/compare_profiles.py

# View results in ParaView (if installed)
paraview glacier_solution_model.vtu
```

---

## 📊 Key Results

### Convergence Analysis

| Regularization ε | Iterations to Converge | Physical Interpretation |
|------------------|------------------------|-------------------------|
| 10⁻¹             | 3                      | Nearly linear regime |
| 10⁻⁵             | 3                      | Moderate non-linearity |
| 10⁻¹⁰            | 10                     | Full non-linear behavior |

### Linear vs Non-Linear Comparison

**Key Observations:**

1. **Velocity Magnitude**: Non-linear model produces velocities ~1000× larger than linear model
   - Linear: ~10⁻⁶ m/s
   - Non-linear: ~10⁻³ m/s

2. **Velocity Profile**: 
   - Linear (n=1): Parabolic vertical profile
   - Non-linear (n=3): "Plug flow" with deformation concentrated at base

3. **Physical Realism**: Non-linear model with Glen's law viscosity matches observed glacier behavior

See `plots/profile_comparison.png` for detailed visualization.

---

## 📖 Mathematical Proofs & Theory

The project includes rigorous mathematical analysis (see `solutions.pdf`):

### Galerkin Orthogonality
Proves that FEM error is orthogonal to the approximation space: `a(u - uₕ, vₕ) = 0`

### Best Approximation Property
Shows FEM solution minimizes energy norm error: `||u - uₕ||_E ≤ ||u - vₕ||_E` for all vₕ

### A Priori Error Estimate
Demonstrates mesh convergence: `||u - uₕ||_E ≤ Ch||u||_{H²}`

---

## � Research Paper

This repository includes a **complete research paper** (`solutions.pdf`) that presents the work in professional academic format:

- **10 pages** of rigorous analysis and results
- **Abstract** summarizing objectives, methods, and key findings  
- **3 formal theorems** with complete mathematical proofs
- **Professional figures** and convergence analysis tables
- **Bibliography** with proper academic citations
- **Structured sections**: Introduction, Mathematical Analysis, Numerical Methods, Results, Discussion, Conclusions

Read [`PAPER_TRANSFORMATION.md`](PAPER_TRANSFORMATION.md) for details on the paper structure.

---

## �🔍 Code Highlights

### Weak Form Implementation (UFL)

```python
# Define function spaces
V = FunctionSpace(mesh, ("Lagrange", 1))
u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = (2 * eta * ufl.Dx(u, 0) * ufl.Dx(v, 0) + 
     0.5 * eta * ufl.Dx(u, 1) * ufl.Dx(v, 1)) * ufl.dx

# Linear form
L = -rho * g * ufl.Dx(h, 0) * v * ufl.dx
```

### Picard Iteration Core

```python
def picard_iteration(mesh, epsilon=1e-10, max_iter=50, tol=1e-6):
    u_k = Function(V)  # Initial guess
    for k in range(max_iter):
        eta_k = compute_viscosity(u_k, epsilon)
        u_next = solve_linear_problem(eta_k)
        
        error = norm(u_next.vector - u_k.vector)
        if error < tol:
            return u_next, k+1
        
        u_k.assign(u_next)
```

---

## 📚 References & Attribution

This project is based on theoretical foundations from:

1. **Greve, R. & Blatter, H.** (2009). *Dynamics of Ice Sheets and Glaciers*. Springer.
2. **Jouvet, G. & Rappaz, J.** (2011). "A cut finite element method for non-Newtonian free surface flows in 2D—application to glacier modelling." *Computers & Fluids*.
3. **Glen, J.W.** (1955). "The creep of polycrystalline ice." *Proceedings of the Royal Society A*.

**Origin**: This work builds upon a group project from a computational glaciology course at **Stockholm University**.

**Individual contributions by Shek Lun Leung**:
- Complete FEniCSx reimplementation and code architecture
- Mathematical derivations and formal proofs (solutions.tex)
- Non-linear Picard solver development
- Comprehensive convergence analysis and visualization pipeline
- Professional documentation and repository structure
- Transformation into standalone research project

---

## 📈 Future Extensions

Potential areas for further development:

- [ ] **3D Full Stokes**: Extend to complete 3D ice flow model
- [ ] **Time Evolution**: Add transient simulation capabilities
- [ ] **Mesh Adaptation**: Implement adaptive mesh refinement
- [ ] **Parallel Computing**: Leverage MPI for large-scale simulations
- [ ] **Uncertainty Quantification**: Stochastic parameter studies
- [ ] **Inverse Problems**: Parameter estimation from field data

---

## 📧 Contact

**Shek Lun Leung**  
Email: sheklunleung.qai@proton.me  
GitHub: [github.com/alanspace](https://github.com/alanspace)  
LinkedIn: [linkedin.com/in/shek-lun-leung-alan](https://www.linkedin.com/in/shek-lun-leung-alan/)

---

## 📄 License

This project is available for academic and educational purposes. Please cite appropriately if used in research.

---

## 🙏 Acknowledgments

- **FEniCS Project** for the excellent finite element framework
- **Glacier data** provided by [data source]
- Original course instructors and teaching assistants

---

*Last Updated: January 2026*