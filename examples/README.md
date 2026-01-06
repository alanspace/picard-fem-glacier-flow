# Examples Directory

This directory contains tutorial examples demonstrating fundamental FEM concepts used in the glacier flow simulations.

---

## Overview

These examples serve as:
1. **Learning Resources**: Step-by-step introduction to FEniCSx
2. **Testing Tools**: Verify your installation works correctly
3. **Reference Implementations**: Simple templates for similar problems

---

## Files

### `poisson2D.py`
**Basic 2D Poisson Equation Solver**

Solves the canonical PDE:
```
-∇²u = f  in Ω
u = 0    on ∂Ω
```

**Concepts Covered**:
- Mesh generation with built-in geometries
- Defining function spaces
- Weak formulation
- Dirichlet boundary conditions
- Basic visualization

**Usage**:
```bash
python examples/poisson2D.py
```

**Expected Output**: 
- Console: Solution computed successfully
- Plot: Heat distribution on unit square

---

### `poisson_project.py`
**Extended Poisson Problem**

Similar to `poisson2D.py` but with:
- Custom source term
- Different boundary conditions
- More detailed plotting
- VTU file export for ParaView

**Usage**:
```bash
python examples/poisson_project.py
```

---

### `2D_Hat.py`
**Hat Function Demonstration**

Illustrates:
- Basis function construction
- Finite element shape functions
- Nodal interpolation

**Concepts**:
- Piecewise linear (P1) elements
- Lagrange basis functions
- Element-wise assembly

**Usage**:
```bash
python examples/2D_Hat.py
```

---

## Learning Path

**Recommended Order**:

1. **Start Here**: `poisson2D.py`
   - Understand basic FEM workflow
   - Get familiar with FEniCSx syntax

2. **Extend Concepts**: `poisson_project.py`
   - See how to customize problems
   - Learn about different output formats

3. **Theory Deep-Dive**: `2D_Hat.py`
   - Explore underlying approximation theory
   - Understand shape functions

4. **Apply to Glaciers**: `src/main.py` and `src/nonlinear_solver.py`
   - Use concepts in realistic simulation
   - Handle non-linear problems

---

## Modifications

These examples can be modified to explore:

### Different PDEs
```python
# Heat equation: ∂u/∂t - ∇²u = f
# Advection-diffusion: ∂u/∂t + v·∇u - ∇²u = f
# Helmholtz: -∇²u + k²u = f
```

### Different Domains
```python
# Circle, L-shape, custom geometries
from dolfinx.mesh import create_circle, create_rectangle
mesh = create_circle(MPI.COMM_WORLD, radius=1.0, refinement=2)
```

### Different Elements
```python
# Quadratic elements (P2)
element = basix.ufl.element("Lagrange", "triangle", 2)

# Discontinuous Galerkin (DG)
element = basix.ufl.element("DG", "triangle", 1)
```

---

## Debugging Tips

If examples fail:

1. **Check Python Environment**:
   ```bash
   conda activate fenicsx-env
   python -c "import dolfinx; print(dolfinx.__version__)"
   ```

2. **Run from Project Root**:
   ```bash
   # Correct
   python examples/poisson2D.py
   
   # May fail due to imports
   cd examples && python poisson2D.py
   ```

3. **Verify Display**:
   - For headless servers, use `matplotlib.use('Agg')`
   - For SSH X-forwarding, enable X11 forwarding

---

## Further Reading

- **FEniCSx Tutorial**: https://jsdokken.com/dolfinx-tutorial/
- **Finite Element Method Theory**: *The Finite Element Method* by Hughes
- **Glacier Modeling**: See `solutions.pdf` for domain-specific background

---

*These examples are adapted from standard FEniCSx tutorials and simplified for educational purposes.*
