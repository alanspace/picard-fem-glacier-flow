"""
Linear Glacier Flow Solver using FEniCSx
=========================================

This module implements the Linear First-Order Stokes approximation for glacier flow
using the Finite Element Method (FEM) via FEniCSx.

Physical Model:
--------------
Solves the steady-state linear problem:
    2∂ₓ(η∂ₓu) + ½∂ᵤ(η∂ᵤu) = ρg∂ₓh

where:
    η : constant effective viscosity (Pa·s)
    u : horizontal velocity component (m/s)
    ρ : ice density (910 kg/m³)
    g : gravitational acceleration (9.81 m/s²)
    h : surface elevation (m)

Boundary Conditions:
-------------------
- Bedrock (Γ_b): No-slip condition, u = 0 (Dirichlet)
- Surface (Γ_s): Stress-free condition (Natural BC)

Implementation:
--------------
- Element Type: Continuous Lagrange (P1)
- Solver: Direct LU decomposition
- Mesh: Arolla glacier geometry (data/arolla.xdmf)

Usage:
------
    $ python main.py

Output:
-------
- Velocity field visualization (plots/velocity_field.png)
- VTX file for ParaView (data/velocity_linear.bp)

Author: [Your Name]
Date: January 2026
"""

import sys
import os
import numpy as np
import ufl
from dolfinx import fem, io
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc
import matplotlib.pyplot as plt

# Ensure src is in path to import helpfunctions locally if running from root
sys.path.append(os.path.abspath("src"))
try:
    import helpfunctions
except ImportError:
    # Fallback if running from src directory
    sys.path.append(os.path.abspath("."))
    import helpfunctions

def solve_linear_glacier():
    # 1. Load Mesh
    try:
        # User mesh is in data/arolla.xdmf
        mesh_path = "data/arolla.xdmf"
        if not os.path.exists(mesh_path):
             # Fallback if I got path wrong
             mesh_path = "arolla.xdmf"
        
        msh = helpfunctions.loadmesh(mesh_path)
        print(f"Mesh loaded: {msh.name}")
    except Exception as e:
        print(f"Error loading mesh: {e}")
        return

    # 2. Parameters (From Ahlkrona et al. 2017 / Project Description)
    # Using SI units: Pa, m, s
    # A = 1e-16 Pa^-3 yr^-1
    years_to_seconds = 31557600.0
    A_yr = 1.0e-16 
    n = 3.0
    # A_sec = A_yr / years_to_seconds # 3.16e-24
    
    # Note: Question 4 asks for a LINEAR version with constant viscosity.
    # It suggests "eta = 1" or similar.
    # If we use physical rho and g, we get huge forcing.
    # Let's use eta = 1.0e14 roughly (ice viscosity) if we want physical numbers,
    # OR follow the prompt "eta to any scalar number e.g. 1".
    # BUT if we use eta=1 and rho*g ~ 9000, velocity will be ~ 10^4 m/s.
    # To get reasonable plots, we should either use physical eta or scaled parameters.
    # Q: "better to use the values in Table 1 ... scaled to more favorable units"
    # Usually this means: MPa, year.
    # rho = 910 kg/m^3 = 910 * 1e-6 / 1e-18 ?? No.
    # Let's stick to SI but pick a physical constant viscosity for the plot to look nice.
    # A typical viscosity is 1/2 A^(-1/n) eps^( (1-n)/n ).
    
    rho = 910.0
    g = 9.81
    
    # Let's use the explicit instruction: "set the viscosity eta to any scalar number for now e.g. 1"
    # But for the plot to be meaningful in the report, meaningful units help. 
    # However, I will strictly follow "e.g. 1" to answer the question's requirement for linearity 
    # and just note the scale in the plot is arbitrary. 
    # WAIT: "compare it to the plots of Figure 5...". That is for Q8 (non-linear).
    # Q4 just says "Include a plot of the velocity field".
    # I will use eta = 1e7 to keep numbers somewhat tame, or 1.0 if strictly following "e.g. 1".
    # I'll use eta = 1e13 to simulate real ice.
    eta = 1.0e13 

    # 3. Function Space
    # Explicit element definition for robustness (dolfinx 0.8+)
    from dolfinx import default_scalar_type
    import basix.ufl
    # P1 element using basix.ufl
    # Assumes 2D triangle mesh as per arolla.xdmf
    element = basix.ufl.element("Lagrange", "triangle", 1)
    # Check if we need to pass shape or similar? usually "triangle" is fine.
    # Note: for dolfinx 0.8+, FunctionSpace might expect an element object from basix.ufl
    # In this version (0.10.0), it seems we must use the factory function 'functionspace' (lowercase).
    try:
        V = fem.functionspace(msh, element)
    except AttributeError:
        # Fallback if functionspace is not found but FunctionSpace expects something else
        # Check if FunctionSpace works with keyword args or different signature
        V = fem.FunctionSpace(msh, element)

    # 4. Boundaries
    # Use helpfunctions to mark boundaries
    # bed_id = 2, surface_id = 1
    facet_tags = helpfunctions.mark_bed_surface(msh, bed_id=2, surface_id=1)
    
    # Dirichlet BC on Bed (id=2) -> u=0
    fdim = msh.topology.dim - 1
    bed_facets = facet_tags.find(2)
    bed_dofs = fem.locate_dofs_topological(V, fdim, bed_facets)
    bcs = [fem.dirichletbc(PETSc.ScalarType(0), bed_dofs, V)]
    
    # 5. Surface Slope (RHS)
    # We need dh/dx.
    # The mesh itself provides geometry. "surface" geometry is x, h(x).
    # We can compute dh/dx from the mesh using helpfunctions? No, helpfunctions has get_h.
    # Ideally, we calculate gradients of the coordinate field or h field.
    
    # Compute h(x)
    # h_fun is a function on the surface extruded down?
    h_fun = helpfunctions.get_h(V, facet_tags, surface_id=1)
    
    # Compute dh/dx
    # We can take derivative of h_fun.
    # Note: h_fun is in V (CG1). Its geometric derivative is discontinuous (DG0).
    # For the weak form `rho * g * dx(h) * v * dx`, we can use `Dx(h_fun, 0)`.
    
    # 6. Variational Problem (Linear)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    # Equation 15: 2 dx(eta dx u) + 0.5 dz(eta dz u) = rho g dx h
    # Weak form:
    # - int (2 eta dx u dx v + 0.5 eta dz u dz v) = int rho g dx h v
    # (assuming natural BCs vanish or stress free matches this)
    # Note signs: Partial integration of LHS moves derivative to v with minus sign.
    # LHS term: int (terms) v.
    # -> - int (2 eta u_x v_x + 0.5 eta u_z v_z) ...
    # So:
    # a(u, v) = int (2 eta u_x v_x + 0.5 eta u_z v_z)
    # L(v) = - int (rho g h_x v)
    # Use ufl.Dx for spatial derivatives
    
    a = (2 * eta * ufl.Dx(u, 0) * ufl.Dx(v, 0) + 
         0.5 * eta * ufl.Dx(u, 1) * ufl.Dx(v, 1)) * ufl.dx
    
    L = - rho * g * ufl.Dx(h_fun, 0) * v * ufl.dx
    
    # Solve
    try:
        problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, 
                                petsc_options_prefix="glacier_linear")
    except TypeError:
        # Fallback if prefix is not required in other versions, but error says it is.
        problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    
    uh = problem.solve()
    uh.name = "Velocity (Linear)"

    # 7. Plotting
    try:
        # Save solution to file
        with io.VTXWriter(msh.comm, "data/velocity_linear.bp", [uh], engine="BP4") as vtx:
            vtx.write(0.0)
            
        # Matplotlib plot
        # Extract values on nodes
        # Triangulation
        x = msh.geometry.x[:, 0]
        z = msh.geometry.x[:, 1]
        
        # We need to map degrees of freedom to vertex indices if they differ, 
        # but for CG1 they are usually aligned or we use fem.Function.eval.
        # Simplest reference: uh.x.array corresponds to dof indices.
        # We need a map from dof to vertex.
        
        # For simplicity in this script, let's use a scatter or tricontour if we can align arrays.
        # Dolfinx V.tabulate_dof_coordinates() gives coordinates of dofs.
        dof_coords = V.tabulate_dof_coordinates()
        x_dofs = dof_coords[:, 0]
        z_dofs = dof_coords[:, 1]
        u_vals = uh.x.array.real
        
        plt.figure(figsize=(12, 6))
        # Use tricontourf
        plt.tricontourf(x_dofs, z_dofs, u_vals, levels=50, cmap='viridis')
        plt.colorbar(label='Horizontal Velocity u (m/s)')
        plt.title(f'Velocity Field (Linear, eta={eta:.1e})')
        plt.xlabel('x (m)')
        plt.ylabel('z (m)')
        plt.axis('equal')
        output_plot = 'plots/velocity_field.png'
        plt.savefig(output_plot)
        print(f"Plot saved to {output_plot}")
        
    except Exception as e:
        print(f"Plotting failed: {e}")

if __name__ == "__main__":
    solve_linear_glacier()
