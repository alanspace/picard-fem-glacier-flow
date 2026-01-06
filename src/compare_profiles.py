
import sys
import os
import numpy as np
import ufl
import basix.ufl
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
from petsc4py import PETSc
import matplotlib.pyplot as plt

# Import helpfunctions
sys.path.append(os.path.abspath("src"))
try:
    import helpfunctions
except ImportError:
    pass

def get_nonlinear_viscosity(u, eps, A_sec):
    ux = ufl.Dx(u, 0)
    uz = ufl.Dx(u, 1)
    strain_term = ux**2 + uz**2 + eps
    eta = 0.25 * (A_sec**(-1/3)) * (strain_term**(-1.0/3.0))
    return eta

def compare_results():
    # Load Mesh
    mesh_path = "data/arolla.xdmf"
    if not os.path.exists(mesh_path):
        print("Mesh not found")
        return
    msh = helpfunctions.loadmesh(mesh_path)
    
    # Setup Spaces
    element = basix.ufl.element("Lagrange", "triangle", 1)
    try:
        V = fem.functionspace(msh, element)
    except AttributeError:
        V = fem.FunctionSpace(msh, element)
        
    facet_tags = helpfunctions.mark_bed_surface(msh, bed_id=2, surface_id=1)
    bed_facets = facet_tags.find(2)
    bed_dofs = fem.locate_dofs_topological(V, msh.topology.dim - 1, bed_facets)
    bcs = [fem.dirichletbc(PETSc.ScalarType(0), bed_dofs, V)]
    
    # Parameters
    rho = 910.0
    g = 9.81
    
    # LINEAR CASE
    # We used eta = 1e13 in main.py
    eta_lin = 1.0e13
    u_lin = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    h_fun = helpfunctions.get_h(V, facet_tags, surface_id=1)
    
    a_lin = (2 * eta_lin * ufl.Dx(u_lin, 0) * ufl.Dx(v, 0) + 
             0.5 * eta_lin * ufl.Dx(u_lin, 1) * ufl.Dx(v, 1)) * ufl.dx
    L_lin = - rho * g * ufl.Dx(h_fun, 0) * v * ufl.dx
    
    try:
        problem_lin = LinearProblem(a_lin, L_lin, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="lin_comp")
    except TypeError:
        problem_lin = LinearProblem(a_lin, L_lin, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh_lin = problem_lin.solve()
    max_lin = np.max(np.abs(uh_lin.x.array))
    print(f"Max Linear Velocity: {max_lin:.2e} m/s")

    # NON-LINEAR CASE
    # Parameters from Q8
    sec_per_yr = 31557600.0
    A_yr = 1.0e-16
    A_sec = A_yr / sec_per_yr
    eps = 1.0e-10
    
    # Solve (simplified Picard, just enough to get result)
    uh_non = fem.Function(V)
    uh_non.x.array[:] = 0.0 # Start from 0
    
    # We need to loop
    for i in range(15):
        u_trial = ufl.TrialFunction(V)
        eta_k = get_nonlinear_viscosity(uh_non, eps, A_sec)
        a_non = (2 * eta_k * ufl.Dx(u_trial, 0) * ufl.Dx(v, 0) + 
                 0.5 * eta_k * ufl.Dx(u_trial, 1) * ufl.Dx(v, 1)) * ufl.dx
        L_non = - rho * g * ufl.Dx(h_fun, 0) * v * ufl.dx
        
        try:
             problem = LinearProblem(a_non, L_non, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix=f"nonlin_comp_{i}")
             u_new = problem.solve()
        except TypeError:
             problem = LinearProblem(a_non, L_non, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
             u_new = problem.solve()
        except Exception as e:
             print(f"Error {e}")
             break
        
        diff = np.linalg.norm(u_new.x.array - uh_non.x.array) / np.linalg.norm(u_new.x.array)
        uh_non.x.array[:] = u_new.x.array[:]
        if diff < 1e-3:
            break
            
    max_non = np.max(np.abs(uh_non.x.array))
    print(f"Max Non-Linear Velocity: {max_non:.2e} m/s")
    
    # COMPARE at a vertical cross section
    # Find a point in the middle x = 2500 roughly
    x_coords = V.tabulate_dof_coordinates()[:, 0]
    z_coords = V.tabulate_dof_coordinates()[:, 1]
    
    # Filter points near x=2500
    mid_mask = (x_coords > 2400) & (x_coords < 2600)
    
    z_mid = z_coords[mid_mask]
    u_lin_mid = uh_lin.x.array[mid_mask]
    u_non_mid = uh_non.x.array[mid_mask]
    
    # Sort by z
    idx_lin = np.argsort(z_mid)
    
    plt.figure(figsize=(6, 8))
    plt.plot(u_lin_mid[idx_lin], z_mid[idx_lin], 'b.', label='Linear (eta=const)')
    plt.plot(u_non_mid[idx_lin], z_mid[idx_lin], 'r.', label='Non-Linear (eps=1e-10)')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Height z (m)')
    plt.title('Vertical Velocity Profile Comparison')
    plt.legend()
    plt.grid(True)
    plt.savefig('plots/profile_comparison.png')
    print("Saved plots/profile_comparison.png")

if __name__ == "__main__":
    compare_results()
