#!/usr/bin/env python3

# --- 1. Import necessary libraries ---
from mpi4py import MPI
from dolfinx import mesh, fem, plot
import numpy
import ufl
from dolfinx import default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import pyvista

# --- 2. Create the Mesh and Function Space ---
domain = mesh.create_unit_square(MPI.COMM_WORLD, 16, 16)

# --- PARAMETER TO CHANGE FOR EXPERIMENTS ---
# Change the '1' to a '2' for second-order elements (P2)
element_order = 1 
# ---

V = fem.functionspace(domain, ("Lagrange", element_order))

# --- 3. Define the Boundary Conditions and RHS ---

# --- PARAMETER TO CHANGE FOR EXPERIMENTS ---
# Change this to "sine" for the problem in Ex. 11
problem_type = "quadratic" 
# ---

if problem_type == "quadratic":
    # --- Problem Data for Ex. 10 (Quadratic Solution) ---
    uD = fem.Function(V)
    uD.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)
    f = fem.Constant(domain, default_scalar_type(-6))
    # For error calculation, the exact solution needs a space rich enough to hold it.
    # P2 elements can hold a quadratic function, so V is sufficient if order is 2.
    V_ex = fem.functionspace(domain, ("Lagrange", 2))
    uex = fem.Function(V_ex)
    uex.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)
    plot_filename = f"quadratic_solution_p{element_order}.png"

elif problem_type == "sine":
    # --- Problem Data for Ex. 11 (Sine Solution) ---
    uD = fem.Function(V)
    uD.interpolate(lambda x: 0.0 * x[0])
    x_coords = ufl.SpatialCoordinate(domain)
    f = 2 * numpy.pi**2 * ufl.sin(numpy.pi * x_coords[0]) * ufl.sin(numpy.pi * x_coords[1])
    # For error calculation, use a higher-order space for the exact solution.
    V_ex = fem.functionspace(domain, ("Lagrange", 3))
    uex = fem.Function(V_ex)
    uex.interpolate(lambda x: numpy.sin(numpy.pi * x[0]) * numpy.sin(numpy.pi * x[1]))
    plot_filename = f"sine_solution_p{element_order}.png"

# --- Apply the Boundary Conditions ---
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs)

# --- 4. Define and Solve the Variational Problem ---
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    petsc_options_prefix="poisson_solver" # <-- THE FIX IS HERE
)
uh = problem.solve()

# --- 5. Compute the Error ---
L2_error = fem.form(ufl.inner(uh - uex, uh - uex) * ufl.dx)
error_local = fem.assemble_scalar(L2_error)
error_L2 = numpy.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

if domain.comm.rank == 0:
    print(f"\n--- Results for {problem_type} problem with P{element_order} elements ---")
    print(f"L2 Error: {error_L2:.3e}")

# --- 6. Plotting the Solution ---
if domain.comm.rank == 0:
    print(f"Generating plot: {plot_filename}")

u_grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(V))
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")
warped = u_grid.warp_by_scalar()

# Save the 3D surface plot
plotter = pyvista.Plotter(off_screen=True)
plotter.add_mesh(warped, show_edges=True, show_scalar_bar=True)
plotter.view_isometric()
plotter.screenshot(plot_filename)