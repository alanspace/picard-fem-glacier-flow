#!/usr/bin/env python3

# --- 1. Import necessary libraries ---
from mpi4py import MPI
from dolfinx import mesh, fem, plot
import numpy
import ufl
from dolfinx import default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import pyvista
# Import the custom helper functions
import helpfunctions as hf

# --- 2. Load the Mesh and Define Function Space ---

# Load the Arolla glacier mesh using the helper function
domain = hf.loadmesh("arolla.xdmf")

# Mark the boundaries of the glacier mesh using the helper function
# bed_id=2 (bottom) and surface_id=1 (top) are the default tags
bed_id, surface_id = 2, 1
facet_tags = hf.mark_bed_surface(domain, bed_id=bed_id, surface_id=surface_id)

# Define the finite element function space V on this new mesh.
# Using P1 (first-order linear) elements for this example.
V = fem.functionspace(domain, ("Lagrange", 1))

# --- 3. Define the Boundary Conditions ---

# Define a constant value of 0 for the Dirichlet BC on the glacier bed.
# This represents a "no-slip" or frozen bed condition.
u_D = fem.Constant(domain, default_scalar_type(0.0))

# Find the facets (edges) that correspond to the bed_id tag.
fdim = domain.topology.dim - 1
bed_facets = facet_tags.find(bed_id)

# Locate the degrees of freedom on these specific bed facets.
bed_dofs = fem.locate_dofs_topological(V, fdim, bed_facets)

# Create the Dirichlet boundary condition object.
bc = fem.dirichletbc(u_D, bed_dofs, V)

# --- 4. Define the Variational Problem ---

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Define the physical constants for the problem.
# These values should be checked against the project description.
rho = 910.0  # Ice density (kg/m^3)
g = 9.81     # Gravity (m/s^2)
eta = 1e13   # A constant viscosity (Pa s), typical value for ice.

# The right-hand side (f) is the gravitational driving stress.
# It is proportional to the surface slope, d(s)/dx.
s = hf.get_h(V, facet_tags, surface_id=surface_id) # Get surface elevation 's'
f = rho * g * ufl.Dx(s, 0) # f = rho * g * (surface slope)

# Define the bilinear form a(u,v) for the linear First-Order (FO) model.
# This corresponds to the weak form of the term -div(eta * grad(u))
# For the FO model, this is integral( 4*eta*du/dx*dv/dx + eta*du/dy*dv/dy )
# Note: In the project, the vertical coordinate is often 'z' instead of 'y'.
# ufl.Dx(u, 0) is the derivative w.r.t. x.
# ufl.Dx(u, 1) is the derivative w.r.t. y (or z).
a = (4 * eta * ufl.Dx(u, 0) * ufl.Dx(v, 0) + eta * ufl.Dx(u, 1) * ufl.Dx(v, 1)) * ufl.dx
L = f * v * ufl.dx

# --- 5. Assemble and Solve the Linear System ---
problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    petsc_options_prefix="Glacier"
)
uh = problem.solve()
uh.name = "u" # Give the solution a name

# --- 6. Plotting and Saving the Solution ---

if MPI.COMM_WORLD.rank == 0:
    print("\n--- Generating plot for glacier geometry ---")

# Create the geometric data for plotting from the FEniCS function space.
u_grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(V))

# Add the solution data to the grid.
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")

# Create a "warped" version of the grid for 3D surface plotting.
warped = u_grid.warp_by_scalar()

# --- Save the 3D plot to a .png file for the report ---
plotter_3D_save = pyvista.Plotter(off_screen=True)
plotter_3D_save.add_mesh(warped, show_edges=True, show_scalar_bar=True)

# Enable lighting that follows the camera to ensure the object is lit.
plotter_3D_save.enable_lightkit()

# Set a very generous clipping range to avoid the object being cut out by the camera.
plotter_3D_save.camera.clipping_range = (1e-2, 1e9)

# Automatically move the camera to frame the object.
plotter_3D_save.reset_camera() 

# Set a nice, standard 3D viewpoint after resetting.
plotter_3D_save.view_isometric() 

# Save the screenshot.
filename_3D = "glacier_solution.png"
plotter_3D_save.screenshot(filename_3D)

if MPI.COMM_WORLD.rank == 0:
    print(f"Saved 3D plot for report to '{filename_3D}'")