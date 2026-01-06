import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# --- Define the Mesh (from the exercise) ---

# Node coordinates: P matrix from the previous question
# The columns correspond to nodes N1 to N9
P = np.array([
    [0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0],  # x-coordinates
    [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]   # y-coordinates
])

# Connectivity matrix: T matrix (nodes are 0-indexed in Python)
# Each column is a triangle [N_i, N_j, N_k]
T = np.array([
    [0, 0, 1, 1, 3, 3, 4, 4],  # 1st vertex index of each triangle
    [1, 4, 2, 5, 4, 7, 5, 8],  # 2nd vertex index
    [4, 3, 5, 4, 7, 6, 8, 7]   # 3rd vertex index
]).T # Transpose to get 8x3 shape

# Create a matplotlib Triangulation object to work with
mesh_triangulation = Triangulation(P[0, :], P[1, :], T)

# --- Function to Define the Hat Basis Functions ---

def get_phi_values(triangulation, home_node_index):
    """
    Calculates the values of the hat function phi_i for a given mesh.
    
    The value of the hat function phi_i is 1 at its home node i, and 0 at all other nodes.
    It is piecewise linear over the mesh.
    """
    z = np.zeros(triangulation.x.shape)
    z[home_node_index] = 1.0
    return z

# --- Generate the plots ---

# Create a figure with two subplots side-by-side
fig, (ax1, ax2) = plt.subplots(
    1, 2, 
    figsize=(14, 6), 
    subplot_kw={'projection': '3d'}
)

# --- Plot for Hat Function phi_1 ---

# Get the values of phi_1 at each node
phi1_values = get_phi_values(mesh_triangulation, home_node_index=0) # N1 is at index 0

# Create the 3D surface plot
ax1.plot_trisurf(mesh_triangulation, phi1_values, cmap='viridis', edgecolor='black', linewidth=0.2)
ax1.set_title('Hat Function $\phi_1$ (for Node $N_1$)', fontsize=16)
ax1.set_xlabel('$x_1$')
ax1.set_ylabel('$x_2$')
ax1.set_zlim(0, 1.2) # Set z-limit for consistent view

# --- Plot for Hat Function phi_5 ---

# Get the values of phi_5 at each node
phi5_values = get_phi_values(mesh_triangulation, home_node_index=4) # N5 is at index 4

# Create the 3D surface plot
ax2.plot_trisurf(mesh_triangulation, phi5_values, cmap='cividis', edgecolor='black', linewidth=0.2)
ax2.set_title('Hat Function $\phi_5$ (for Node $N_5$)', fontsize=16)
ax2.set_xlabel('$x_1$')
ax2.set_ylabel('$x_2$')
ax2.set_zlim(0, 1.2) # Set z-limit for consistent view

# Display the final plot
plt.tight_layout()
plt.show()