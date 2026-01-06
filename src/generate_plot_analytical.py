
import numpy as np
import matplotlib.pyplot as plt

def generate_analytical_plot():
    # Parameters
    rho = 910.0
    g = 9.81
    eta = 1.0
    L = 5000.0
    H = 200.0
    alpha = 0.03 # slope roughly
    
    # Analytical solution for parallel flow assumption
    # 2 eta dx(dx u) term is 0 if u only depends on z.
    # 0.5 eta dz(dz u) = rho g dx h
    # dz(dz u) = 2 rho g (-alpha) / eta
    # u(z) = (rho g alpha / eta) * (2*H*z - z**2)
    
    x = np.linspace(0, L, 100)
    z = np.linspace(0, H, 50)
    X, Z = np.meshgrid(x, z)
    
    # Velocity depends mostly on Z in this simple approximation
    C = (rho * g * alpha) / eta
    # Scaling C to get reasonable numbers for visualization if needed, 
    # but let's keep physical or normalized. 
    # With eta=1, numbers will be huge. 
    # The prompt said "set viscosity to any scalar number e.g. 1".
    
    U = C * (2*H*Z - Z**2)
    
    plt.figure(figsize=(10, 5))
    contour = plt.contourf(X, Z, U, levels=20, cmap='viridis')
    plt.colorbar(contour, label='Velocity u (m/s)')
    plt.title('Analytical Velocity Field (Approximation for Q4)')
    plt.xlabel('Distance along flowline (m)')
    plt.ylabel('Height above bedrock (m)')
    plt.tight_layout()
    plt.savefig('plots/velocity_field.png')
    print("Analytical plot saved to plots/velocity_field.png")

if __name__ == "__main__":
    generate_analytical_plot()
