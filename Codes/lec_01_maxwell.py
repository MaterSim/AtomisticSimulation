import numpy as np
import matplotlib.pyplot as plt

def generate_velocities(num_particles, temperature, mass, k_B):
    # Boltzmann constant in appropriate units (e.g., J/K or eV/K)

    # Standard deviation of the velocity distribution
    sigma_v = np.sqrt(k_B * temperature / mass)

    # Generate velocities from a normal distribution for each component
    velocities = np.random.normal(0, sigma_v, (num_particles, 3))

    # Calculate the speed for each particle
    speeds = np.linalg.norm(velocities, axis=1)

    return velocities, speeds

# Example usage
num_particles = 10000       # Number of particles
temperature = 300          # Temperature in Kelvin
mass = 1.67e-27            # Mass of a particle in kg (e.g., proton)
k_B = 1.380649e-23  # J/K

velocities, speeds = generate_velocities(num_particles, temperature, mass, k_B)

# Plot the histogram of speeds
plt.hist(speeds, bins=50, density=True, alpha=0.6, color='g')

# Overlay the Maxwell-Boltzmann theoretical distribution
v = np.linspace(0, np.max(speeds), 200)
distribution = (4 * np.pi * v**2) * (mass / (2 * np.pi * k_B * temperature))**(3/2) * np.exp(-mass * v**2 / (2 * k_B * temperature))
plt.plot(v, distribution, linewidth=2, color='r')
plt.xlabel('Speed (m/s)')
plt.ylabel('Probability Density')
plt.title('Maxwell-Boltzmann Speed Distribution')
plt.tight_layout()
plt.savefig('lec_01_MaxWell.pdf')
