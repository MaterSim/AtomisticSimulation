"""
Here is a script to run MD under different thermostat techniques.
You can run the code with 500-1000 steps to monitor the evolution of Kinetic energies.
So far, NVE, Anderson and Nose-Hoover methods work properly.
However, the Langevin method does not work well.
You job is to debug the code and let it work in function.
"""
import numpy as np
from numba import njit     # Hint: used for speeping up the run
import matplotlib.pyplot as plt

# Parameters
KB = 1.3806452 * 1e-23      # Boltzmann constant in J/K
TEMPERATURE = 94.4          # in K
EPSILON = 120.0 * KB        # in J (epsilon = 120 * k_B)
SIGMA = 3.4 * 1e-10         # in meters (3.4 Å)
MASS = 39.95*1.6747*1e-27   # mass of Argon atom in kg
L = 10.229 * SIGMA          # cubic box side length
N = 864                     # Number of atoms
TIMESTEP = 1.0 * 1e-14      # time step in seconds

def initialize_position(L):
    """
    Initialize positions in a 6*6*6 of the FCC unit cell
    864 atoms to a cubic.
    In a each FCC unit cell, there are 4 atoms. 4*6*6*6 = 864

    Args:
        L (float): unit length of the cubic box

    Returns:
        R (float) : (4*6*6*6, 3) array
    """
    # FCC unit cell fractional positions
    r = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ])
    a = L/6

    r *= a
    R = []
    for i in range(6):
        for j in range(6):
            for k in range(6):
                R.extend(r + np.array([i, j, k])*a)
    return np.array(R)


def initialize_velocity(N):
    """
    Initialize velocities using Maxwell-Boltzmann distribution

    Args:
        N (int): Number of atoms

    Returns:
        V (float) : (N, 3) array
    """
    # Standard deviation of the velocity distribution
    sigma = np.sqrt(TEMPERATURE * KB / MASS)  # kb*T = mv^2
    V = np.random.normal(0, sigma, (N, 3))   # V: (N, 3)

    # Center the velocities
    V -= np.mean(V, axis=0)
    return V

# Computation of Energy and Forces


@njit
def LJ(r):
    """
    Lennard-Jones potential function
    Two-body potential, for two particles
    """
    return 4 * EPSILON * ((SIGMA / r) ** 12 - (SIGMA / r) ** 6)


@njit
def LJ_energy_forces(R, L):
    """
    Compute the energy and forces from the given system

    Args:
        R (float) : (N, 3) array
        L (float): unit length of the cubic box

    Returns:
        PE (float): total energy
        F (float): atomic forces [N, 3] array
    """
    N = len(R)                # Number of atoms as a scalor
    F = np.zeros_like(R)      # forces [N, 3] as a 2D array
    PE = 0.0                  # total potential energy as a scalor

    for i in range(N-1):
        for j in range(i + 1, N):
            # Compute R between (i, j)
            r_vec = R[i] - R[j]
            r_vec -= np.round(r_vec / L) * L  # Peridoic boundary condition
            r_dist = np.linalg.norm(r_vec)    # sqrt(x**2 + y**2 + z**2)

            # Compute the potential Energy
            PE += LJ(r_dist)
            # Compute and update forces
            force_magnitude = EPSILON * \
                (48 * (SIGMA / r_dist) ** 12 - 24 *
                 (SIGMA / r_dist) ** 6) / r_dist**2  # dV/dr
            force_vec = r_vec * force_magnitude  # dV/dr * dr/d(x,y,z)
            F[i] += force_vec   # Make sure direction is consistent
            F[j] -= force_vec   # Make sure direction is consistent

    return PE, F


def verlet_integration(R, V, F, L):
    """
    Intergration based on the Verlet Velocity algorithm

    Args:
        R (float) : (N, 3) array
        V (float) : (N, 3) array
        F (float) : (N, 3) array
        L (float) : Cell length

    Returns:
        R (float) : (N, 3) array
        V (float) : (N, 3) array
        F (float) : (N, 3) array
    """
    # Update R
    R += V * TIMESTEP + 0.5 * F/MASS * TIMESTEP ** 2
    R = R % L

    # Compute F_new
    PE, F_new = LJ_energy_forces(R, L)

    # Update V
    V += 0.5 * (F + F_new) / MASS * TIMESTEP

    # Update F
    F = F_new

    return R, V, F


def langevin_thermostat(R, V, F, L, gamma):
    """
    Langevin thermostat
    """
    R += V * TIMESTEP + 0.5 * F/MASS * TIMESTEP ** 2
    R = R % L

    # Update velocities with deterministic part
    V += 0.5 * F * TIMESTEP / MASS

    # Apply friction and random force (stochastic part)
    V -= gamma * V * TIMESTEP
    sigma = np.sqrt(2 * gamma * KB * TEMPERATURE / MASS)
    V += np.random.randn(len(V), 3) * np.sqrt(TIMESTEP) * sigma

    # Update forces
    PE, F_new = LJ_energy_forces(R, L)
    V += 0.5 * (F + F_new) / MASS * TIMESTEP
    F = F_new

    return R, V, F


def anderson_thermostat(V, nu=0.5):
    """
    Anderson thermostat
    """
    sigma = np.sqrt(KB * TEMPERATURE / MASS)
    # Randomly assign new velocities
    lists = []
    for i in range(len(V)):
        if np.random.rand() < nu:
            V[i] = np.random.normal(0, sigma, 3)
            lists.append(i)
    # print(f"anderson_thermostat, reassign V for {len(lists)} particles")

    return V


def Nose_Hoover_thermostat(R, V, F, xi, L, Q):
    """
    Nose-Hoover thermostat
    """
    R += V * TIMESTEP + 0.5 * F/MASS * TIMESTEP ** 2
    R = R % L
    PE, F_new = LJ_energy_forces(R, L)
    V += 0.5 * (F + F_new) / MASS * TIMESTEP
    V *= (1 - 0.5 * xi * TIMESTEP) / (1 + 0.5 * xi * TIMESTEP)
    KE = 0.5 * np.sum(V**2) * MASS
    xi += TIMESTEP * (2 * KE / (3 * len(R) * KB * TEMPERATURE) - 1) / (Q*MASS)
    return R, V, F, xi


def MD(thermostat=None, nu=0.1, gamma=1e-13, Q=1.0, num_steps=500):
    """
    Run MD simulation

    Args:
        thermostat (str): "Langevin", "Anderson", "Nose-Hoover" or None
        nu (float): Anderson thermostat parameter
        gamma (float): Langevin thermostat parameter
        Q (float): Nose-Hoover thermostat parameter
        num_steps (int): Number of steps to simulate
    """
    # Initialize system
    R = initialize_position(L)
    V = initialize_velocity(N)
    E, F = LJ_energy_forces(R, L)
    KEs = []
    PEs = []
    TEs = []
    xi = 0.0  # Used by Nose-Hoover

    # MD propogation
    for step in range(num_steps):
        if thermostat == "Nose-Hoover":
            R, V, F, xi = Nose_Hoover_thermostat(R, V, F, xi, L, Q)
        elif thermostat == "Langevin":
            R, V, F = langevin_thermostat(R, V, F, L, gamma)
        else:
            R, V, F = verlet_integration(R, V, F, L)
            if thermostat == "Anderson":
                V = anderson_thermostat(V, nu)

        # Compute PE, KE, and TE
        PE, _ = LJ_energy_forces(R, L)
        KE = 0.5 * np.sum(V**2) * MASS
        KEs.append(KE)
        PEs.append(PE)
        TEs.append(PE+KE)
        if step % 10 == 0:
            s = np.exp(-xi * TIMESTEP)
            print(f"Step {step:6d}, PE: {PE:.5e} KE: {KE:.5e} E: {PE+KE:.5e} s: {s:.4f}")
    return KEs, PEs, TEs


if __name__ == "__main__":
    results = []
    for thermostat in ["Langevin"]:
        print(f"Simulation with {thermostat} thermostat")
        KEs, PEs, TEs = MD(thermostat=thermostat)
        results.append((thermostat, KEs, PEs, TEs))

    for result in results:
        thermostat, KEs, PEs, TEs = result
        if thermostat is None: thermostat = "NVE"
        plt.plot(KEs, label=thermostat)
    plt.legend()
    plt.savefig('KE-comparison.png')
