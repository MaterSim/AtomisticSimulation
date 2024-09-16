"""
Here is a script to run NVE MD for 864 argon atoms.
You can run the code with 500 steps to monitor the evolution of total energies.
The current code does not generate the conserved total in each step.
You job is to debug the code and let it work in function.
"""

import numpy as np
from numba import njit     # Hint: used for speeping up the run
from time import time

def initialize_position(L):
    """
    Initialize positions in a 6*6*6 of the FCC unit cell

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
    sigma_v = np.sqrt(TEMPERATURE * KB / MASS)
    V = np.random.normal(0, sigma_v, (N, 3))
    # Center the velocities
    V -= np.mean(V, axis=0)
    return V

@njit
def LJ(r):
    """
    Lennard-Jones potential function
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
    N = len(R)
    F = np.zeros_like(R)
    PE = 0.0

    for i in range(N - 1):
        for j in range(i + 1, N):
            # Compute R between (i, j)
            rvec = R[i] - R[j]
            rvec -= np.round(rvec / L) * L
            rmag = np.linalg.norm(rvec)
            # Compute the potential Energy
            PE += LJ(rmag)
            # Compute and update forces
            force_mag = 24 * EPSILON * (2 * (SIGMA / rmag) ** 12 - (SIGMA / rmag) ** 6) / rmag
            forec_vec = force_mag * rvec
            F[i] -= forec_vec
            F[j] += forec_vec

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
    R += V * TIMESTEP + 0.5 * F / MASS * TIMESTEP ** 2
    R = np.mod(R, L)
    # Compute F_new
    PE, F_new = LJ_energy_forces(R, L)
    # Update V
    V += 0.5 * (F + F_new) / MASS * TIMESTEP

    # Update F
    F = F_new

    return R, V, F, PE


if __name__ == "__main__":
    np.random.seed(0)

    # System Parameters
    KB = 1.380649e-23     # Boltzmann constant in J/K
    TEMPERATURE = 94.4    # in K
    EPSILON = 120*KB      # in J (epsilon = 120 * k_B)
    SIGMA = 3.4e-10       # in meters (3.4 Å)
    L = 10.229 * SIGMA    # cubic box side length
    N = 864               # Number of atoms
    MASS = 39.948 * 1.66053906660e-27    # mass of Argon atom in kg

    # MD Parameters
    TIMESTEP = 1.0 * 1e-14  # time step in seconds
    num_steps = 200         # Number of steps to simulate

    # Initialize the system
    R = initialize_position(L)#; print(R[:2])
    V = initialize_velocity(N)#; print(V[:2])
    E, F = LJ_energy_forces(R, L)#; print(E, F)

    KEs = []
    PEs = []
    # MD propogation
    t0 = time()
    for step in range(500): #num_steps):
        R, V, F, PE = verlet_integration(R, V, F, L)
        KE = 0.5 * MASS * np.sum(V ** 2)
        KEs.append(KE)
        PEs.append(PE)
        if step % 10 == 0:
            print(f"Step {step:6d}, PE: {PE:.5e} KE: {KE:.5e} E_total: {PE+KE:.5e} CPU_time: {time()-t0:.2e}")
