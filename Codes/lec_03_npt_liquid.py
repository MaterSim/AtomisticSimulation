"""
Here is a script to run MD under different barostat techniques.
You can run the code with 500 steps to monitor the evolution of Volumes.
You job is to debug the code and let it work as expected.
"""
import numpy as np
from numba import njit     # Hint: used for speeping up the run
import matplotlib.pyplot as plt

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


def Nose_Hoover_thermostat(R, V, F, xi, L, Q):
    """
    Nose-Hoover thermostat
    """
    # Update R
    R += V * TIMESTEP + 0.5 * F/MASS * TIMESTEP ** 2
    R = R % L

    # Update forces
    PE, F_new = LJ_energy_forces(R, L)
    V += 0.5 * (F + F_new) / MASS * TIMESTEP
    V *= (1 - 0.5 * xi * TIMESTEP) / (1 + 0.5 * xi * TIMESTEP)

    # Update xi
    KE = 0.5 * np.sum(V**2) * MASS
    xi += TIMESTEP * (2 * KE / (3 * len(R) * KB * TEMPERATURE) - 1) / (Q*MASS)

    # Update forces
    F = F_new
    return R, V, F, xi

def compute_pressure(R, V, volume):
    """
    Compute the scalor pressure using the virial equation.

    Args:
        R (np.array): N x 3 array of particle positions.
        V (np.array): N x 3 array of particle velocities.
        volume (float): volume of the simulation box

    Returns:
        float: Computed pressure of the system.
    """
    N = len(R)
    L = volume ** (1.0 / 3.0)
    # kinetic contribution to the pressure
    P_KE = MASS * np.sum(V**2)#; print("KE", 0.5 * np.sum(V**2) * MASS, "Volume", volume)
    P_KE /= (3 * volume)
    PE, _ = LJ_energy_forces(R, L)#; print("PE", PE)

    # virial contribution to the pressure
    P_virial = 0.0
    for i in range(N-1):
        for j in range(i+1, N):
            r_vec = R[i] - R[j]
            r_vec -= np.round(r_vec / L) * L  # Peridoic boundary condition
            r_dist = np.linalg.norm(r_vec)    # sqrt(x**2 + y**2 + z**2)

            f_magnitude = EPSILON * \
                (48 * (SIGMA / r_dist) ** 12 - 24 *
                 (SIGMA / r_dist) ** 6) / r_dist**2  # dV/dr
            f_vec = r_vec * f_magnitude  # dV/dr * dr/d(x,y,z)

            P_virial += np.dot(r_vec, f_vec)
    P_virial /= (3 * volume)
    #print(f"Current: P_KE: {P_KE:.2E}, P_virial: {P_virial:.2E}, {P_KE + P_virial:.2E}")#; import sys; sys.exit()
    return P_KE + P_virial

def Berendsen_barostat(R, V, L, tau_P):
    """
    Adjust volume and rescale positions to maintain constant pressure.
    Compute the scalor pressure using the virial equation.

    Args:
        R (np.array): N x 3 array of particle positions.
        V (np.array): N x 3 array of particle velocities.
        L (float): volume of the simulation box

    Returns:
        Updated R, V, L
    """
    volume = L ** 3
    P = compute_pressure(R, V, volume)
    dP = P - PRESSURE
    scale_factor = 1.0 + (dP / tau_P) * TIMESTEP

    # Rescale positions and velocities
    rescale_factor = scale_factor ** (1.0 / 3.0)
    R *= rescale_factor
    V *= rescale_factor
    L *= rescale_factor
    # print('update', L, rescale_factor, TIMESTEP, dP)
    return R, V, L

def MD(L0, barostat=None, Q=1.0, tau_P=0.1, num_steps=500):
    """
    Run MD simulation

    Args:
        barostat (str): "Berendsen", "" or None
        Q (float): Nose-Hoover thermostat parameter
        num_steps (int): Number of steps to simulate
    """
    # Data to monitor
    KEs, PEs, TEs = [], [], []
    Pressures, Volumes = [], []

    # Initialize system
    L = L0
    R = initialize_position(L)
    V = initialize_velocity(N)
    E, F = LJ_energy_forces(R, L)

    xi = 0.0  # Used by Nose-Hoover

    # MD propogation
    for step in range(num_steps):
        R, V, F, xi = Nose_Hoover_thermostat(R, V, F, xi, L, Q)
        if barostat == 'Berendsen':
            R, V, L = Berendsen_barostat(R, V, L, tau_P)

        # Compute PE, KE, and TE
        PE, _ = LJ_energy_forces(R, L)
        KE = 0.5 * np.sum(V**2) * MASS
        volume = L**3
        P = compute_pressure(R, V, volume)
        KEs.append(KE)
        PEs.append(PE)
        TEs.append(PE+KE)
        Volumes.append(volume)
        Pressures.append(P)

        if step % 2 == 0:
            s = np.exp(-xi * TIMESTEP)
            print(f"Step {step:6d}, PE: {PE:.5e} KE: {KE:.5e} E: {PE+KE:.5e} Vol: {L**3:.5e} s: {s:.4f}")
    return KEs, PEs, TEs, Volumes, Pressures


if __name__ == "__main__":

    # Parameters
    KB = 1.3806452 * 1e-23      # Boltzmann constant in J/K
    TEMPERATURE = 94.4          # in K
    PRESSURE = 101325           # in pascal
    EPSILON = 120.0 * KB        # in J (epsilon = 120 * k_B)
    SIGMA = 3.4 * 1e-10         # in meters (3.4 Å)
    MASS = 39.95*1.6747*1e-27   # mass of Argon atom in kg
    L0 = 10.229 * SIGMA         # cubic box side length
    N = 864                     # Number of atoms
    TIMESTEP = 1.0 * 1e-14      # time step in seconds

    results = []
    for barostat in ["Berendsen", None]:
        print(f"Simulation with {barostat} barostat")
        KEs, PEs, TEs, Vs, Ps = MD(L0, barostat=barostat, tau_P=0.01, num_steps=10)
        results.append((barostat, KEs, PEs, TEs, Vs, Ps))

