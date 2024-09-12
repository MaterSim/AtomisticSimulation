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
        P (float): P_virial
    """
    N = len(R)                # Number of atoms as a scalor
    F = np.zeros_like(R)      # forces [N, 3] as a 2D array
    PE = 0.0                  # total potential energy as a scalor
    P_virial = 0.0            # if stress is on

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
            P_virial += np.dot(r_vec, force_vec)

    return PE, F, P_virial


def Nose_Hoover_thermostat(R, V, F, xi, L, Q):
    """
    Nose-Hoover thermostat
    """
    # Update R
    R += V * TIMESTEP + 0.5 * F/MASS * TIMESTEP ** 2
    R = R % L

    # Update forces
    PE, F_new, P_virial = LJ_energy_forces(R, L)
    V += 0.5 * (F + F_new) / MASS * TIMESTEP
    V *= (1 - 0.5 * xi * TIMESTEP) / (1 + 0.5 * xi * TIMESTEP)

    # Update xi
    KE = 0.5 * np.sum(V**2) * MASS
    xi += TIMESTEP * (2 * KE / (3 * len(R) * KB * TEMPERATURE) - 1) / (Q*MASS)

    # Update forces
    F = F_new
    return R, V, F, xi, PE, P_virial

def compute_pressure(R, V, volume, P_virial):
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
    # kinetic contribution to the pressure
    P_KE = MASS * np.sum(V**2)#; print("KE", 0.5 * np.sum(V**2) * MASS, "Volume", volume)
    P_KE /= (3 * volume)
    P_virial /= (3 * volume)
    #print(f"P_KE: {P_KE:.2E}, P_virial: {P_virial:.2E}, {P_KE + P_virial:.2E}")#; import sys; sys.exit()
    return P_KE + P_virial

def Berendsen_barostat(R, V, L, P_virial, tau_P):
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
    P = compute_pressure(R, V, volume, P_virial)
    dP = P - PRESSURE
    scale_factor = 1.0 + (dP / tau_P) * TIMESTEP

    # Rescale positions and velocities
    rescale_factor = scale_factor ** (1.0 / 3.0)
    R *= rescale_factor
    V *= rescale_factor
    L *= rescale_factor
    return R, V, L, P

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
    E, F, P_virial = LJ_energy_forces(R, L)

    xi = 0.0  # Used by Nose-Hoover

    # MD propogation
    for step in range(num_steps):
        R, V, F, xi, PE, P_virial = Nose_Hoover_thermostat(R, V, F, xi, L, Q)
        if barostat == 'Berendsen':
            R, V, L, P = Berendsen_barostat(R, V, L, P_virial, tau_P)
        else:
            P = compute_pressure(R, V, L**3, P_virial)

        # Compute PE, KE, and TE
        KE = 0.5 * np.sum(V**2) * MASS
        volume = L**3

        KEs.append(KE)
        PEs.append(PE)
        TEs.append(PE+KE)
        Volumes.append(volume)
        Pressures.append(P)

        if step % 10 == 0:
            s = np.exp(-xi * TIMESTEP)
            print(f"Step {step:6d}, PE: {PE:.5e} KE: {KE:.5e} P: {P:.5e} Vol: {volume:.5e} s: {s:.4f}")
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
    for params in [("Right Volume", L0, None, None),
                   ("Large Volume", 1.05*L0, None, None),
                   ("Right Volume", L0, "Berendsen", 0.001),
                   ("Large Volume", 1.05*L0, "Berendsen", 0.001)]:

        (tag, L, barostat, tau_P) = params
        print(f"Simulation with {barostat} barostat {tag}")
        KEs, PEs, TEs, Vs, Ps = MD(L, barostat=barostat, tau_P=tau_P, num_steps=2000)
        results.append((tag, barostat, KEs, PEs, TEs, Vs, Ps))

    fig, axs = plt.subplots(2, len(results), figsize=(16, 6))
    for i, result in enumerate(results):
        (tag, barostat, KEs, PEs, TEs, Vs, Ps) = result
        if barostat is None: barostat = "NVT"
        axs[0, i].set_title(tag + '-' + barostat)
        axs[0, i].plot(KEs, label='KE')
        axs[0, i].plot(PEs, label='PE')
        axs[0, i].plot(TEs, label='Total')
        axs[0, i].set_ylim([-1.0e-17, 0.3e-17])
        #axs[0, i].set_xlabel("Timesteps")

        ax_vol = axs[1, i]
        ax_pre = ax_vol.twinx()  # Create a twin y-axis for KEs
        ax_vol.plot(Vs, color='b') #, label=barostat+'_Volume')
        ax_pre.plot(Ps, color='r') #, label=barostat+'_Pressure')
        #ax_vol.set_xlabel("Timesteps")

        if i == 0:
            axs[0, i].set_ylabel("Energy")
            ax_vol.set_ylabel("Volume", color='b')
            ax_vol.tick_params(axis='y', labelcolor='b')
        else:
            ax_vol.set_yticks([])
            axs[0, i].set_yticks([])

        if i + 1 == len(results):
            ax_pre.set_ylabel("Pressure", color='r')
            ax_pre.tick_params(axis='y', labelcolor='r')
        else:
            ax_pre.set_yticks([])

        ax_vol.set_ylim([3.6e-26, 5.1e-26])
        ax_pre.set_ylim([-3.0e+8, 1.2e+8])
        ax_vol.spines['left'].set_color('b')  # Color left y-axis (volume) green
        ax_pre.spines['right'].set_color('r')  # Color right y-axis (KE) blue

    for ax in axs.flat:
        ax.legend()
    plt.savefig('results.png')

    """
    Below are the codes to check the correctness based on ASE.
    You are welcome to check the results

    # ASE
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    R = initialize_position(L0) * 1e+10  # in angstrom
    cell = np.eye(3) * L0 * 1e+10
    argon = Atoms('Ar864', positions=R, cell=cell, pbc=True)
    lj_calculator = LennardJones(sigma=3.4, epsilon=EPSILON/1.60218e-19)
    argon.set_calculator(lj_calculator)
    MaxwellBoltzmannDistribution(argon, temperature_K=TEMPERATURE)

    KE = argon.get_kinetic_energy() * 1.60218e-19
    PE = argon.get_potential_energy() * 1.60218e-19
    stress = argon.get_stress(include_ideal_gas=True) * 1.60218e+11
    print(f"KE: {KE:.2E} J")
    print(f"PE: {PE:.2E} J")
    #print("Stress: ", stress)

    volume = argon.get_volume()
    P_KE = 2/3*argon.get_kinetic_energy()/volume * 1.60218e+11 # Pa
    P_virial = -np.trace(argon.get_stress(voigt=False).reshape(3, 3)) / 3 * 1.60218e+11
    P = P_KE + P_virial
    print(f"Final {P_KE:.2E}, {P_virial:.2E}, {P:.2E}, Vol {volume:.2E}")
    """
