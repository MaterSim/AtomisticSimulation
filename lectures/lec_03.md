# Week 3: Barostat under the NPT ensemble

## 3.1 Moltivation
The introduction of temperature control in the previous lecture allows us to simulation the system to maintain a constant temperature through a heat exchange with the reservior. While this *canonical ensemble* enables more realistic modelling of the real world, it still has some limitation. Just imagine that you want to model a periodic system under different temperatures, it is natural to think about the volume of the system is subject to change due to the thermal expansion effects.
Within the NVT ensemble, you need to manually adjust the volume of the system when initializing the positions for each different temperture. This can be very tedious and requires a trial and error iterations. Thus, we would like to seek a better solution to let the system adjust the volume by itself during the MD simulation. 

Similar to how a thermostat maintains constant temperature, a barostat adjusts the simulation box dimensions and/or particle positions to ensure that the system stays at the desired pressure. This is crucial for simulating systems in ensembles like NPT, where pressure fluctuations need to be controlled. Below, I will introduce two barostat techniques.

## 3.2 Barostat Techniques
### 3.2.1 Berendsen Barostat
This is a simple barostat that rescales the simulation box gradually toward the target pressure. To implement a barostat, the key idea is to adjust the simulation box size in response to the difference between the current pressure and the target pressure. This is done by scaling the box dimensions and particle positions, and updating the system’s volume accordingly.

1. The instantaneous pressure $P$ in an MD simulation is calculated using the virial equation. It includes contributions from the kinetic energy and the virial of the system (related to particle interactions):

$$
P = \frac{N k_B T}{V} + \frac{1}{3V} \sum_{i \lt j} r_{ij} \cdot \mathbf{F}_{ij}
$$

Where:

- $N$ is the number of particles.
- $k_B$ is the Boltzmann constant.
- $T$ is the temperature.
- $V$ is the volume of the simulation box.
- $\mathbf{r}_{ij}$ is the position vector between particles $i$ and $j$.
- $\mathbf{F}_{ij}$ is the force acting on particle $i$ due to particle $j$.

2. Compute pressure difference. At each time step, calculate the difference between the current and target pressures.

3. Adjust the simulation box volume. For isotropic pressure control (same scaling in all directions), the new volume is updated by:

$$
V_{new} = V_{old} \left( 1 + \frac{\Delta P}{\tau_P} \cdot dt \right)
$$

Where:
- $\tau_P$ is a time constant controlling the pressure coupling strength.
- $dt$  is the time step.
- $V_{old}$ is the current volume.

4. Rescale the positions. The positions of all particles are scaled accordingly to maintain their relative distances within the simulation box. For isotropic scaling, each position  $\mathbf{r}$ is rescaled:

$$
r_{new} = r_{old} \cdot \left( \frac{V_{new}}{V_{old}} \right)^{1/3}
$$

The pseudo code should look like the following
```python
import numpy as np

P_target = 1.0  # Target pressure (arbitrary units)
tau_P = 0.5  # Pressure coupling constant (barostat relaxation time)
kb = 1.38e-23  # Boltzmann constant in J/K

def compute_temperature(velocities, masses):
    """
    Compute the temperature of the system from particle velocities.

    Parameters:
    velocities (np.array): N x 3 array of particle velocities.
    masses (np.array): Array of particle masses.

    Returns:
    float: Temperature of the system.
    """
    KE = 0.5 * np.sum(masses[:, None] * velocities**2)  # Sum of kinetic energies
    temperature = (2 * KE) / (3 * len(velocities))  # Ideal gas temperature relation
    return temperature

def compute_virial_pressure(positions, velocities, forces, V):
    """
    Compute the pressure using the virial equation in the Berendsen barostat.

    Parameters:
    positions (np.array): N x 3 array of particle positions.
    velocities (np.array): N x 3 array of particle velocities.
    forces (np.array): N x 3 array of interatomic forces.
    V (float): Volume of the simulation box.
    
    Returns:
    float: Computed pressure of the system.
    """
    
    N = len(positions)  # Number of particles
    
    # Step 1: Compute the temperature from velocities
    T = compute_temperature(velocities, masses)
    
    # Step 2: Compute the kinetic contribution to the pressure
    P_KE = (N * kB * T) / V
    
    # Step 3: Compute the virial contribution to the pressure
    P_virial = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r_ij = positions[i] - positions[j]  # Displacement vector between particles i and j
            F_ij = forces[i]  # Force on particle i due to particle j (assumed already calculated)
            P_virial += np.dot(r_ij, F_ij)  # Dot product r_ij · F_ij
    
    P_virial /= (3 * V)  # Virial term divided by 3V
    
    # Total pressure is the sum of kinetic and virial contributions
    P = P_kinetic + P_virial
    
    return P

# Barostat function (Berendsen type)
def apply_barostat(positions, velocities, forces, V, P_target, tau_P, dt):
    """Adjust volume and rescale positions to maintain constant pressure."""
    # Calculate current pressure
    P = compute_virial_pressure(positions, velocities, V, T)
    
    # Calculate volume scaling factor
    dP = P - P_target
    scale_factor = 1.0 + (dP / tau_P) * dt
    
    # Update volume
    V_new = V * scale_factor

    # Rescale positions according to the new volume
    rescale_factor = (V_new / V) ** (1.0 / 3.0)
    positions *= rescale_factor

    # Optionally rescale velocities
    velocities *= rescale_factor

    return positions, velocities, V_new
```

This is a relatively simple method, where the system’s volume is gradually rescaled to match the target pressure. It does not rigorously conserve the ensemble, but it is computationally efficient and often used for equilibration runs for the simulation of isotropic systems like liquid.


### 3.2.2 Parrinello-Rahman Barostat 
The Parrinello-Rahman barostat is a more advanced method for controlling pressure in MD simulations, particularly useful when the system undergoes anisotropic volume changes. Unlike Berendsen barostat that scales the simulation box isotropically, the Parrinello-Rahman barostat allows both the shape and size of the simulation box to change. This is especially important in simulations of materials under strain, phase transitions, or when dealing with anisotropic systems like crystals.

1. To enable this barostate, we first represent simulation box by a matrix $\mathbf{h}$ that defines the three lattice vectors of the simulation box. This matrix allows for changes in both the box size and shape.

$$
\mathbf{h} =
\begin{pmatrix}
\mathbf{a}_x & \mathbf{b}_x & \mathbf{c}_x \\
\mathbf{a}_y & \mathbf{b}_y & \mathbf{c}_y \\
\mathbf{a}_z & \mathbf{b}_z & \mathbf{c}_z
\end{pmatrix}
$$

Here, $\mathbf{a}$, $\mathbf{b}$, and $\mathbf{c}$ are the lattice vectors.

2. Pressure Tensor: The barostat works with the full pressure tensor $\mathbf{P}$, which describes how pressure acts differently along different directions. The pressure tensor $\mathbf{P}$ can be computed from the system’s kinetic energy and virial:

$$
\mathbf{P} = \frac{1}{V} \left( \sum_i m_i \mathbf{v}_i \otimes \mathbf{v}_i + \sum_i \mathbf{r}_i \otimes \mathbf{F}_i \right)
$$

Where:

- $V$ is the volume of the simulation box.
- $\mathbf{v}_i$ and $\mathbf{F}_i$  are the velocity and force on particle $i$ , respectively.
- $m_i$ is the mass of particle $i$.


3. Additional variable strain Rate Tensor \mathbf{W}. The time evolution of the box matrix $\mathbf{h}$ is governed by the equation:

$$
\dot{\mathbf{h}} = \mathbf{h} \cdot \mathbf{W}
$$

**Where $\mathbf{W}$ is the strain rate tensor, which determines how the box evolves over time.** 

$$
\frac{d \mathbf{W}}{dt} = \frac{1}{Q} \left( \mathbf{P} - P_{\text{target}} \mathbf{I} \right)
$$

Here, $Q$ is the fictitious barostat mass, and $P_{\text{target}} \mathbf{I}$  is the target pressure tensor.

4. Update Particle Positions and Box Matrix

Once $\mathbf{h}$ is updated, the particle positions need to be rescaled by the new box matrix. The rescaled positions are calculated as:

$$
\mathbf{r}_i = \mathbf{h} \cdot \mathbf{s}_i
$$


In short, this approach introduces a few additional variables:
- $\mathbf{h}$ : The simulation box matrix, which evolves over time and controls both the size and shape of the box.
- $Q$ : The fictitious mass associated with the barostat, controlling the rate of volume and shape changes.
- $\mathbf{W}$ : The strain rate tensor, which governs how the box matrix changes over time.

These variables allow the Parrinello-Rahman barostat to apply pressure anisotropically, enabling the box to deform naturally while maintaining the target pressure in the system.


Below is a pseduo code to achieve Parrinello-Rahman barostat.

```python
import numpy as np

# Constants
dt = 0.001  # Time step
P_target = np.eye(3) * 1.0  # Target pressure tensor (isotropic)
Q = 10.0  # Fictitious mass for the barostat, analogical to the Q in Nose-Hoover thermostat
N = 100  # Number of particles
kb = 1.38e-23  # Boltzmann constant
T = 300  # Temperature
V = 1.0  # Initial volume
positions = np.random.randn(N, 3)
velocities = np.random.randn(N, 3)
forces = np.zeros_like(positions)  # Placeholder for forces
h = np.eye(3) * V ** (1/3)  # Initial box matrix (cubic box)

def compute_pressure_tensor(positions, velocities, forces, V, kB):
    """
    Compute the internal pressure tensor using the virial equation.
    
    Parameters:
    positions (np.array): N x 3 array of particle positions.
    velocities (np.array): N x 3 array of particle velocities.
    forces (np.array): N x 3 array of interatomic forces.
    V (float): Volume of the simulation box.
    
    Returns:
    np.array: 3 x 3 pressure tensor.
    """
    
    # Number of particles
    N = len(positions)
    
    # Compute the kinetic energy contribution to the pressure
    KE = np.sum(0.5 * masses[:, None] * velocities**2)
    T = (2 * kinetic_energy) / (3 * N * kB)
    P_kE = N * kB * T / V
    
    # Compute the virial contribution to the pressure tensor
    P_virial = np.zeros((3, 3))
    for i in range(N):
        for j in range(i + 1, N):
            r_ij = positions[i] - positions[j]  # Displacement vector
            F_ij = forces[i]  # Force on particle i
            P_virial += np.outer(r_ij, F_ij)
    
    # Average pressure tensor by dividing by the volume
    P_virial /= volume
    
    # Total pressure tensor is the sum of kinetic and virial contributions
    P_total = P_KE* np.eye(3) + P_virial
    return P_total

# Parrinello-Rahman barostat step
def parrinello_rahman_barostat(h, positions, velocities, forces, P_target, Q, dt):
    """Update the box and rescale positions using Parrinello-Rahman barostat."""
    # Compute current pressure tensor
    V = np.linalg.det(h)  # Current volume
    P = calculate_pressure_tensor(positions, velocities, forces, V)

    # Compute strain rate tensor (dW/dt)
    W_dot = (P - P_target) / Q

    # Update the box matrix h
    h_new = h + h @ W_dot * dt

    # Rescale positions
    s_positions = np.linalg.inv(h) @ positions.T  # Transform positions to fractional coordinates
    positions = h_new @ s_positions  # Transform back to new coordinates

    return h_new, positions.T
```


A more complete discussion can be found [here](https://computecanada.github.io/molmodsim-md-theory-lesson-novice/08-barostats/index.html).

## 3.3 Full code to run NPT simulation
Complete the codes in [Colab](https://colab.research.google.com/drive/1x8FFEDrvmThUQGhfVCJd9LZka0aO1zhe?usp=sharing)


