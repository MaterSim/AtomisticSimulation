# Week 6: MD Analysis II: Transport Processes

Transport properties describe how particles, energy, and momentum move within a system. In the lecture, we’ll discuss two key transport properties: diffusion (particle transport) and thermal conductivity (heat transport). We’ll again use argon gas as an example, but the methods are generalizable to other systems.

## 6.1 Diffusion

Diffusion constant $D$ describes how fast particles diffuse through the medium. In a simple random walk model of diffusion (based on Brownian motion), we assume that a particle is undergoing random, uncorrelated displacements over time. The mean square displacement (MSD) of a particle is a function of time, given by:

$$
\langle \Delta r^2(t) \rangle = 2d D t
$$

Where:

- $\langle \Delta r^2(t) \rangle$ is the mean square displacement (MSD) of the particle after time  t .
- $d$ is the number of spatial dimensions (e.g.,  d = 3  for 3D).
- $t$ is time.

Now, suppose we apply a small external force $F$ to the system (e.g., an electric field $E$ acting on a charged particle with charge $q$). The particle will respond to the applied force by moving with a drift velocity $v_d$ .

From Newton’s second law $F = m \frac{dv}{dt}$, causing the particle to accelerate. However, in a system with thermal motion (e.g., Brownian motion), the particle reaches a steady-state drift velocity $v_d$ due to the balance between the external force and the opposing forces from random collisions with the surrounding particles,

$$
v_d = \mu F
$$

Where $\mu$ is the mobility of the particle, defined as the proportionality constant between the $v_d$ and $F$.

At thermal equilibrium, the particle’s motion is governed by the thermal energy  $k_B T$. The energy imparted by the external force $F$ to move the particle is balanced by the random thermal energy. The thermal energy of a particle at temperature $T$  is given by the equipartition theorem:

$$
E_{\text{thermal}} = \frac{1}{2} m \langle v^2 \rangle = \frac{1}{2} k_B T
$$

Since thermal motion drives diffusion, and temperature is the main source of this motion, we expect the diffusion constant  D  to be related to temperature $T$. At steady state, when the system is in thermal equilibrium, the diffusion constant is related to the particle’s mobility through the **Einstein relation**:

$$
D = \mu k_B T = \frac{1}{2d} \frac{d}{dt} \langle \Delta r^2(t) \rangle
$$

In MD simulations, the MSD is computed as the time-averaged square of the particle displacements from their initial positions.

```python
def compute_msd(positions):
    """
    Compute the mean square displacement (MSD) from atomic positions over time.
    positions: array of shape (num_timesteps, num_atoms, 3)
    """
    num_timesteps = positions.shape[0]
    num_atoms = positions.shape[1]

    msd = np.zeros(num_timesteps)
    for t in range(1, num_timesteps):
        displacements = positions[t] - positions[0]  # Displacement from initial positions
        squared_displacements = np.sum(displacements**2, axis=1)  # Sum of squared displacements per atom
        msd[t] = np.mean(squared_displacements)  # Average over all atoms

    return msd

# Example usage
positions = np.random.randn(1000, 100, 3)  # Simulated positions: 1000 timesteps, 100 atoms in 3D
msd = compute_msd(positions)
time = np.linspace(0, 100, len(msd))

# Fit MSD to time to calculate diffusion constant
D = np.polyfit(time, msd, 1)[0] / 6  # Slope of MSD vs time gives 6*D for 3D system

plt.plot(time, msd, label='MSD')
plt.xlabel('Time (ps)')
plt.ylabel('MSD (Å²)')
plt.title(f'Diffusion Constant: {D:.3e} cm²/s')
plt.show()
```

## 6.2 The Green-Kubo Relation
In addition to Einstein relation, $D$ can also be derived from the velocity autocorrelation function (VACF).

$$
D = \frac{1}{d} \int_0^\infty \langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle \, dt
$$

It comes from linear response theory, which relates macroscopic transport coefficients to time correlations of microscopic quantities under a spontaneous fluctuations in thermal equilibrium. More generally, many transport processes can be described by the [Green-Kubo relation](https://en.wikipedia.org/wiki/Green–Kubo_relations). For any transport coefficient $\lambda$ (such as diffusion coefficient, thermal conductivity, or viscosity):

$$
\lambda = \int_0^\infty \langle J(0) J(t) \rangle \, dt
$$

- Diffusion constant .v.s velocity. The flux of velocity causes diffusion. If a particle “forgets” its initial velocity, the particle loses memory of its initial velocity, leading to faster diffusion.
- Thermal conductivity .v.s heat current. Applying a temperature gradient causes a heat current to flow, leading to thermal conductivity. The faster the decay of heat current autocorrelation, the lower the thermal conductivity.
- Viscosity .v.s stress tensors. In a fluid, applying a shear stress leads to a flow, related to viscosity. A high viscosity fluid (like honey) will have a slower decay of the stress autocorrelation function, meaning the fluid “remembers” shear stresses for a longer time compared to a low viscosity fluid (like water).

These relations are derived from linear response theory, which states that the response of a system to a small perturbation is proportional to the equilibrium fluctuations of microscopic quantities (e.g., heat current, velocity, or stress)

## 6.3 Thermal Conductivity 
Thermal conductivity $\kappa$ is a measure of a material’s ability to conduct heat. The Green-Kubo relation is used to calculate $\kappa$ from the heat current autocorrelation function (HCACF):

$$
\kappa = \frac{1}{k_B T^2 V} \int_0^\infty \langle J_Q(t) J_Q(0) \rangle dt
$$

where $J_Q(t)$ is the heat current, $k_B$ is the Boltzmann constant,  $T$ is temperature, and $V$ is the system volume.


```python
def compute_hcacf(heat_current):
    """
    Compute the heat current autocorrelation function (HCACF).
    heat_current: array of shape (num_timesteps, 3) containing heat current at each time step
    """
    num_timesteps = heat_current.shape[0]
    hcacf = np.zeros(num_timesteps)

    # Compute autocorrelation
    for t in range(num_timesteps):
        hcacf[t] = np.sum(np.dot(heat_current[0], heat_current[t]))  # Dot product of heat current at t=0 and t

    return hcacf

# Example usage
heat_current = np.random.randn(1000, 3)  # Simulated heat current data (1000 timesteps, 3D)
hcacf = compute_hcacf(heat_current)

# Integrate HCACF to calculate thermal conductivity
k_B = 1.380649e-23  # Boltzmann constant (J/K)
T = 300  # Temperature (K)
V = 1e-24  # Volume (m^3)
thermal_conductivity = np.trapz(hcacf) / (k_B * T**2 * V)

plt.plot(hcacf)
plt.xlabel('Time (ps)')
plt.ylabel('HCACF')
plt.title(f'Thermal Conductivity: {thermal_conductivity:.3e} W/mK')
plt.show()
```

## 6.4 Further disscusions

- Discuss how the MSD should scale linearly with time at longer time scales if the system has reached diffusive behavior.
- Explain the process of calculating the HCACF.
- Integrating the HCACF gives the thermal conductivity.
- Discuss the physical meaning of the results and how the thermal conductivity compares with known values for argon.

