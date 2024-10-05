# 6 MD Analysis II: Transport Processes
Transport properties describe how particles, energy, and momentum move within a system. In the lecture, we’ll discuss two key transport properties: diffusion (particle transport) and thermal conductivity (heat transport). We’ll again use argon gas as an example, but the methods are generalizable to other systems.

## 6.1 Diffusion

According to [Fick's 2nd law](https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion), the diffusion equation is a partial differential equation that describes how a substance spreads over time. For a one-dimensional system under equilibrium, it is given by:

$$
\frac{\partial C(x,t)}{\partial t} = D \frac{\partial^2 C(x,t)}{\partial x^2}
$$

where:

- $C(x, t)$ is the concentration of particles at position $x$ and time $t$.
- $D$ is the diffusion constant (or diffusivity), which characterizes how fast particles diffuse.

The equation states that the rate of change of concentration over time $\frac{\partial C}{\partial t}$ is proportional to the second spatial derivative of the concentration $\frac{\partial^2 C}{\partial x^2}$.

The solution to the diffusion equation gives the probability distribution for the particle’s position over time. In the case of one-dimensional diffusion starting from a point, the solution is a Gaussian distribution:

$$
C(x,t) = \frac{1}{\sqrt{4\pi D t}} \exp\left( -\frac{x^2}{4Dt} \right)
$$

The Gaussian form suggests that:
1. The peak of the distribution (the center, at $x = 0$) remains at the origin because particles are assumed to start there.
2. The spread of the distribution increases with time, meaning that particles are more likely to be found farther away from the origin as time goes on.

To find the mean squared displacement (MSD), we calculate the expected value of $x^2$ with respect to this distribution. The MSD in one dimension is:

$$
\langle x^2(t) \rangle = \int_{-\infty}^{\infty} x^2 C(x,t) \, dx = 2Dt
$$

Hence, we end up with the following relation:

$$
\langle \Delta r^2(t) \rangle = 2d D t
$$

where:

- $\langle \Delta r^2(t) \rangle$ is the mean squared displacement (MSD) of the particle after time $t$.
- $d$ is the number of spatial dimensions (e.g., $d = 3$ for 3D).

It is also convenient to write it in the derivative form:

$$
\frac{\partial \langle \Delta r^2(t) \rangle}{\partial t} = 2d D
$$

<!Now, suppose we apply a small external force $F$ to the system (e.g., an electric field $E$ acting on a charged particle with charge $q$). The particle will respond to the applied force by moving with a drift velocity $v_d$ .

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
D = \mu k_B T = \frac{\langle \Delta r^2(t) \rangle}{2dt}
$$
|>

In MD simulations, the MSD is computed as the time-averaged square of the particle displacements from their initial positions. From the MSD, we can then perform linear regression to find the slope of the MSD($t$) curve to determine the diffusion constant.

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
### 6.2.1 The Alternative Expression of MSD

Let's introduce an alternative definition of the displacement of an atom.

$$
\Delta x(t) = \int_0^t dt\prime v_x(t\prime)
$$

Hence, the MSD can be expressed as:

$$
\begin{align*}
\langle \Delta x^2(t) \rangle &= \left\langle \left( \int_0^t dt{\prime}  v_x(t{\prime}) \right)^2 \right\rangle \\
&= \left\langle\int_0^t dt{\prime} v_x(t{\prime}) \int_0^t dt{\prime}{\prime} v_x(t{\prime}{\prime})\right\rangle\\
& = \int_0^t dt{\prime} \int_0^t dt{\prime}{\prime}  \langle v_x(t{\prime}) v_x(t{\prime}{\prime}) \rangle \\
&= 2 \int_0^t dt{\prime} \int_0^{t{\prime}} dt{\prime}{\prime}  \langle v_x(t{\prime}) v_x(t{\prime}{\prime}) \rangle
\end{align*}
$$

Note that the factor of 2 arises because the integral over the entire $t$-range is equivalent to twice the integral over the *half-space* where $t{\prime} > t{\prime}{\prime}$ (or vice versa, due to symmetry). Hence, $D$ can be computed by the VACF as introduced in the previous lecture.

$$
D = \frac {\partial \langle \Delta r^2(t) \rangle}{2d \partial t} = \frac{1}{d} \int_0^\infty \langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle  dt
$$

## 6.2.2 The General Green-Kubo Relation

In fact, the previous expression between $D$ and VACF is a special case. According to the linear response theory, which relates macroscopic transport coefficients to time correlations of microscopic quantities under spontaneous fluctuations in thermal equilibrium, many transport processes can be described by the [Green-Kubo relation](https://en.wikipedia.org/wiki/Green–Kubo_relations). For any transport coefficient $\lambda$ (such as diffusion coefficient, thermal conductivity, or viscosity):

$$
\lambda = \int_0^\infty \langle J(0) \cdot J(t) \rangle  dt
$$
- **Diffusion constant vs. velocity**: The flux of velocity causes diffusion. If a particle quickly *forgets* its initial velocity, it leads to faster diffusion as the particle loses memory of its initial velocity.
- **Thermal conductivity vs. heat current**: Applying a temperature gradient causes a heat current to flow, leading to thermal conductivity. The faster the decay of heat current autocorrelation, the greater the thermal conductivity.
- **Viscosity vs. stress tensors**: In a fluid, applying a shear stress leads to a flow, related to viscosity. A high viscosity fluid (like honey) will have a slower decay of the stress autocorrelation function, meaning the fluid *remembers* shear stresses for a longer time compared to a low viscosity fluid (like water).
These relations are derived from linear response theory, which states that the response of a system to a small perturbation is proportional to the equilibrium fluctuations of microscopic quantities (e.g., heat current, velocity, or stress).

## 6.3 Thermal Conductivity 
Thermal conductivity $\kappa$ is a measure of a material’s ability to conduct heat. The Green-Kubo relation is used to calculate $\kappa$ from the heat current autocorrelation function (HCACF):

$$
\kappa = \frac{1}{k_B T^2 V} \int_0^\infty \langle J(0) \cdot J(t) \rangle dt
$$

where $J(t)$ is the heat current, $k_B$ is the Boltzmann constant,  $T$ is temperature, and $V$ is the system volume.


The heat current is defined as the rate at which energy crosses a unit area per unit time. Thue, $J$ can be computed using the [following relation](https://docs.lammps.org/compute_heat_flux.html)

$$
\mathbf{J} = \frac{1}{V} \left[ \sum_i e_i \mathbf{v}_i - \sum_i S_i \mathbf{v}_i \right]
$$



where 
- $e_i$ is the per-atom energy (potential and kinetic). 
- $S_i$ is the per-atom stress tensor
- $v_i$ is the velocity

This relation can be understood as a microscopic version of the derivative form of 1st law ($\dot{Q}=\dot{E}-\dot{W}$).

And the HCACF can be numerically evaluated as

$$
\langle J(t) \cdot J(0) \rangle_{\text{eq}} = \frac{1}{N} \sum_{i=1}^{N} \langle J_i(t) \cdot J_i(0) \rangle
$$

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
## 6.4 Further Discussions

Finally, it is important to note that the techniques to compute diffusion coefficients or thermal conductivity are indirect. They rely on the statistical properties of the system. They are essentially rooted in the [Fluctuation-Dissipation Theorem](https://en.wikipedia.org/wiki/Fluctuation–dissipation_theorem), a principle underlying the connection between the fluctuations in a system at equilibrium and its response to perturbations. The theorem states that the response of a system to a small perturbation can be related to the equilibrium fluctuations of observables. Rather than measuring these properties directly, which can be difficult or impossible in practice, we can derive them from measurable fluctuations in the system.

That said, it is possible to use so-called nonequilibrium MD to directly measure transport properties like diffusion coefficients and thermal conductivity, but there are limitations and challenges associated with this approach. The system might take a long time to reach a steady state, or the steady state might not be easy to achieve in practice.

- Discuss how the MSD should scale linearly with time at longer time scales if the system has reached diffusive behavior.
- Explain the process of calculating the HCACF.
- Integrating the HCACF gives the thermal conductivity.
- Discuss the physical meaning of the results and how the thermal conductivity compares with known values for argon.

