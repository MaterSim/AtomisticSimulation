# Week 2: Thermostat under the NVT ensemble

## 2.1 Moltivation
### 2.1.1 The limitation of NVE ensemble
So far, we have learned how to run a NVE MD simulation for a periodic system from both programing and application points of view. In such simulations, the total energy E should be a constant with the propagation of time. This is called *microcanonical ensemble* in statistical physics. 
However, this set up is not really the truth for many practical simulations. It is more likely that the system would interact with the surrounding environment and exchange heat over the boundaries. 

### 2.1.2 Extension to NVE by allowing the heat exchange.
In a real scenario, we often want to study the system under a constant temperature condition, instead of constant energy.  This method is particularly useful for simulating systems in the *canonical ensemble* (constant number of particles, volume, and temperature, often denoted as NVT).

To maintain the system at a desired temperature by coupling it to an external *heat bath*. There would be several thermostat techniques. In real life, temperature is a measure of how fast the particles are moving on average. But in a computer simulation, you need a way to control this “temperature” to make sure it stays at the desired level throughout the simulation. A thermostat in molecular dynamics is a tool that helps you keep the temperature of your simulated system steady, just like how a thermostat in your house keeps the room temperature stable.

## 2.2 Different types of thermostats

To introduce the thermostat to a MD system, the trick is to modify the integrator. Currently, there exist several flavors of thermostat techniques. Perhaps the easiest way to rescale velocities to force the total kinetic energy  energy to be equal to $\frac{3NkT}{2}$ at every few steps. However, this kind of rescaling can definitely perturb the MD trajectory strongly and thus not recommended. Below I will show you a few better strategies. 

### 2.2.1 The Anderson thermostat

The main idea of Anderson thermostat is inspired by the observation of physical collisions between particles in the system and particles in the surrounding environment. After collisions, the particles in the system would change the velocities. These “collisions” ensure that the system exchanges energy with the environment, maintaining thermal equilibrium.

Thus we could periodically pick some particles and randomizes the velocities of some particles in the system. These randomizations mimic the effect of an external environment interacting with the particles, ensuring that the system’s temperature remains constant. Hence, we allow two types of MD integration rules in the actual code.

1. **Random particle selection and velocity assignment**. With a certain probability $\nu$, the velocity of each particle is reassigned by sampling from a Maxwell-Boltzmann distribution corresponding to the desired temperature $T$.
2. **Ordinary update**. If the particle’s velocity is not reassigned, it evolves according to the usual equations of motion (e.g., using the Verlet integration method).

In this technique, **collision Frequency** $$\nu$$ determines how often the particle velocities are randomized (following a [Poisson Distribution](https://en.wikipedia.org/wiki/Poisson_distribution)). A higher $\nu$ means more frequent collisions (interaction) with the heat bath, leading to stronger coupling to the temperature bath. We should choose a $$\nu$$ so that velocity reassignment happens at an appropriate rate to maintain the desired temperature without overly disrupting the natural dynamics of the system.

The Anderson thermostat is relatively simple to implement, requiring only the addition of random velocity reassignment at each time step. However, it may not reflect the real dynamics. Since velocities are randomly reassigned, the resulting particle trajectories may not correspond to those in a real physical system where energy exchange occurs through physical interactions. This is particularly true for a periodic  system without the clear definition of boundary. In addition, one needs to play with the $\nu$ values.


The algorithm can be described as follows.

```python
import numpy as np

def anderson_thermostat(V, T, nu):
    sigma_v = np.sqrt(k_B * T / MASS)
    for i in range(num_particles):
        if np.random.rand() < nu * TIMESTEP:
            V[i] = np.random.normal(0, sigma_v, 3)
    return V

# Initialization
T = 300 # in K
R = Initialize_positions()
V = Initialize_velocities()
F = calculate_forces(R)

# Main MD loop
for step in range(num_steps):

    # Update R, V, F using Verlet integration
    R += V * TIMESTEP + 0.5 * F * TIMESTEP**2 / MASS
    F_new = calculate_forces(R)
    V += 0.5 * (F + F_new) * TIMESTEP / MASS
    F = F_new

    # Apply Anderson thermostat to randomly reassign V
    V = anderson_thermostat(V, T, nu)
```

### 2.2.2 The Langevin thermostat 
The Langevin thermostat maintains the temperature of a system while also modeling the effects of friction and random forces, similar to those that might be encountered in a viscous fluid. The basic idea it to modify the equations of motion by adding two additional terms to the standard Newtonian dynamics:

$$
m_i \frac{d\mathbf{v}_i}{dt} = \mathbf{F}_i - \gamma m_i \mathbf{v}_i + \mathbf{R}_i(t)
$$

* Frictional Force $\gamma m_i \mathbf{v}_i$: the damping effect of the environment, which tends to slow down the particles.

* Random Force $\mathbf{R}_i(t)$ : the random collisions with particles from the heat bath, which cause the particles to move randomly, maintaining the system’s temperature. These kicks help maintain the temperature of the system by continuously injecting energy into the system.

A typical value of $\gamma$ used in many MD simulations is around $\gamma = 0.1 \text{ps}^{-1}$. This value provides a good balance between maintaining temperature control and preserving realistic dynamics. The system is weakly coupled to the heat bath, ensuring that it can sample the canonical ensemble without heavily damping the natural motion of the particles.

In MD simulations, the Langevin equation is integrated with the Velocity Verlet algorithm, modified to include the Langevin terms. A simple update rule for the velocity might look like this:

$$
\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i(t) + \frac{\mathbf{F}_i(t)}{m_i} \Delta t - \gamma \mathbf{v}_i(t) \Delta t + \mathbf{R}_i \sqrt{\Delta t}
$$

```python
import numpy as np

def langevin_thermostat(V, F, gamma):
	# Update R, V, F using Verlet integration
    R += V * TIMESTEP + 0.5 * F * TIMESTEP **2 / MASS

	# Update velocities with deterministic part
    V += 0.5 * F * TIMESTEP / MASS

    # Apply friction and random force (stochastic part)
    V += -gamma * V * TIMESTEP
	sigma = np.sqrt(2 * gamma * kB * T / MASS)
    V += np.random.normal(0, sigma, V.shape) * np.sqrt(dt)

	# Update forces and velocities
	F_new = calculate_forces(positions)
	V += 0.5 * (F + F_new) / MASS * TIMESTEP
	F = F_new
    return R, V, F

# Initialization
T = 300 # in K
R = Initialize_positions()
V = Initialize_velocities()
F = calculate_forces(R)

# Main MD loop
for step in range(num_steps):
	R, V, F = langevin_thermostat(R, V, F, L, gamma)
```


### 2.2.3 The Nosé-Hoover thermostat
If we want to avoid the use of brute-force velocity reassignment, a gentler approach is to control the temperature by coupling the system to an additional degree of freedom, which acts as a “thermal reservoir” that exchanges energy with the system.
Nosé introduced a method where the system’s Hamiltonian is extended by adding an artificial variable that represents the thermal reservoir. This approach ensures that the system’s temperature fluctuates around a desired value, allowing it to correctly sample the canonical ensemble. Then, Hoover reformulated the equations of motion to include a friction term ($\xi$) that dynamically adjusts the particle velocities to maintain the target temperature, simplified Nosé’s method, making it more efficient for MD simulations. See an extended discussion in the Appendix.

In the Nosé-Hoover thermostat, the velocity is updated via the following term

$$
\frac{d\mathbf{v}_i}{dt} = \frac{\mathbf{F}_i}{m_i} - \xi \mathbf{v}_i
$$

where $\mathbf{v}_i$ is the velocity of particle  $i$ , $\mathbf{F}_i$ is the force acting on particle $i$ , $m_i$ is the mass of the particle, and $\xi$ is the friction coefficient or thermostat variable.
	
Then $\xi$ is updated as follows

$$
\frac{d\xi}{dt} = \frac{1}{Q} \left(\sum_i \frac{m_i \mathbf{v}_i^2}{3Nk_BT} - 1\right)
$$

where $Q$ is the “thermal inertia” parameter (which controls how strongly the system is coupled to the thermostat),  $T$ is the target temperature, $N$ is the number of particles, and  $k_B$  is the Boltzmann constant.


```python
import numpy as np

# Thermostat variables
Q = 100.0                      # Thermal inertia parameter (tune this for stability)
xi = 0.0

# Initialization
R = Initialize_positions()
V = Initialize_velocities()
F = calculate_forces(R)

# Main MD loop
for step in range(num_steps):

    # Verlet-like integration
    R += V * TIMESTEP + 0.5 * F * TIMESTEP**2 / MASS
    F_new = calculate_forces(R)

    # Update velocities
    V += 0.5 * (F + F_new) * TIMESTEP / MASS
    V *= (1 - 0.5 * xi * TIMESTEP) / (1 + 0.5 * xi * TIMESTEP)

    # Update the Nosé-Hoover thermostat variable
    kE= 0.5 * np.sum(mass * V**2)
    xi += dt * (2 * kE / (3 * N * k_B * T) - 1) / Q
```

### 2.2.4 Extended Discussions in Ensembles and Nosé-Hoover thermostat
Complete the reading in [Appendix-W2](https://github.com/qzhu2017/AtomisticSimulation/blob/main/Appendix/W2_NoseHoover.pdf).

## 2.3 Full code to run NVT simulation
1. Complete the codes in [Colab](https://colab.research.google.com/drive/1x8FFEDrvmThUQGhfVCJd9LZka0aO1zhe?usp=sharing)
2. Modify your code with different $nu$, $gamma$ and $Q$ parameters for different thermostats and monitor the progress of kinetic energies.
3. Debug the code [lec_02_langevin_debug.py](https://github.com/qzhu2017/AtomisticSimulation/blob/main/Codes/lec_02_langevin_debug.py)

