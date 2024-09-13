# MEGR 7090/8090: Atomistic Simulation in Materials Modeling
 
## Course Introduction
This is a 3-credit course requires three hours of classroom or direct faculty instruction and six hours of out-of-class student work for the equivalent of approximately 15 weeks. Out-of-class work may include but is not limited to: required reading; homework; studying for quizzes and exams; research; written assignments; and project design, simulation, testing and demonstration.

## Instructor: Qiang Zhu (Battcave 114, qzhu8@uncc.edu)

## Textbooks
- *Understanding molecular simulation from algorithms to applications*, By Daan Frankel and Berend Smit, [3rd Edition](https://shop.elsevier.com/books/understanding-molecular-simulation/frenkel/978-0-323-90292-2)
- *Electronic structure*, By Richard. Martin, [2nd Edition](https://www.cambridge.org/core/books/electronic-structure/DDFE838DED61D7A402FDF20D735BggC63A)

The lecture notes were made based on these two excellent books. However, the hard copies of textbooks are not strictly required. We will also keep updating this lecture notes and provide more open access video or text tutorials throughout the course.

## Course Description
This course aims to use the atomistic computer simulation to model and understand the properties of real materials and their accompanying process and phenomena. It will primarily focus on two approaches: molecular dynamics and electronic structure calculation based on density functional theory. Some typical examples, codes, analytical tools will be also covered in this course. 

The expected outcomes include: 
- Understand the fundamental of Molecular dynamics simulation and its connection with statistical physics
- Apply the molecular dynamics simulation technique to model the physical process in real materials
- Understand the concept of electronic structure simulation based on density functional theory 
- Use the available software LAMMPS and VASP to compute material’s properties

## Tenative schedules

### I: Molecular dynamics simulation
- Week 1: Motivating example 1: Argon under the NVE ensemble
- Week 2: Motivating example 2: Thermostat NVT ensemble
- Week 3: Motivating example 3: Simulation of solids under the NPT ensemble 
- Week 4: Introduction to the LAMMPS package
- Week 5: MD Analysis I: structural characterization (RDF), degree of order
- Week 6: MD Analysis II: transport processes
- Week 7: Selected presentations from students

### II: Electronic structure calculation
- Week 8: Gentle introduction to Density functional theory 
- Week 9: Motivating example 4: DFT treatment of H2 molecule
- Week 10: From molecule to the periodic system
- Week 11: Introduction to VASP
- Week 12: Band structure analysis
- Week 13: Phonon calculation from both classical force field and DFT
- Week 14: Selected presentations from students 


# Week 1: Simulation under the NVE ensemble

## 1.1 Overview
### 1.1.1 Prehistory of Computer Simulation:
* The Los Alamos MANIAC (Mathematical Analyzer, Numerical Integrator, and Computer) became operational in 1952. This event marks a significant milestone in the history of computing. Nicholas Metropolis, as the most notable early user, developed the Monte Carlo method, a statistical technique that utilizes random sampling to solve complex mathematical problems. 

* The launch of computer also opened the door for the study of many fundamental problems. Most of the material systems consist of many atoms or molecules, how can we infer the properties of such systems? In the past, people have to to do it either analytically (e.g., thermodynamics and statiscal mechanics have been developed to study some classical systems such as ideal gas, Ising Model, ferromagentic phase transition and alloys. Some analytic solutions can be derived). They are very intelligent but lacks the atomic detail. An alterative approach is directly model the system (straightforward but very time consuming and tedious). Notable examples include. 
1. [Buffon's needle experiment to compute $\pi$](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem), 
2. [Bernal's ball-bearing model](https://iopscience.iop.org/article/10.1088/0953-8984/26/46/463102), 
3. [Kitaigorodskii’s structure seeker](https://pubs.acs.org/doi/10.1021/acs.cgd.8b00972).

* Molecular dynamics Simulation is generally a technical to directly study the atomic evolution of atoms or molecules in the material system based on the simple law of Newtonian Dynamics. Immediately after the invention of computers, MD simulations have been quickly applied to several systems
1. First ever, 1956, Alder and Wainwright (Livermore), [Phase transition of hard spheres](https://gibbs.ccny.cuny.edu/teaching/s2021/labs/HardDiskSimulation/Alders&Wainwright1957.pdf)
2. Real material, 1959, Gibson, Goland, Milgram, Vineyard (Brookhaven), [Dynamics of radiation damage](https://journals.aps.org/pr/abstract/10.1103/PhysRev.120.1229)
3. Liquid, 1964, Rahman, [Correlations in Liquid Argon](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405)

### 1.1.2 Why do we need such kind of direct simulation method like MD
* Go beyond the experimental capability
* Gain some physical insights

### 1.1.3 Homework (W1M1): 
1. Read the Alder & Wainwright's [1956 paper]([Phase transition of hard spheres](https://gibbs.ccny.cuny.edu/teaching/s2021/labs/HardDiskSimulation/Alders&Wainwright1957.pdf) and understand [hard sphere potential](https://en.wikipedia.org/wiki/Hard_spheres)
2. Read the [Rahman](https://en.wikipedia.org/wiki/Aneesur_Rahman)'s [1964 paper](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405). We will try to reproduce some results from this work in this week. 

## 1.2 A first MD simulation of liquid argon under NVE ensemble.

### 1.2.1 A Simple MD workflow

```pseudo
I. Initialization
    * Set simulation parameters (e.g., number of particles, temperature, time step)
    1. particle positions (R) randomly or on a lattice
    2 particle velocities (V) Maxwell-Boltzmann distribution

II. Computation of Energy and Forces
    * Compute potential energy using the chosen potential
    * Compute forces on each particle by differentiating the potential energy

III. Integration (Time Evolution)
    * For each time step:
        1. Update particle positions using the integration algorithm (e.g., Verlet)
        2. Update particle velocities based on the computed forces
        3. Apply periodic boundary conditions (if necessary)
        4. Recompute forces on all particles
        5. Update potential energy

IV. Termination
    * Check if the simulation time has reached the desired number of steps
```
### 1.2.2 Interatomic potential

[Lennard Jones Potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) express the assumes that all particles interact with each other via pairwise interactions (i.e., the interaction energy only depends on the distance).

$$
V(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]
$$

It consist of two components:
* $1/r^{12}$ term to account for short range repulsion
* $-1/r^6$ term to account for the long range attraction (also called [London dispersion force](https://en.wikipedia.org/wiki/London_dispersion_force))
* ε and σ control the equilibrum distance and shape of the energy well. 

This form has been widely used to model the essential features of interactions between simple atoms and molecules.

#### Questions:
1. Why were the 12 and 6-powered terms were chosen? Any other choices?
2. How does the LJ potential decay with respect to r?
3. The limitations of LJ potential?
4. Can we use them to model metal, ceramics or others?

### 1.2.3 The computation of energy and forces

After knowing the energy model, we can proceed to compute the total energy and forces for each particle in the given system.
Assuming the system consists of $N$ atoms, and the positions ($R$) are recorded by an array of [N, 3], we can use the following pseudo Python code for the given task.
```python
import numpy as np

E = 0                    # Total Energy
F = np.zeros([N, 3])     # Atomic forces

for atom i in range(N-1):
    for atom j in range(i+1, N):
        Compute the distance: R = Ri - Rj
        Compute the energy: V = V(R_ij)
        E += V
        Compute the derivative w.r.t R: dE/dR_i and dE/dR_j
        F(i) -= dE/dR_i
        F(j) -= dE/dR_j
```

### 1.2.4 Notes: derivation of Forces
The force  $\mathbf{F}(r)$ between two particles is the negative gradient of the potential:

$$
\mathbf{F}(r) = -\frac{dV(r)}{dr}
$$

Taking the derivative of the potential with respect to  r :

$$
\begin{align*}
\frac{dV(r)}{dr} &= 4\epsilon \left[ -12\left(\frac{\sigma}{r}\right)^{12} \frac{1}{r} + 6\left(\frac{\sigma}{r}\right)^{6} \frac{1}{r} \right]\\
                 &= \epsilon \left[ \frac{48\sigma^{12}}{r^{13}} - \frac{24\sigma^{6}}{r^{7}} \right]
\end{align*}
$$

In practice, the $r$ is a 3-vector $(x, y, z)$, to compute the force component on like $F_x(r)$, there is an addition work.
The distance $r$ between two particles can be written as:

$$
r = \sqrt{x^2 + y^2 + z^2}
$$

$$
\begin{align*}
F_x(r) &= -\frac{dV(x, y, z)}{dx} = -\frac{dV(x, y, z)}{dr}\frac{dr}{dx} \\
       &= - \frac{dV(x, y, z)}{dr} \frac{x}{r} \\
       &= x \epsilon \left[ \frac{48\sigma^{12}}{r^{14}} - \frac{24\sigma^{6}}{r^{8}} \right]
\end{align*}
$$


### 1.2.5 Initialization
* Atoimic positions
If we study a system that mimics solid, we can just create the position on a lattice compatible with the structure that we aim to simulate. You must avoid the case of geometry consisting of two atoms with very short distances.

* Velocities
Ideally, we should generate the initial velocities to follow the [Maxwell-Boltzmann distribution](https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution).

$$
p(v) = 4\pi \left( \frac{m}{2\pi k_B T} \right)^{3/2} v^2 \exp\left(-\frac{mv^2}{2k_B T}\right)
$$


To achieve this, the idea is to sample velocities from a normal (Gaussian) distribution, where the standard deviation is related to the temperature and the mass of the particles.

```python
import numpy as np

def generate_velocities(num_particles, temperature, mass):
    # Boltzmann constant in appropriate units (e.g., J/K or eV/K)
    k_B = 1.380649e-23  # J/K

    # Standard deviation of the velocity distribution
    sigma_v = np.sqrt(k_B * temperature / mass)

    # Generate velocities from a normal distribution
    velocities = np.random.normal(0, sigma_v, (num_particles, 3))

    return velocities

# Example usage
num_particles = 1000       # Number of particles
temperature = 300          # Temperature in Kelvin
mass = 1.67e-27            # Mass of a particle in kg (e.g., proton)

velocities = generate_velocities(num_particles, temperature, mass)

# Subtract the mean velocity to ensure zero net momentum
mean_velocity = np.mean(velocities, axis=0)
velocities -= mean_velocity

print("Velocities adjusted for zero total momentum (m/s):")
print(velocities[:5])
```
<!---In addition the above script, one can also use [scipy.stats.maxwell](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.maxwell.html#scipy.stats.maxwell) to generate random samples.

Note, for the current task of NVE simulation, we would compute velocity from $\mathbf{r}(t)$ trajectory. Therefore, this setup would only impact the first step of integration. If the simulation converges, it won't have impacts on the simulation.
--->

### 1.2.6 Integrator (updating rule)
After knowing the forces, we can proceed to update the velocities (V) and positions (R) for the next time step:

$$
\begin{align*}
\mathbf{r}(t+dt) &= \mathbf{r}(t) + \mathbf{v}(t)dt\\
\mathbf{v}(t+dt) &= \mathbf{v}(t) + \mathbf{a}(t)dt
\end{align*}
$$

However, this update will suffer from a rapid propagation of error. To reduce the error propogation, we use the so called [Velocity Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration). 

$$
\begin{align*}
\mathbf{r}(t + dt) &= \mathbf{r}(t) + \mathbf{v}(t)dt + 0.5 \mathbf{a}(t) dt^2  \\
\mathbf{v}(t + dt) &= \mathbf{v}(t) + 0.5 [\mathbf{a}(t) + \mathbf{a}(t+dt)]dt
\end{align*}
$$

According to Taylor expansion, this algorithm is accurate to $O(dt^3)$ in position and $O(dt^2)$ in velocity. 

## 1.3 Full code
Complete the codes in [Colab](https://colab.research.google.com/drive/1lB3R0N_s2gP-IhjrxBWq2mDW2VlqIE_c#scrollTo=KDtmzZIA2kvp)

### 1.3.1 Summary of Code Implementation
1. Make sure you have go over all equations and finish the pseudo code before writing the real code
2. Split the entire workflow into several subtasks.
3. For each subtask, make sure you have some mechanisms to validate and debug your code.
4. Validate the final results with some physical guidance (in NVE simulation, ensure you check if the total energy is conserved).

Hopefully, you are able to write a basic code for MD simulation after this practice. You are expected to reinforce your understanding by writing your own code. 

### 1.3.2 Cross-validation with other MD codes.
Of course, there are many excellent open-source MD codes with more functional support. For productive research project, you would probably use those codes. In this course, we recommend the use of [LAMMPS](https://github.com/lammps/lammps), one of the most popular code for materials modelling.

For students who already have LAMMPS experience, there is a bonus credit opportunity. 
Please rerun the simulation in LAMMPS with the same parameter setup. Post your LAMMPS script to [our forum](https://github.com/qzhu2017/AtomisticSimulation/issues/1)


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

1. With a certain probability  $\nu$  (collision frequency), the velocity of each particle is reassigned by sampling from a Maxwell-Boltzmann distribution corresponding to the desired temperature T.
1. If the particle’s velocity is not reassigned, it evolves according to the usual equations of motion (e.g., using the Verlet integration method).

In this technique, collision Frequency $$\nu$$ determines how often the particle velocities are randomized (following a [Poisson Distribution](https://en.wikipedia.org/wiki/Poisson_distribution)). A higher $\nu$ means more frequent collisions (interaction) with the heat bath, leading to stronger coupling to the temperature bath. We should choose a $$\nu$$ so that velocity reassignment happens at an appropriate rate to maintain the desired temperature without overly disrupting the natural dynamics of the system.

The Anderson thermostat is relatively simple to implement, requiring only the addition of random velocity reassignment at each time step. However, it may not reflect the real dynamics. Since velocities are randomly reassigned, the resulting particle trajectories may not correspond to those in a real physical system where energy exchange occurs through physical interactions. This is particularly true for a periodic  system without the clear definition of boundary. In addition, one needs to play with the $\nu$ values.


The algorithm can be described as follows.

```python

def anderson_thermostat(V, T, mass, nu, dt):
    sigma_v = np.sqrt(k_B * temperature / mass)
    for i in range(num_particles):
        if np.random.rand() < collision_frequency * time_step:
            V[i] = np.random.normal(0, sigma_v, 3)
    return V

# Initialization
R = Initialize_positions()
V = Initialize_velocities()
F = calculate_forces(R)

# Main MD loop
for step in range(num_steps):

    # Update R, V, F using Verlet integration
    R += V * dt + 0.5 * F * dt**2 / mass
    F_new = calculate_forces(positions)
    V += 0.5 * (F + F_new) * dt / mass
    F = F_new

    # Apply Anderson thermostat to randomly reassign V
    V = anderson_thermostat(V, T, mass, nu, dt)
```

### 2.2.2 The Langevin thermostat 
The Langevin thermostat maintains the temperature of a system while also modeling the effects of friction and random forces, similar to those that might be encountered in a viscous fluid. The basic idea it to modify the equations of motion by adding two additional terms to the standard Newtonian dynamics:

$$
m_i \frac{d\mathbf{v}_i}{dt} = \mathbf{F}_i - \gamma m_i \mathbf{v}_i + \mathbf{R}_i(t)
$$

* Frictional Force  $\gamma m_i \mathbf{v}_i$  : represents the damping effect of the environment, which tends to slow down the particles.

* Random Force  $\mathbf{R}_i(t)$ : represents the random collisions with particles from the heat bath, which cause the particles to move randomly, maintaining the system’s temperature. These kicks help maintain the temperature of the system by continuously injecting energy into the system.

A typical value of $\gamma$ used in many MD simulations is around $\gamma = 0.1 \, \text{ps}^{-1}$. This value provides a good balance between maintaining temperature control and preserving realistic dynamics. The system is weakly coupled to the heat bath, ensuring that it can sample the canonical ensemble without heavily damping the natural motion of the particles.

In MD simulations, the Langevin equation is integrated with the Velocity Verlet algorithm, modified to include the Langevin terms. A simple update rule for the velocity might look like this:

$$
\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i(t) + \frac{\mathbf{F}_i(t)}{m_i} \Delta t - \gamma \mathbf{v}_i(t) \Delta t + \mathbf{R}_i \sqrt{\Delta t}
$$

```python

import numpy as np


def langevin_thermostat(V, F, gamma, sigma, dt):
    # Update velocities with deterministic part
    V += 0.5 * F * dt / mass

    # Apply friction and random force (stochastic part)
    V += -gamma * V * dt
    V += np.random.normal(0, sigma, V.shape) * np.sqrt(dt)

    return V

# Initialization
R = Initialize_positions()
V = Initialize_velocities()
F = calculate_forces(R)

# Main MD loop
for step in range(num_steps):

    # Update R, V, F using Verlet integration
    R += V * dt + 0.5 * F * dt**2 / mass

    # Apply Langevin thermostat to update velocities
    V = langevin_thermostat(V, F, gamma, sigma, dt)

	# Update F and velocity using Verlet
    F_new = calculate_forces(positions)
    V += 0.5 * (F + F_new) * dt / mass
    F = F_new
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
    R += velocities * time_step + 0.5 * forces * time_step**2 / mass
    F_new = calculate_forces(R)

    # Update velocities
    V += 0.5 * (F + F_new) * dt / mass
    V *= (1 - 0.5 * xi * dt) / (1 + 0.5 * xi * dt)

    # Update the Nosé-Hoover thermostat variable
    kE= 0.5 * np.sum(mass * V**2)
    xi += dt * (2 * kE / (3 * N * k_B * T) - 1) / Q
```

### 2.2.4 Extended Discussions in Ensembles and Nosé-Hoover thermostat
Complete the reading in [Appendix-W2](https://github.com/qzhu2017/AtomisticSimulation/blob/main/Appendix/W2_NoseHoover.pdf).

## 2.3 Full code to run NVT simulation
Complete the codes in [Colab](https://colab.research.google.com/drive/1x8FFEDrvmThUQGhfVCJd9LZka0aO1zhe?usp=sharing)

### 2.3.1 Summary of Code Implementation
1. Make sure you have go over all equations and finish the pseudo code before writing the real code
2. Split the entire workflow into several subtasks.
3. Validate the final results with some physical guidance (in NVT simulation, ensure you check if the final kinetic energy fluctuate around the desired value).

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
P = \frac{N k_B T}{V} + \frac{1}{3V} \sum_{i < j} r_{ij} \cdot \mathbf{F}_{ij}
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


# Week 4: Realistic MD simulation with LAMMPS

## 4.1 Introduction to LAMMPS
LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is an open-source software tool designed for performing classical molecular dynamics (MD) simulations. It can be used to model an array of particle interactions, ranging from simple atomic systems to complex materials and biomolecular systems. As one of the most population materials simulation package, LAMMPS is specifically optimized for large-scale simulations, which involve millions to billions of particles, making it suitable for high-performance computing (HPC) environments. Its versatility allows for simulations of a variety of systems such as metals, polymers, proteins, granular media, and even non-equilibrium phenomena like crack propagation, heat transfer, and more. Some typical applications include

- Crystal Defects and Deformation: Investigating dislocation motion, grain boundary evolution, and fracture in materials.
- Thermal Conductivity: Computing thermal transport properties of materials via Green-Kubo or direct methods.
- Phase Transitions: Simulating phase changes in metals and alloys, such as melting, solidification, or the formation of microstructures.
- Mechanical Properties: Studying stress-strain relationships, elasticity, and plasticity at the atomic scale.
- Protein Folding: Understanding the folding and dynamics of proteins or peptides.
- Drug Interactions: Simulating how drug molecules interact with proteins or other biological targets.
- Membrane Simulations: Investigating lipid bilayers, ion channels, or membrane proteins.

## 4.2 Why is LAMMPS efficient and popular?
LAMMPS is designed to perform large-scale simulations by taking advantage of parallel computing architectures, including multi-core processors and high-performance computing (HPC) clusters. In particular, LAMMPS domain decomposition technique to divide the simulation space into subdomains. Each processor or core is assigned a subdomain, and they work together by communicating boundary conditions and interacting forces. LAMMPS also uses MPI to handle communication between processors, ensuring minimal overhead and efficient data transfer during the simulation. Thanks to these considerations, LAMMPS has demonstrated excellent scalability across thousands of processors, which makes it suitable for simulating systems with millions to billions of particles. This scalability is crucial for studying complex, large-scale systems like polymers, metals, or biomolecules over long time scales.

- Support for Multiple Interatomic Potentials (LJ, Embedded Atom Method, .etc)
- Flexible Input Script: LAMMPS uses a simple, flexible scripting language for setting up simulations. This language allows users to define new potential types, set conditions, and run specific simulation protocols.

## 4.3 Input and Output in LAMMPS

After you compiled the LAMMPS code into an executable (often called `lmp_serial` or `lmp_mpi` by default). 

```
path_to_lmp_mpi < lmp.in > lmp.out
```

### LAMMPS Input
LAMMPS simulations are controlled by input scripts, which consist of a series of commands written in a simple text format. These scripts define the simulation parameters, system setup, and specific instructions for running the simulation.

Structure of LAMMPS Input Scripts

A typical LAMMPS input script is organized into several key sections:
```
- Initialization Section (Global settings, such as the units, boundary conditions, and atom styles).
units real: Defines the units for physical quantities (e.g., distance in angstroms, energy in kcal/mol).
boundary p p p: Sets periodic boundary conditions in all three spatial dimensions.
atom_style atomic: Defines how atoms are represented (e.g., atomic, charge, molecular).

- Atom Definition Section (Atomic coordinates, and initial velocities).
read_data data.file: Reads the atomic configuration and other system properties from a file .
velocity all create 300.0 12345: Assigns random initial velocities to atoms at a temperature of 300 K with a random seed.

- Force Field Definition Section:
pair_style lj/cut 2.5: Defines a Lennard-Jones potential with a cutoff distance of 2.5.
pair_coeff * * 0.1 3.0: Sets the Lennard-Jones coefficients (epsilon and sigma) for the interactions between atom types.

- Simulation Parameters Section:

timestep 1.0: Sets the time step for integration to 1.0 (in the units defined by the units command).
fix 1 all nve: Applies a constant energy (NVE) ensemble to all atoms.
fix 2 all temp/rescale 100 300.0 300.0 0.02 1.0: Rescales the temperature every 100 timesteps to maintain a temperature of 300 K.

- Output Control Section (frequency and format of the simulation output).
thermo 100: Prints thermodynamic data (e.g., temperature, pressure, energy) every 100 timesteps.
dump 1 all atom 1000 dump.atom: Outputs the atomic positions to a file (dump.atom) every 1000 timesteps.

- Run Section:
run 10000: Runs the simulation for 10,000 timesteps.
minimize 1.0e-4 1.0e-6 1000 10000: Performs energy minimization on the system with specified tolerance and iteration limits.
```

### LAMMPS Output

- `lammps.log`. The log file is generated automatically for every LAMMPS run and records all the commands executed in the input script. It also contains thermodynamic data (e.g., energy, pressure, temperature) at intervals specified by the thermo command.
```
Step Temp E_pair E_mol TotEng Press Volume
0 300.0 -143.53 0 -123.46 1.35 1000.0
100 310.2 -142.11 0 -122.45 1.40 1000.0
```

- Dump Files: Dump files store detailed trajectory information about the system’s atomic coordinates, velocities, and forces. These files are typically used for post-processing to analyze system configurations, create visualizations, or calculate structural properties.

```
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
1000
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
3 2 0.5 0.5 0.5
```

- Restart Files: store the entire state of a LAMMPS simulation, allowing users to pause and later continue a simulation from where it left off. These files contain information about the atom positions, velocities, forces, and other system properties.


## 4.3 Simulation Process

LAMMPS begins with setting up the system and defining parameters. This involves specifying the atomic system geometry, simulation box, interatomic potentials, and initial conditions like atomic velocities and temperature.

```
units real
atom_style atomic
read_data argon.data
pair_style lj/cut 2.5
pair_coeff * * 0.238 3.4   # Argon-specific parameters for Lennard-Jones potential
velocity all create 300.0 12345

fix 1 all nvt temp 300.0 300.0 100.0   # NVT ensemble to control temperature
timestep 1.0                           # Time step of 1.0 fs
run 50000                              # Run simulation for 50,000 timesteps
```

Post-Processing

After a simulation, the results need to be visualized and analyzed. LAMMPS produces several types of output files, which contain thermodynamic data, atom positions, velocities, and forces. VMD (Visual Molecular Dynamics) and OVITO (Open Visualization Tool) are popular tools for visualizing molecular dynamics simulations.

$$ 4.4 Running LAMMPS on HPC
While LAMMPS can be ran at many different platforms. For most research projects, running LAMMPS on a high-performance computing (HPC) environment is essential for large-scale simulations that require significant computational resources. Most modern supercomputers use job schedulers like SLURM to manage computational tasks.

```
#!/bin/bash
#SBATCH --job-name=lammps_job          # Job name
#SBATCH --nodes=4                      # Number of nodes
#SBATCH --ntasks-per-node=32           # Number of tasks (processes) per node
#SBATCH --time=24:00:00                # Max time limit (HH:MM:SS)
#SBATCH --partition=compute            # Partition or queue to submit to
#SBATCH --output=job_output.log        # Output log file

module load lammps/3Mar2020            # Load LAMMPS module
mpirun -np 128 lmp_mpi -in input_file.in   # Run LAMMPS in parallel across 128 processes
```

SBATCH options are used to specify the number of nodes, tasks, job name, and time limit.

- mpirun -np 128 launches LAMMPS across 128 processes in parallel, ensuring that the simulation scales across multiple cores.
- lmp_mpi is the parallel version of LAMMPS used for multi-node execution.




