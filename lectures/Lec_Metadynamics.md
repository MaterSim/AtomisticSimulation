# Week 8: Enhanced Sampling with Metadynamics

In standard MD, the system evolves under Newtonian dynamics, meaning that it can get stuck in local minima for long periods. 
This limits the exploration of the free energy surface, making it difficult to observe transitions between different states. 
To address this challenge, one may consider the use of enhanced sampling techniques. In this lecture, we will focus on a particular type 
of technique called Metadynamics.


## 8.1 What is Metadynamics?

Metadynamics is a MD-based sampling technique to explore free energy landscapes and overcome energy barriers. 
It is particularly useful for systems where traditional MD struggles to sample rare events due to high energy barriers or complex transitions. 
not only accelerates the sampling of rare events but also provides a mechanism for reconstructing the free energy surface. 
This is achieved by systematically filling the energy wells with biasing potentials, eventually allowing the system to overcome energy barriers 
and explore alternative configurations. 

Compared to the traditional MD simulation, Metadynamics adds the following concepts:

1. **Biasing Potential**: The fundamental idea is to add a bias potential to the system that changes over time. 
By adding Gaussian-shaped potentials periodically to the explored regions, the system is pushed out of energy wells, allowing for better exploration. 
The accumulation of these Gaussians eventually results in a bias that compensates for the underlying free energy, 
enabling the system to explore new configurations freely.

2. **Collective Variables (CVs)**: Metadynamics works by adding the bias in the space of carefully chosen collective variables (CVs). 
CVs are reduced-dimensional representations of the system, such as distances between atoms, angles, or other descriptors capturing essential behavior. 
The choice of CVs is crucial, as it directly influences the efficiency and success of the metadynamics simulation. 
Properly chosen CVs can significantly accelerate the exploration of relevant conformational space, while poorly chosen CVs may lead to incomplete sampling.

3. **Free Energy Estimation***: As the bias potential accumulates, it fills the wells of the free energy landscape. 
The accumulated bias potential can be used to reconstruct the underlying free energy surface, providing valuable insights into the 
thermodynamic properties of the system. This reconstruction allows researchers to identify stable states, transition states, 
and the pathways connecting them, which is critical for understanding the behavior of molecular systems.

## 8.2 MD and Metadynamics in 1D potential well

To demonstrate the concept of Metadynamics, let's consider simple model system with a particle moving in a double-well potential as follows

$$
V(x) = x^4 - 2x^2.
$$

```python
import numpy as np
import matplotlib.pyplot as plt

# Define the double well potential
def double_well_potential(x):
    return x**4 - 2*x**2

# Derivative of the potential (force)
def potential_force(x):
    return -4*x**3 + 4*x

# Time evolution of the particle using Langevin dynamics 
def md(steps=5000, dt=0.01, gamma=0.1, temp=0.1):
    x = 0.0  # initial position
    x_positions = [x]

    for step in range(steps):
        # Langevin dynamics with force from the double well potential
        force = potential_force(x)
        thermal_force = np.sqrt(2 * gamma * temp / dt) * np.random.normal()

        # Update position with Langevin equation
        x += force * dt - gamma * x * dt + thermal_force * np.sqrt(dt)
        x_positions.append(x)

    return np.array(x_positions), centers

# Run the simulation
x_positions = md()

# Plot the results
x_vals = np.linspace(-2, 2, 100)
potential_vals = double_well_potential(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, potential_vals, label="Double Well Potential", color="blue")
plt.xlabel('Position (x)')
plt.ylabel('Potential Energy')
plt.legend()
plt.title('Double Well Potential with Langevin dynamics')
plt.grid(True)
plt.show()
```

Clearly, one can find that the particle would just oscillate around the energy well in this simulation. 

If one is interested in sampling ***. 

$$
V_{\text{biased}}(s,t) = V_{\text{system}}(s) + \sum_{t{\prime} \leq t} W \exp\left( -\frac{(s - s(t{\prime}))^2}{2 \sigma^2} \right)
$$

Metadynamics adds a bias to the system to prevent the particle from getting stuck in one well, 
allowing it to explore the other regions of the potential landscape.

```python
def gaussian_bias(x, centers, width=0.1, height=0.1):
    bias = 0
    for c in centers:
        bias += height * np.exp(-0.5 * (x - c)**2 / width**2)
    return bias

# Time evolution of the particle using Langevin dynamics with Metadynamics
def metadynamics_simulation(steps=5000, dt=0.01, gamma=0.1, temp=0.1):
    x = 0.0  # initial position
    x_positions = [x]
    centers = []  # store the positions where bias is added

    for step in range(steps):
        # Langevin dynamics with force from the double well potential
        force = potential_force(x)
        thermal_force = np.sqrt(2 * gamma * temp / dt) * np.random.normal()

        # Apply bias from metadynamics
        bias_force = -np.gradient(gaussian_bias(x, centers))

        # Update position with Langevin equation
        x += (force + bias_force) * dt - gamma * x * dt + thermal_force * np.sqrt(dt)
        x_positions.append(x)

        # Add Gaussian bias every 100 steps
        if step % 100 == 0:
            centers.append(x)

    return np.array(x_positions), centers

# Run the simulation
x_positions, centers = metadynamics_simulation()

# Plot the results
x_vals = np.linspace(-2, 2, 100)
potential_vals = double_well_potential(x_vals)
bias_vals = [gaussian_bias(x, centers) for x in x_vals]

plt.figure(figsize=(10, 6))
plt.plot(x_vals, potential_vals, label="Double Well Potential", color="blue")
plt.plot(x_vals, potential_vals + bias_vals, label="Metadynamics Bias", color="orange")
plt.xlabel('Position (x)')
plt.ylabel('Potential Energy')
plt.legend()
plt.title('Double Well Potential with Metadynamics Bias')
plt.grid(True)
plt.show()
```

This simple illstrate the main idea behind metadynamics. 



## 8.3 How Metadynamics Works?

In practice, metadynamics simulations can be applied to higher-dimensional systems and more complex potentials 
such as chemical reactions, protein folding and phase transitions. 
It typically involves the following steps:

1. Initialization: Select appropriate CVs that describe the transition of interest. 
The choice of CVs is critical as they determine how effectively metadynamics can enhance sampling. 
CVs should capture the essential degrees of freedom involved in the transition, and their proper selection often requires prior knowledge of the system.

2. Bias Addition: During the simulation, small Gaussian potentials are periodically added along the CVs to the current position of the system. 
This gradually discourages the system from revisiting already visited states. 
The height and width of these Gaussians are important parameters that control the rate of exploration; 
a careful balance must be struck to ensure efficient sampling without overshooting important regions of the free energy landscape.

3. Exploration of the Free Energy Surface: By continuously adding Gaussians, the system is encouraged to move out of energy minima, 
eventually allowing the exploration of the entire relevant free energy landscape. The bias potential effectively smooths out the energy barriers, 
enabling the system to transition between states more easily. Over time, the system visits all accessible regions of the free energy surface, 
and the accumulated bias provides an estimate of the free energy differences between states.

The bias potential is given by:

Where:

 is the height of the Gaussian.

 is the width of the Gaussian.

 represents the CVs at time .

Gaussian Parameters: The height and width of the Gaussians must be carefully chosen to balance between exploration speed and accuracy. 
If the Gaussians are too high or too wide, the system may be forced to explore unphysical regions, 
whereas too small Gaussians may lead to slow convergence.

## 8.4 Choice of CVs

The efficiency of metadynamics heavily depends on the proper choice of CVs. Poor selection can lead to ineffective sampling, 
as the system may not be driven along the most relevant pathways. 
The process of selecting suitable CVs often requires trial and error or prior knowledge of the system.


## 8.5 Well-Tempered Metadynamics

In traditional metadynamics, the constant addition of Gaussian potentials can lead to excessive bias accumulation, 
which may result in poor sampling or an inaccurate free energy landscape. 
Well-tempered metadynamics addresses this by gradually reducing the height of the Gaussians added over time, 
which helps prevent oversampling and ensures that the system does not accumulate too much bias in any particular region.
 As such, it is expected to improve convergence and enhance the accuracy of the free energy surface estimation. 


The key idea behind well-tempered metadynamics is to scale the bias deposition rate according to the amount of bias already present. 
This scaling is achieved by introducing a parameter called the bias factor , 
which controls how much the bias potential decreases as the simulation progresses. 
The bias factor is related to a fictitious temperature that effectively dictates how smoothly the bias is added. 
The bias potential in well-tempered metadynamics evolves as:

Where:

 is the Boltzmann constant.

 is the bias factor that determines the rate of decrease in bias deposition.

The bias factor  effectively controls the level of exploration versus exploitation in the simulation. 
A larger bias factor results in slower reduction of the bias potential, allowing the system to continue exploring new regions of the free energy surface. 
Conversely, a smaller bias factor leads to faster convergence, which is beneficial for accurately reconstructing the free energy surface without excessive biasing.

Well-tempered metadynamics has several advantages over traditional metadynamics. By reducing the rate of bias deposition over time, 
it ensures that the system can focus more on the most relevant regions of the free energy surface, leading to a more accurate reconstruction. 
Additionally, well-tempered metadynamics helps maintain a balance between exploration of new configurations and refinement of previously visited regions, 
which is crucial for systems with multiple metastable states or complex free energy landscapes.

Another significant benefit of well-tempered metadynamics is that it provides a natural mechanism for achieving convergence of the free energy surface. 
As the system becomes more thoroughly explored, the bias added to the system decreases, ultimately reaching a point where it no longer significantly alters the free energy landscape. This gradual reduction in bias allows the system to settle into the correct free energy minima, providing a reliable estimate of the underlying free energy differences between states.

Well-tempered metadynamics has been successfully applied to a wide range of systems, from simple model potentials to complex biomolecular and materials 
science applications. Its ability to adaptively control the bias potential makes it a versatile tool for studying processes such as protein folding, 
chemical reactions, and phase transitions. By providing a more controlled and convergent approach to free energy estimation, 
well-tempered metadynamics has become a preferred method for enhanced sampling in many challenging molecular simulations.

Moreover, well-tempered metadynamics can be used in conjunction with other enhanced sampling methods, 
further improving the efficiency and accuracy of molecular simulations. For instance, combining well-tempered metadynamics 
with umbrella sampling or replica exchange can provide additional sampling power and improve the characterization of complex energy landscapes. 
This versatility makes well-tempered metadynamics a valuable tool in the toolbox of computational chemists and molecular modelers.

Well-tempered metadynamics also benefits from relatively simple implementation, 
which makes it accessible for integration into various molecular dynamics software packages. 
Popular MD software such as GROMACS, LAMMPS, and PLUMED offer built-in support for well-tempered metadynamics, 
allowing researchers to easily incorporate this technique into their workflows. 
This accessibility has contributed to the widespread adoption of well-tempered metadynamics in both academic research and industrial applications, 
where the ability to accurately estimate free energy surfaces is of great interest.


Metadynamics is a powerful enhanced sampling technique to overcome the limitations of traditional MD simulation. 
Its ability to reconstruct free energy surfaces makes it an invaluable tool for studying complex molecular systems, phase transitions, and reaction mechanisms. 
However, its success depends on careful selection of collective variables and parameters. 
By systematically adding bias potentials, metadynamics allows for the study of rare events and provides detailed insights into the thermodynamics and 
kinetics of molecular systems. Its versatility and adaptability make it a widely used technique in various fields of molecular simulation, 
from chemistry to materials science.

References

Laio, A., & Parrinello, M. (2002). Escaping free-energy minima. Proceedings of the National Academy of Sciences, 99(20), 12562-12566.

Barducci, A., Bussi, G., & Parrinello, M. (2008). Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method. Physical Review Letters, 100(2), 020603.