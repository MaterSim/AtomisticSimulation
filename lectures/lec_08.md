# 8 Enhanced Sampling with Metadynamics
In standard MD, the system evolves under Newtonian dynamics, meaning that it can get stuck in local minima for long periods. 
This limits the exploration of the free energy surface, making it difficult to observe transitions between different states. 
To address this challenge, one may consider the use of enhanced sampling techniques. In this lecture, we will focus on a particular type 
of technique called Metadynamics.

## 8.1 What is Metadynamics?

Metadynamics is an MD-based sampling technique to explore free energy landscapes and overcome energy barriers. 
It is particularly useful for systems where traditional MD struggles to sample rare events due to high energy barriers or complex transitions. 
Metadynamics not only accelerates the sampling of rare events but also provides a mechanism for reconstructing the free energy surface. 
This is achieved by systematically filling the energy wells with biasing potentials, eventually allowing the system to overcome energy barriers 
and explore alternative configurations. 

Compared to the traditional MD simulation, Metadynamics adds the following concepts:

1. **Biasing Potential**: The fundamental idea is to add a bias potential to the system that changes over time. By adding Gaussian-shaped potentials periodically to the explored regions, the system is pushed out of energy wells, allowing for better exploration. The accumulation of these Gaussians eventually results in a bias that compensates for the underlying free energy, enabling the system to explore new configurations freely.

2. **Collective Variables (CVs)**: Metadynamics works by adding the bias in the space of carefully chosen collective variables (CVs). CVs are reduced-dimensional representations of the system, such as distances between atoms, angles, or other descriptors capturing essential behavior. The choice of CVs is crucial, as it directly influences the efficiency and success of the metadynamics simulation. Properly chosen CVs can significantly accelerate the exploration of relevant conformational space, while poorly chosen CVs may lead to incomplete sampling.

3. **Free Energy Estimation**: As the bias potential accumulates, it fills the wells of the free energy landscape. The accumulated bias potential can be used to reconstruct the underlying free energy surface, providing valuable insights into the thermodynamic properties of the system. This reconstruction allows researchers to identify stable states, transition states, and the pathways connecting them, which is critical for understanding the behavior of molecular systems.

## 8.2 MD and Metadynamics in 1D potential well

To demonstrate the concept of Metadynamics, let's consider simple model system with a particle moving in a double-well potential as follows

$$
V(x) = x^4 - 3x^2
$$

```python
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(font_scale=1.2)

# Define the double well potential
def double_well_potential(x):
    return x**4 - 3*x**2

# Derivative of the potential (force)
def potential_force(x):
    return -4*x**3 + 6*x

# Time evolution of the particle using Langevin dynamics
def md(steps=20000, dt=1e-2, gamma=0.02, temp=0.01):
    x = 0.5  # initial position
    x_positions = [x]

    for step in range(steps):
        # Langevin dynamics with force from the double well potential
        force = potential_force(x)
        thermal_force = np.sqrt(2 * gamma * temp / dt) * np.random.normal()

        # Update position with Langevin equation
        x += force * dt - gamma * x * dt + thermal_force * np.sqrt(dt)
        x_positions.append(x)

    return np.array(x_positions)

# Run the simulation
xs = md()

# Create the figure with two subplots with shared x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

# Double well potential with Langevin dynamics
x_vals = np.linspace(-2, 2, 100)
potential_vals = double_well_potential(x_vals)
ax1.plot(x_vals, potential_vals, '--', lw=1.0, label="Double Well Potential", color="k")
sc = ax1.scatter(xs, double_well_potential(xs), c=range(len(xs)), s=2, cmap='viridis', alpha=0.5)
ax1.set_ylabel('Potential Energy')
ax1.legend()
ax1.set_title('Double Well Potential with Langevin Dynamics', fontweight='bold')
ax1.grid(True)
cbar = plt.colorbar(sc, ax=ax1, orientation='horizontal', pad=0.1)
cbar.set_label('Time Step')

# Histogram of x values
ax2.hist(xs, bins=50, color='skyblue', edgecolor='black')
ax2.set_xlabel('Position (x)')
ax2.set_ylabel('Frequency')
ax2.set_title('Histogram of Particle Positions')
ax2.grid(True)

plt.tight_layout()
plt.savefig('1D-MD.png')
```

Clearly, one can find that the particle would just oscillate around the energy well in this simulation. 
   <img src="https://github.com/qzhu2017/AtomisticSimulation/blob/main/Codes/lec_08_MD.png" alt="MD" height="400">

If one is interested in sampling more of the phase space, there must be a way to escape from the local minima. Metadynamics adds a bias to the system to prevent the particle from getting stuck in one well, allowing it to explore other regions of the potential landscape.


$$
V_{\text{bias}}(s,t) = V_{\text{system}}(s) + \sum_{t{\prime} \leq t} W \exp\left( -\frac{(s - s(t{\prime}))^2}{2 \sigma^2} \right)
$$

Where:

- $W$ is the height of the Gaussian, which controls the magnitude of the bias added at each time step.
- $\sigma$ is the width of the Gaussian, determining how localized the bias is in the CV space.
- $s(t)$ represents the value of the CVs at time $t$.

The bias potential $V_{\text{bias}}(s,t)$ accumulates Gaussians placed at the positions the system has visited in the CV space, thus gradually filling the wells in the free energy landscape and pushing the system to explore other areas.

Below shows how to implement the bias potential in the context of metadynamics.

```python
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(font_scale=1.2)

# Derivative of the potential (force)
def potential_force(x):
    return -4*x**3 + 6*x

def gaussian_bias(x, centers, width=0.1, height=0.1):
    bias = 0
    for c in centers:
        bias += height * np.exp(-0.5 * (x - c)**2 / width**2)
    return bias

# Derivative of the Gaussian bias (force due to bias)
def bias_force(x, centers, width=0.1, height=0.1):
    force = 0
    for c in centers:
        force += height * (x - c) / (width**2) * np.exp(-0.5 * (x - c)**2 / width**2)
    return force

# Time evolution of the particle using Langevin dynamics
def metaD(steps=20000, dt=1e-2, gamma=0.02, temp=0.01):
    x = 0.5  # initial position
    x_positions = [x]
    centers = []  # store the positions where bias is added

    for step in range(steps):
        # Langevin dynamics with force from the double well potential
        force = potential_force(x)
        thermal_force = np.sqrt(2 * gamma * temp / dt) * np.random.normal()

        # Apply bias from metadynamics
        force += bias_force(x, centers)

        # Update position with Langevin equation
        x += force * dt - gamma * x * dt + thermal_force * np.sqrt(dt)
        x_positions.append(x)

        # Add Gaussian bias every 100 steps
        if step % 100 == 0:
            centers.append(x)

    return np.array(x_positions)

# Run the simulation
xs = metaD()

# Create the figure with two subplots and shared x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

# Double well potential with Langevin dynamics
x_vals = np.linspace(-2, 2, 100)
potential_vals = double_well_potential(x_vals)
ax1.plot(x_vals, potential_vals, '--', lw=1.0, label="Double Well Potential", color="k")
sc = ax1.scatter(xs, double_well_potential(xs), c=range(len(xs)), s=2, cmap='viridis', alpha=0.5)
ax1.set_ylabel('Potential Energy')
ax1.legend()
ax1.set_title('Double Well Potential with MetaDynamics', fontweight='bold')
ax1.grid(True)
cbar = plt.colorbar(sc, ax=ax1, orientation='horizontal', pad=0.1)
cbar.set_label('Time Step')

# Histogram of x values
ax2.hist(xs, bins=50, color='skyblue', edgecolor='black')
ax2.set_xlabel('Position (x)')
ax2.set_ylabel('Frequency')
ax2.set_title('Histogram of Particle Positions')
ax2.grid(True)

plt.tight_layout()
plt.savefig('1D-MetaD.png')
```
<img src="https://github.com/qzhu2017/AtomisticSimulation/blob/main/Codes/lec_08_MetaD.png" alt="MD" height="400">

This simple illustration demonstrates how metadynamics escape from the local minima by using the bias potential. 


## 8.3 How Metadynamics Works?

In practice, metadynamics simulations can be applied to higher-dimensional systems and more complex potentials such as chemical reactions, protein folding, and phase transitions. It typically involves the following steps:

1. **Initialization**: Select appropriate CVs that describe the transition of interest. The choice of CVs is critical as they determine how effectively metadynamics can enhance sampling. CVs should capture the essential degrees of freedom involved in the transition, and their proper selection often requires prior knowledge of the system.

2. **Bias Addition**: During the simulation, small Gaussian potentials are periodically added along the CVs to the current position of the system. This gradually discourages the system from revisiting already visited states. The height and width of these Gaussians are important parameters that control the rate of exploration; a careful balance must be struck to ensure efficient sampling without overshooting important regions of the free energy landscape. 

3. **Exploration of the Free Energy Surface**: By continuously adding Gaussians, the system is encouraged to move out of energy minima, eventually allowing the exploration of the entire relevant free energy landscape. The bias potential effectively smooths out the energy barriers, enabling the system to transition between states more easily. Over time, the system visits all accessible regions of the free energy surface, and the accumulated bias provides an estimate of the free energy differences between states.


## 8.4 Choice of CVs

The efficiency of metadynamics heavily depends on the proper choice of CVs. Poor selection can lead to ineffective sampling, as the system may not be driven along the most relevant pathways. The process of selecting suitable CVs often requires trial and error or prior knowledge of the system.


## 8.5 Well-Tempered Metadynamics

In traditional metadynamics, the constant addition of Gaussian potentials can lead to excessive bias accumulation, which may result in poor sampling or an inaccurate free energy landscape. Well-tempered metadynamics addresses this by gradually reducing the height of the Gaussians added over time, which helps prevent oversampling and ensures that the system does not accumulate too much bias in any particular region. As such, it is expected to improve convergence and enhance the accuracy of the free energy surface estimation.

The key idea behind well-tempered metadynamics is to scale the bias deposition rate according to the amount of bias already present. This scaling is achieved by introducing a parameter called the bias factor $\gamma$, which controls how much the bias potential decreases as the simulation progresses. The bias factor is related to a fictitious temperature that effectively dictates how smoothly the bias is added.

The bias potential in well-tempered metadynamics evolves as:

$$
V_{\text{bias}}(s,t) = \sum_{t{\prime} \leq t} W \exp\left( -\frac{V_{\text{bias}}(s,t{\prime})}{k_B \Delta T} \right) \exp\left( -\frac{(s(t) - s(t{\prime}))^2}{2 \sigma^2} \right) 
$$

Compared to the previous equation, an additional exponential form was applied to scale the Gaussian potential, where 

- $k_B$ is the Boltzmann constant.
- $\Delta T$ is the fictitious temperature, defined as  $\Delta T = T_{\text{system}} (\gamma - 1)$, where $T_{\text{system}}$ is the real temperature of the system.

In the standard metadynamics, the height of Gaussian has a fixed value of $W$. In the Well-tempered metadynamics, the height becomes $W\exp\left[-{V_{\text{bias}}(s,t{\prime})}/{k_B \Delta T}\right].$ Recall that $V_{\text{bias}}$ is always positive and motonically increases with time. So $\exp \left[-V_{\text{bias}}(s,t{\prime})/{k_B \Delta T} \right]$ tends to gradually decay from 1 to a very small number as timestep increases. This would prevent a very large bias being deposited at one location. 

Thus, $\Delta T$ can be used to adjust how fast the decay should be. 
- When $\Delta T \rightarrow 0$, $\exp\left[-V_{\text{bias}}(s,t{\prime})/{k_B \Delta T} \right] \rightarrow 0$, which means a zero Gaussian height. Thus, the whole simulation returns to a standard MD simulation with a zero bias. In this scenario, one cannot see the barrier crossing. 
- When $\Delta T \rightarrow \infty$ , $\exp\left[-V_{\text{bias}}(s,t{\prime})/{k_B \Delta T} \right] \rightarrow 1$, which means a constant Gaussian height W is added at every update. This returns to a standard Metadynamics simulation. In this scenario, one can see the barrier crossing. But this leads to a even sampling of all $s$ values if one run the simulation for a very long time. 
- By choosing a suitable $\Delta T$ between 0 and infinity, one still keeps despositing nonzero Gussians. But the $s$ region with low $V$ values won't be filled quickly. Therefore, it ensures that low $V$ regions would be visited more frequently than those high $V$ regions.
- Hypothetically, one can imagine that the atomic coordinates $\mathbf{R}$ are weakly coupled with a set of collective variables $\mathbf{s}$, where $\mathbf{R}$ fluctuates around the system temperature $T$ wheras $\mathbf{s}$ fluctuatates around $T+\Delta T$. Therefore, the strong fluctuatation of $\mathbf{s}$ can help drag the system away from the minima while maintaining the same partition function for $\mathbf{R}$ under the given system temperature $T$.

In many applications, we are interested in knowing the free energy difference between different energy minima. One can infer the $F(s)$ by counting the histogram via $F(s) = -T \ln N(s, t)$, one eventually finds $F(s)$ is equal for any $\mathbf{s}$ in a standard metadynamics simulation. However, you will find $F(s)$ values converge to some distinct values in a Well-tempered metadynamics.
  
In summary, the Well-tempered metadynamics has several advantages over traditional metadynamics. 

1. By reducing the rate of bias deposition over time, it ensures that the system can **focus more on the most relevant regions of the free energy surface**, leading to a more accurate reconstruction, which is crucial for systems with multiple metastable states or complex free energy landscapes.

2. It provides a natural mechanism for achieving convergence of the free energy surface. As the system becomes more thoroughly explored, the bias added to the system decreases, ultimately reaching a point where it no longer significantly alters the free energy landscape. This gradual reduction in bias allows the system to settle into the correct free energy minima, providing a reliable estimate of the underlying free energy differences between states.

Well-tempered metadynamics can be applied to a wide range of systems, from simple model potentials to complex biomolecular and materials science applications. Its ability to adaptively control the bias potential makes it a versatile tool for studying processes such as protein folding, chemical reactions, and phase transitions. By providing a more controlled and convergent approach to free energy estimation, well-tempered metadynamics has become a preferred method for enhanced sampling in many challenging molecular simulations.

## 8.6 Further Discussions

Metadynamics is a powerful enhanced sampling technique to overcome the limitations of traditional MD simulation. 
Its ability to reconstruct free energy surfaces makes it an invaluable tool for studying complex molecular systems, phase transitions, and reaction mechanisms. However, its success depends on careful selection of collective variables and parameters. By systematically adding bias potentials, metadynamics allows for the study of rare events and provides detailed insights into the thermodynamics and kinetics of molecular systems. 

- Discuss how the choice of Gaussian parameters can impact the ordinary metadynamics simulation.
- Discuss the impact of $\gamma$ on the well tempered metadynamics simulation
- Explore the choice of CV in different kinds of simulations.

<!----References

1. Laio, A., & Parrinello, M. (2002). Escaping free-energy minima. Proceedings of the National Academy of Sciences, 99(20), 12562-12566.
2. Barducci, A., Bussi, G., & Parrinello, M. (2008). Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method. Physical Review Letters, 100(2), 020603.
---!>
