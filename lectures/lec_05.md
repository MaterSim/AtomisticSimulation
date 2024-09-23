# Week 5: MD Analysis I: Structural Characterization

So far, we have learned some fundamentals about how to develop a MD simulation to model the atomistic process of materials. By running the simulation, we expect to generate a set of time-dependent atomic trajectories by solving Newton’s equations of motion for a system of particles. Next, it is important to understand the simulation results. Indeed, post-analysis is essential for interpreting these trajectories to extract meaningful physical properties. In this lecture, we will cover several fundamental post-analysis techniques: evolution of observables, visualization of MD trajectories, the radial distribution function (RDF), and the vibrational spectrum.

## 5.1 Evolution of Macroscopic Observables 
In the previous lecture and coding exercises, we have frequently mentioned the tracking of observables like total energy, volume, and temperature throughout an MD simulation to ensure the system reaches equilibrium and conserves energy (in appropriate ensembles).

- Constant Energy (NVE): Total energy should remain constant.
- Constant Temperature (NVT): Temperature is maintained constant via a thermostat (a fluctuation of Temperature or kinetic energy around some values).
- Constant Pressure (NPT): Both Temperature and Volume fluctuate to maintain temperature and pressure.

In your MD simulation, it is advised to save these values every some time interval. Many codes would print out these information to some file. You must go over these results to ensure your simulation results make sense!

```python
import numpy as np
import matplotlib.pyplot as plt

# Load MD data (time, energy, volume, temperature)
time = np.loadtxt('time.dat')  # Time data
energy = np.loadtxt('energy.dat')  # Total energy
volume = np.loadtxt('volume.dat')  # System volume
temperature = np.loadtxt('temperature.dat')  # Temperature

# Plot evolution of observables
plt.figure(figsize=(12, 6))
plt.subplot(3, 1, 1)
plt.plot(time, energy, label='Total Energy', color='b')
plt.xlabel('Time (ps)')
plt.ylabel('Energy (eV)')
plt.title('Evolution of Total Energy')

plt.subplot(3, 1, 2)
plt.plot(time, volume, label='Volume', color='g')
plt.xlabel('Time (ps)')
plt.ylabel('Volume (Å³)')
plt.title('Evolution of Volume')

plt.subplot(3, 1, 3)
plt.plot(time, temperature, label='Temperature', color='r')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('Evolution of Temperature')
```

## 5.2 MD Trajectory Visualization
Visualization allows researchers to visually inspect the dynamics of the system, spot abnormalities, and better understand atomic movements. Tools like VMD and OVITO are commonly used. Please refer to OVITO pape to find more functions.

## 5.3 Radial Distribution Function (RDF)
The [RDF](https://en.wikipedia.org/wiki/Radial_distribution_function) $g(r)$ measures the probability of finding a particle at a distance $r$ from a reference particle.

$$
g(r) = \frac{V}{N^2} \left\langle \sum_{i=1}^{N} \sum_{j \neq i}^{N} \delta(r - r_{ij}) \right\rangle \cdot \frac{1}{4 \pi r^2 \Delta r}
$$

Where:

- $V$  is the volume of the system.
- $N$  is the number of particles.
- $r_{ij}$  is the distance between particles $i$ and $j$ .
- $\delta(r - r_{ij})$ is the Dirac delta function ensuring that only pairs with separation $r_{ij}$ equal to $r$ contribute.
- $4\pi r^2 \Delta r$ is the volume of a spherical shell at distance $r$ with thickness $\Delta r$.

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/3b/Rdf_schematic.svg/1280px-Rdf_schematic.svg.png" alt="Alt text" width="300"/>
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/31/Lennard-Jones_Radial_Distribution_Function.svg/1920px-Lennard-Jones_Radial_Distribution_Function.svg.png" alt="rdf" width="400"/>
</p>

It is commonly used to characterize the short-range order in liquids and gases. In a RDF, we are interested the location of Peaks and their spreads, as they indicate common interatomic distances (e.g., nearest-neighbor distance).

To compute $g(r)$, one should
1. **Compute the distance pair** $r_{ij}$ between every pair of particles $i$ and $j$ .
2. **Group the distances** into bins of width $\Delta r$ to accumulate how many pairs of particles have a separation within each bin.
3. **Normalize the number of pairs** in each bin by:
- The number of particle pairs $N(N-1)/2$.
- The volume of the spherical shell at distance $r$, given by $4\pi r^2 \Delta r$.
- The average particle density $N/V$.

```Python
def compute_rdf(positions, num_bins=100, r_max=10.0):
    N = len(positions)  # Number of atoms
    rdf = np.zeros(num_bins)
    dr = r_max / num_bins
    for i in range(N):
        for j in range(i+1, N):
            r = np.linalg.norm(positions[i] - positions[j])
            if r < r_max:
                bin_index = int(r / dr)
                rdf[bin_index] += 2  # Each pair counted twice

    # Normalize RDF
    r = np.linspace(0, r_max, num_bins)
    rdf /= (4 * np.pi * r**2 * dr * N)
    return r, rdf

# Example usage
positions = np.random.rand(100, 3) * 10  # Generate random positions
r, g_r = compute_rdf(positions)

# Plot RDF
plt.plot(r, g_r)
plt.xlabel('r (Angstrom)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.show()
```

## 5.4 Vibration Spectrum
In addition to RDF, another commond characterization is to understand how particles in a system vibrate. In experiment, such information can be measured from Infrared (IR) or Raman Spectroscopy, and Inelastic Neutron/X-ray Scattering. An analogical measurement in MD is the **Vibrational Density of States** (VDOS). 

In a MD simulation, the VDOS ($D(\omega)$) is essentially the frequency spectrum of atomic vibrations in the system. The VDOS is closely related to the **Velocity Autocorrelation Function** (VACF). The VDOS provides information about the frequencies at which particles in a system vibrate, while the VACF gives insight into how the velocity of a particle at one time is correlated with its velocity at a later time. It is computed by taking the Fourier transform of the VACF. This relationship is derived from the fact that oscillatory motions (vibrations) in the system are directly reflected in the velocity autocorrelation function, and the Fourier transform allows us to extract the frequency components of these vibrations.

Mathematically, the relationship is:

$$
D(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \text{VACF}(\tau) e^{-i\omega\tau}  d\tau
$$

Where:

- $\text{VACF}(\tau)$ is the velocity autocorrelation function.
- $\omega$ is the angular frequency.
- $\tau$ is the time lag.

The Fourier transform is performed over the time correlation $\tau$ to convert the time-domain information in the VACF into frequency-domain information in the VDOS.

The VACF measures how the velocity of a particle at a given time $t$ correlates with its velocity at some later time $t + \tau$. It is useful for understanding particle dynamics and is related to the vibrational properties and transport coefficients (like diffusion).

$$
\text{VACF}(\tau) = \frac{1}{N} \sum_{i=1}^{N} \left\langle \mathbf{v}_i(0) \cdot \mathbf{v}_i(\tau) \right\rangle
$$

- $N$ is the total number of time steps in the simulation.
- $\mathbf{v}_i(t)$ is the velocity of particle $i$ at time step $t$.

In practice, because simulations are finite and VACF data is computed over a limited time interval, we typically use the discrete Fourier transform (DFT) or fast Fourier transform (FFT) to compute the VDOS numerically:

$$
D(\omega) = \frac{1}{2\pi} \int_0^{\infty} \text{VACF}(\tau) \cos(\omega \tau) d\tau
$$


```python
from scipy.fft import fft

def compute_vacf(velocities):
    vacf = np.correlate(velocities, velocities, mode='full')
    return vacf[vacf.size // 2:]  # Only take positive lag times

# Example: Simulate random velocities
velocities = np.random.randn(1000)

# Compute VACF
vacf = compute_vacf(velocities)

# Compute vibration spectrum via Fourier Transform
vibration_spectrum = np.abs(fft(vacf))

# Plot vibration spectrum
plt.plot(vibration_spectrum[:len(vibration_spectrum)//2])
plt.xlabel('Frequency (THz)')
plt.ylabel('Intensity')
plt.title('Vibration Spectrum')
plt.show()
```

## 5.5 Further discussions

- Interpretation of RDF: Peaks in RDF tells the atomic neighbor counts. Liquid and solid have very different behaviors in their RDF.  
- Interpretation of VDOS: Peaks in the VDOS correspond to characteristic vibrational modes of the system. Try to identify the Low-frequency and high-frequency modes and link them to atomic motions from MD trajectory
