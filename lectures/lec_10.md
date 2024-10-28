# 10 DFT Simulation of a Hydrogen Molecule

Having learned the basic concept of DFT, we will continue to apply the DFT approach to simulate a more complicated system than a single harmonic oscillator, i.e., the H<sub>2</sub> molecule. H<sub>2</sub> is the simplest molecule, consisting of two protons and two electrons. This small size is ideal for demonstrating the DFT method in a manageable way.
By solving for the ground state energy and electron density of H<sub>2</sub>, we hope to better understand the numerical aspects of the DFT method. In addition, we will analyze the results to understand the bonding in H$_2$ and the
electron distribution.

## 10.1 Basic Setup

In a H<sub>2</sub>  molecule, we need to consider the following variables in the KS equation.

- Nuclei: Two protons, located at positions $R_1$ and $R_2$ .
- Electrons: Two electrons interacting with the protons and each other.

We first write down its Hamiltonian,

$$
\hat{H} = \hat{T}_e + \hat{V} _{\text{ext}} + \hat{V} _{\text{ee}} + \hat{V} _{\text{nn}}
$$

Where
- $\hat{T}_e$ : Kinetic energy of the electrons.
- $\hat{V}_{\text{ext}}$ : External potential due to nuclei.
- $\hat{V}_{\text{ee}}$ : Electron-electron interaction.
- $\hat{V}_{\text{nn}}$ : Nucleus-nucleus interaction.

Hence, Kohn-Sham formalism reduces the many-body problem to a series of single-particle equations. For the case of H<sub>2</sub>, there is only one Kohn-Sham equation solved because both electrons occupy the same Kohn-Sham orbital. We assume that both electrons have opposite spins and fill the same orbital, so the total electron density is computed as $\rho(x) = 2 |\phi_0(x)|^2$. **If you want to solve a spin-polarized (or unrestricted) system, you would need to solve two distinct Kohn-Sham equations, one for each electron.**

$$
\left[ -\frac{\hbar^2}{2m} \nabla^2 + V_{\text{eff}}(\mathbf{r}) \right] \phi_i(\mathbf{r}) = \epsilon_i \phi_i(\mathbf{r})
$$

1. The total energy is
   
$$
E[\rho] = T_s[\rho] + E_{\text{ext}}[\rho] + E_H[\rho] + E_{\text{XC}}[\rho]
$$

Where:
- $T_s[\rho]$ : Kinetic energy of non-interacting electrons.
- $E_{\text{ext}}[\rho]$ : Interaction with external potential.
- $E_H[\rho]$ : Hartree energy (electron-electron repulsion).
- $E_{\text{XC}}[\rho]$ : Exchange-correlation energy.

2. The external potential  $V_{\text{ext}}(\mathbf{r})$  for each electron is given by the Coulomb interaction with the two protons:

$$
V_{\text{ext}}(\mathbf{r}) = -\frac{1}{|\mathbf{r} - \mathbf{R}_1|} - \frac{1}{|\mathbf{r} - \mathbf{R}_2|}
$$

## 10.2 Python Implementation

1. Define a spatial grid and use the finite difference method to approximate the kinetic energy operator.
```python
import numpy as np
import matplotlib.pyplot as plt

# Spatial grid parameters
x_min, x_max = -10.0, 10.0  # Grid boundaries
N = 1000  # Number of grid points
x = np.linspace(x_min, x_max, N)
dx = x[1] - x[0]  # Grid spacing

# Define external potential for H2 molecule
R1, R2 = -1.0, 1.0  # Proton positions
V_ext = -1 / np.abs(x - R1) - 1 / np.abs(x - R2)
```

2. Initialize Electron Density
```python
# Initial guess for electron density (uniform)
rho = np.ones(N) * 0.01
```

3. Compute Effective Potential $V_{\text{eff}}(x)$ including 
- external potential
- Hartree potential
- exchange-correlation potential.

```python
def compute_hartree_potential(rho, dx):
    """Compute Hartree potential from the electron density."""
    V_H = np.convolve(rho, 1 / np.abs(x - x[:, None]), mode='same') * dx
    return V_H

# Exchange-correlation potential (local density approximation)
def compute_exchange_correlation_potential(rho):
    return -(3 / np.pi)**(1/3) * rho**(1/3)

def compute_effective_potential(rho, V_ext, dx):
    V_H = compute_hartree_potential(rho, dx)
    V_XC = compute_exchange_correlation_potential(rho)
    return V_ext + V_H + V_XC
```

4. Iterative update with the self-consistent field (SCF) approach

```python
# Kinetic energy operator (finite difference approximation)
T = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2

def solve_kohn_sham(V_eff, T):
    # Hamiltonian
    H = -0.5 * T + np.diag(V_eff)
    energies, orbitals = np.linalg.eigh(H)
    return energies, orbitals

# Solve self-consistently
for iteration in range(100):
    V_eff = compute_effective_potential(rho, V_ext, dx)
    energies, orbitals = solve_kohn_sham(V_eff, T)

    # Update electron density (2 electrons in total, filling the lowest orbital)
    rho_new = 2 * np.abs(orbitals[:, 0])**2

    # Check for convergence
    if np.linalg.norm(rho_new - rho) < 1e-6:
        print(f"Converged after {iteration+1} iterations")
        break

    rho = rho_new
```
The code iteratively solves the KS equations. In each iteration, the electron density is updated from the orbitals, and the effective potential is recalculated using the new density. This process continues until the electron density converges (i.e., when the difference between the new and old density is smaller than a set tolerance).

Here the function ``np.linalg.eigh(H)`` solves the eigenvalue problem for the Hamiltonian, returning the Kohn-Sham energies and orbitals. Once the Hamiltonian is constructed and solved, the Kohn-Sham orbitals  $\phi_i(x)$  are used to update the electron density  $\rho(x)$ . The SCF loop continues until the electron density converges.

5. Compute Total Energy

```python
def compute_total_energy(rho, V_ext, V_eff, T, orbitals, dx):
    # Kinetic energy
    T_s = -0.5 * np.sum(orbitals[:, 0] * np.dot(T, orbitals[:, 0])) * dx

    # External potential energy
    E_ext = np.sum(rho * V_ext) * dx

    # Hartree energy
    E_H = 0.5 * np.sum(rho * compute_hartree_potential(rho, dx)) * dx

    # Exchange-correlation energy (LDA)
    E_XC = np.sum(-(3 / np.pi)**(1/3) * rho**(4/3)) * dx

    return T_s + E_ext + E_H + E_XC

total_energy = compute_total_energy(rho, V_ext, V_eff, T, orbitals, dx)
print(f"Total Energy of H2 molecule: {total_energy:.4f} Hartree")
```

6. Visualize the Electron Density
```python
plt.plot(x, rho)
plt.title('Electron Density for H2 Molecule')
plt.xlabel('x (Bohr)')
plt.ylabel('Electron Density')
plt.show()
```
## 10.3 GGA implementation 
In the previous example, we used Local Density Approximation (LDA) for the exchange-correlation energy $E_{\text{XC}}$.

$$
E_{\text{XC}}^{\text{LDA}}[\rho] = \int \epsilon_{\text{XC}}(\rho(\mathbf{r})) \rho(\mathbf{r})  d\mathbf{r}
$$

The Generalized Gradient Approximation (GGA) improves upon the LDA by including not only the electron density  $\rho(\mathbf{r})$  but also the gradient of the electron density $\nabla \rho(\mathbf{r})$. This allows for a more accurate representation of the exchange-correlation effects, especially in systems where the electron density varies significantly, such as in molecular systems like H<sub>2</sub>.

$$
E_{\text{XC}}^{\text{GGA}}[\rho] = \int f(\rho(\mathbf{r}), \nabla \rho(\mathbf{r})) d\mathbf{r}
$$

where $f(\rho, \nabla \rho)$ is a function that depends on both the local density $\rho$ and its gradient $\nabla \rho$.

To implement GGA, we need to modify the exchange-correlation potential $V_{\text{XC}}$  and the corresponding exchange-correlation energy  $E_{\text{XC}}$. For simplicity, we will use the PBE (Perdew-Burke-Ernzerhof) form of GGA, which is widely used in DFT calculations. In this case, both the exchange-correlation energy and potential are functions of $\rho(\mathbf{r})$ and $|\nabla \rho(\mathbf{r})|$.

```python
def compute_density_gradient(rho, dx):
    # Central difference for gradient
    grad_rho = (np.roll(rho, -1) - np.roll(rho, 1)) / (2 * dx)
    return grad_rho

def compute_exchange_correlation_gga(rho, grad_rho):
    # PBE-like GGA exchange-correlation energy density
    # Simplified for the purpose of this example
    A = -(3 / np.pi)**(1/3)  # Constant for exchange energy
    k = 0.804  # Gradient correction term (simplified)

    # LDA part
    exc_lda = A * rho**(4/3)

    # GGA part (simplified version)
    s = np.abs(grad_rho) / rho**(4/3)  # Reduced gradient
    exc_gga = exc_lda * (1 + k * s**2)

    return exc_gga

def compute_total_energy_gga(rho, V_ext, T, orbitals, dx):
    # Kinetic energy (same as before)
    T_s = -0.5 * np.sum(orbitals[:, 0] * np.dot(T, orbitals[:, 0])) * dx

    # External potential energy (same as before)
    E_ext = np.sum(rho * V_ext) * dx

    # Hartree energy (same as before)
    E_H = 0.5 * np.sum(rho * compute_hartree_potential(rho, dx)) * dx

    # Compute the density gradient for GGA
    grad_rho = compute_density_gradient(rho, dx)

    # GGA exchange-correlation energy
    E_XC = np.sum(compute_exchange_correlation_gga(rho, grad_rho)) * dx

    return T_s + E_ext + E_H + E_XC
```


## 10.4 Further discussions

- Interpret the Ground State Energy
- Compare the computed total energy of the H<sub>2</sub> molecule with known reference values.
- Discuss how the electron density represents the covalent bonding between the two protons.
- Compare the results obtained from LDA and GGA
- Analyze how DFT with GGA provides a more accurate description of the electron density and energy.
- Discuss the numerical aspects when extending DFT to more complicated systems (molecules or crystals)

