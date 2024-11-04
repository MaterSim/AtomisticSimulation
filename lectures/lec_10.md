# 10 DFT Simulation of a Hydrogen Molecule
Having learned the basic concept of DFT, we will continue to apply the DFT approach to simulate a more complicated system than a single electron system, i.e., the H<sub>2</sub> molecule. H<sub>2</sub> is the simplest molecule, consisting of two protons and two electrons. This small size is ideal for demonstrating the DFT method in a manageable way.

By solving for the ground state energy and electron density of H<sub>2</sub>, we hope to better understand the numerical aspects of the DFT method. In addition, we will analyze the results to understand the bonding in H<sub>2</sub> and the electron distribution.

## 10.1 Basic Setup

In an H<sub>2</sub> molecule, we need to consider the following variables in the KS equation.

- Nuclei: two protons, located at positions $R_1$ and $R_2$.
- Electrons: two electrons interacting with the protons and each other.

We first write down its Hamiltonian,

$$
\hat{H} = \hat{T}_e + \hat{V} _{\text{external}} + \hat{V} _{\text{ee}} + \hat{V} _{\text{nn}}
$$

Where
- $\hat{T}_e$ : Kinetic energy of the electrons.
- $\hat{V}_{\text{ext}}$ : External potential due to nuclei.
- $\hat{V}_{\text{ee}}$ : Electron-electron interaction.
- $\hat{V}_{\text{nn}}$ : Nucleus-nucleus interaction.

Hence, the Kohn-Sham formalism reduces the many-body problem to a series of single-particle equations. For the case of H<sub>2</sub>, there is only one Kohn-Sham equation solved because both electrons occupy the same Kohn-Sham orbital. We assume that both electrons have opposite spins and fill the same orbital, so the total electron density is computed as $\rho(\mathbf{r}) = 2 |\phi_0(\mathbf{r})|^2$. **If you want to solve a spin-polarized (or unrestricted) system, you would need to solve two distinct Kohn-Sham equations, one for each electron.**

$$
\left[ -\frac{\hbar^2}{2m} \nabla^2 + V_{\text{eff}}(\mathbf{r}) \right] \phi_i(\mathbf{r}) = \epsilon_i \phi_i(\mathbf{r})
$$

## 10.2 Effective Potentials

### 10.2.1 The external potential $V_{\text{ext}}(\mathbf{r})$

In a hydrogen molecule, the external potential is given by the Coulomb interaction with the two protons:

$$
V_{\text{external}}(\mathbf{r}) = -\frac{1}{|\mathbf{r} - \mathbf{R}_1|} - \frac{1}{|\mathbf{r} - \mathbf{R}_2|}
$$

To avoid the singularity in the Coulomb interaction for each electron with the protons, the external potential  $V_{\text{external}}(\mathbf{r})$ can be softened as follows:

$$
V_{\text{external}}(\mathbf{r}) = -\frac{1}{\sqrt{|\mathbf{r} - \mathbf{R}_1|^2 + \alpha^2}} - \frac{1}{\sqrt{|\mathbf{r} - \mathbf{R}_2|^2 + \alpha^2}}
$$

where $\alpha$ is a small positive softening parameter. This parameter helps to avoid the singularity at points $\mathbf{r} = \mathbf{R}_1$  and  $\mathbf{r} = \mathbf{R}_2$, where the distance between the electron and the nuclei would otherwise be zero, causing a divergence in the potential.

### 10.2.2 Hartree potential and Energy

The Hartree potential represents the classical electrostatic potential at a point $\mathbf{r}$ due to the electron density distribution  $\rho(\mathbf{r})$. It can be calculated as:

$$
V_{\text{Hartree}}(\mathbf{r}) = \int \frac{\rho(\mathbf{r}{\prime})}{|\mathbf{r} - \mathbf{r}{\prime}|}  d\mathbf{r}{\prime}
$$

Direct evaluation of $V_{\text{Hartree}}(\mathbf{r})$  in real space can be computationally expensive, scaling up to $O(N^6)$ for large systems. To address this, the Fourier Space method is often preferred, especially in the periodic systems, as it efficiently computes the potential by solving the Poisson equation in reciprocal space. The Hartree potential can be derived from the Fourier-transformed electron density $\rho(\mathbf{k})$ as:

$$
V_{\text{Hartree}}(\mathbf{k}) = \frac{4 \pi \rho(\mathbf{k})}{|\mathbf{k}|^2}
$$

where $V_{\text{Hatree}}(\mathbf{k})$ is the Fourier-transformed Hartree potential.

This involves three steps:

1.	Fourier transform of electron density:

$$
\rho(\mathbf{k}) = \text{FFT}(\rho(\mathbf{r}))
$$

2.	Calculate the Hartree Potential in Fourier Space (avoid division by zero at  $\mathbf{k}$ = 0):

$$
V_{\text{Hartree}}(\mathbf{k}) = \frac{4 \pi \rho(\mathbf{k})}{|\mathbf{k}|^2}
$$


3.	Inverse Fourier Transform to Real Space:

$$
V_{\text{Hartree}}(\mathbf{r}) = \text{IFFT}(V_{\text{Hartree}}(\mathbf{k}))
$$


Correspondingly, the Hartree Energy is 

$$
E_{\text{Hatree}} = \frac{1}{2} \int \rho(\mathbf{r}) V_{\text{Hatree}}(\mathbf{r}) d\mathbf{r}
$$


### 10.2.3 Exchange and Correlation 

In a hydrogen molecule, the exchange and correlation effects account for the complex quantum interactions between electrons, beyond the simple electrostatic repulsion captured by the Hartree term. These effects are crucial in accurately describing the binding energy, molecular structure, and overall electronic distribution within the molecule.

The exchange-correlation (XC) energy, $E_{\text{XC}}$, combines both exchange effects, which arise from the antisymmetry requirement of the electron wavefunction, and correlation effects, which account for electron-electron repulsion. The exchange-correlation potential, $V_{\text{XC}}(\mathbf{r})$, is defined as the functional derivative of $E_{\text{XC}}$ with respect to the electron density $\rho(\mathbf{r})$:

$$
V_{\text{XC}}(\mathbf{r}) = \frac{\delta E_{\text{XC}}[\rho]}{\delta \rho(\mathbf{r})}
$$


In the Local Density Approximation (LDA), the exchange-correlation energy density at each point is assumed to depend only on the local electron density. For a given density $\rho(\mathbf{r})$, the LDA exchange-correlation energy is:

$$
E_{\text{XC}}^{\text{LDA}}[\rho] = \int \epsilon_{\text{XC}}(\rho(\mathbf{r})) \rho(\mathbf{r}) d\mathbf{r}
$$

where $\epsilon_{\text{XC}}(\rho)$ is the exchange-correlation energy per electron of a uniform electron gas with density  $\rho$.

$$
\epsilon_{\text{X}}(\rho) = -\frac{3}{4} \left( \frac{3}{\pi} \right)^{1/3} \rho^{1/3}
$$

For the case of simplicity, we ignore the correlation term here.

### 10.2.4 The Expression of Effective Potential and Total Energy

Adding all potential terms together, the effective potential is

$$
V_{\text{eff}}(\mathbf{r}) = V_{\text{external}}(\mathbf{r}) + V_{\text{Hartree}}[\rho] + V_{\text{XC}}(\mathbf{r})
$$

The total energy is the sum of potential and kinetic energies. The kinetic energy $T_s$ can be numerically evaluated as

$$
T_s = -\frac{1}{2} \sum_i^{\text{occupied}} \int \phi_i^*(\mathbf{r}) \nabla^2 \phi_i(\mathbf{r})  d\mathbf{r}
$$

where $\phi_i(\mathbf{r})$  are the Kohn-Sham orbitals, and the sum runs over all occupied orbitals.

Thus, the total energy is

$$
E[\rho] = T_s[\rho] + E_{\text{ext}}[\rho] + E_H[\rho] + E_{XC}[\rho] + E_{\text{NN}}
$$

where $E_{\text{NN}}$ is the nuclear-nuclear repulsion energy, given by $\frac{e^2}{|R_2 - R_1|}$ for two protons at $R_1$, $R_2$.


## 10.3 Python Implementation

Below gives an example Python code to compute the ground state of hydrogen molecule.

### 10.3.1 DFT-SCF calculation
```python
import numpy as np
from scipy.sparse import kron, eye
from scipy.sparse import diags, csr_matrix
from scipy.fft import fftn, ifftn
from scipy.sparse.linalg import eigsh  # Sparse eigenvalue solver for faster performance
import matplotlib.pyplot as plt

# Define 3D grid parameters
N = 80 #120  # Grid size along each dimension
L = 6  #40  # Simulation box size in Bohr
dx = L / N
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
z = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

# Positions of the two protons in H2
R1 = np.array([-0.7, 0, 0])
R2 = np.array([0.7, 0, 0])
num_electrons = 2.0
E_nuc = 1.0 / (R2[0]-R1[0])

# Softened Coulomb potential for two protons
softening = 0.1  # To prevent singularity 
V_ext = -1 / np.sqrt((X - R1[0])**2 + (Y - R1[1])**2 + (Z - R1[2])**2 + softening**2)
V_ext += -1 / np.sqrt((X - R2[0])**2 + (Y - R2[1])**2 + (Z - R2[2])**2 + softening**2)

# Initial guess for electron density
rho = np.exp(-1 * ((X - R1[0])**2 + (Y - R1[1])**2 + (Z - R1[2])**2))
rho += np.exp(-1 * ((X - R2[0])**2 + (Y - R2[1])**2 + (Z - R2[2])**2))
total_density = np.sum(rho) * dx**3
rho *= num_electrons / total_density

# Define FFT-based Poisson solver for Hartree potential
def compute_hartree_potential(rho):
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    kz = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1  # Avoid division by zero for the zero-frequency term

    rho_k = fftn(rho)
    V_H_k = 4 * np.pi * rho_k / k2
    V_H_k[0, 0, 0] = 0  # Set the zero-frequency term to zero for neutrality

    V_H = np.real(ifftn(V_H_k))
    return V_H #* 2

# Exchange-correlation potential (Local Density Approximation)
def compute_exchange_correlation_potential(rho):
    c = (3 / np.pi)**(1 / 3)
    return -c * rho**(1 / 3)

def kinetic_energy_operator(N, dx):
    # Define 1D kinetic energy finite difference operator
    main_diag = -2 * np.ones(N)
    side_diag = np.ones(N - 1)
    T_1D = diags([main_diag, side_diag, side_diag], [0, -1, 1], shape=(N, N)) / dx**2

    # Build 3D kinetic energy operator using Kronecker products
    I = eye(N, format='csr')  # Identity matrix for each dimension
    T = kron(kron(T_1D, I), I) + kron(kron(I, T_1D), I) + kron(kron(I, I), T_1D)

    return -0.5 * T  # Scale by -0.5 

# SCF loop parameters
tolerance = 1e-6
max_iterations = 100
damping_factor = 0.5

# Kinetic energy operator (sparse)
T = kinetic_energy_operator(N, dx); print(T.shape)

for iteration in range(max_iterations):
    # Compute Hartree and XC potentials
    V_H = compute_hartree_potential(rho)
    V_XC = compute_exchange_correlation_potential(rho)
    V_eff = V_ext + V_H + V_XC

    # Construct H
    V_eff_flat = V_eff.flatten()
    H = T + diags(V_eff_flat, 0, shape=(N**3, N**3))

    # Solve the Kohn-Sham equation (using sparse eigenvalue solver for the lowest eigenvalue)
    energies, orbitals = eigsh(H, k=1, which='SA')  
    psi = orbitals[:, 0].reshape((N, N, N))  
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx**3)

    # Update electron density
    rho_new = num_electrons * np.abs(psi)**2  # Two electrons in ground state

    # Calculate energies for this iteration
    T_s = 2* np.sum(orbitals[:, 0] * T.dot(orbitals[:, 0])) * dx**3
    E_H = 0.5 * np.sum(rho * V_H) * dx**3
    E_ext = np.sum(rho * V_ext) * dx**3; print(rho.sum()*dx**3, rho.shape)
    E_XC = np.sum(rho * V_XC) * dx**3

    # Print the energies
    print(f"Iteration {iteration + 1}:")
    print(f"Kinetic Energy       (T)    = {T_s:.6f} Hartree")
    print(f"External Energy      (E_ext)= {E_ext:.6f} Hartree")
    print(f"Hartree Energy       (E_H)  = {E_H:.6f} Hartree")
    print(f"Exchange-Correlation (E_XC) = {E_XC:.6f} Hartree")
    print(f"Total Effective Energy      = {T_s + E_ext + E_H + E_XC + E_nuc:.6f} Hartree")
    print(f"Solver KS HOMO Energy       = {energies[0]:.6f} Hartree\n")
    # Damping update
    rho = (1 - damping_factor) * rho + damping_factor * rho_new
    total_density = np.sum(rho) * dx**3

    total_electrons = np.sum(rho) * dx**3
    print(f"Total electrons: {total_electrons}")

    #XC_scaling_factor = min(XC_scaling_factor + 0.1, 1.0)  # Gradually ramp up to 1.0
    # Check convergence
    if np.linalg.norm(rho_new - rho) < tolerance:
        print(f"Converged after {iteration + 1} iterations")
        break
else:
    print("Did not converge within the maximum number of iterations")
```

The code iteratively solves the KS equations. In each iteration, the electron density is updated from the orbitals, and the effective potential is recalculated using the new density. This process continues until the electron density converges (i.e., when the difference between the new and old density is smaller than a set tolerance).

Here the function `eigsh` solves the eigenvalue problem for the Hamiltonian, returning the Kohn-Sham energies and orbitals. Once the Hamiltonian is constructed and solved, the Kohn-Sham orbitals $\phi_i(\mathbf{r})$ are used to update the electron density $\rho(\mathbf{r})$. The SCF loop continues until the electron density converges.

An example output looks like the following

```
Iteration 1:
Kinetic Energy       (T)    = 1.263639 Hartree
External Energy      (E_ext)= -3.557743 Hartree
Hartree Energy       (E_H)  = 0.545992 Hartree
Exchange-Correlation (E_XC) = -0.819774 Hartree
Total Effective Energy      = -1.853602 Hartree
Solver HOMO Energy          = -1.046968 Hartree

Iteration 2:
Kinetic Energy       (T)    = 1.256039 Hartree
External Energy      (E_ext)= -3.599140 Hartree
Hartree Energy       (E_H)  = 0.545473 Hartree
Exchange-Correlation (E_XC) = -0.812328 Hartree
Total Effective Energy      = -1.895671 Hartree
Solver HOMO Energy          = -1.046764 Hartree

Iteration 3:
Kinetic Energy       (T)    = 1.253817 Hartree
External Energy      (E_ext)= -3.614615 Hartree
Hartree Energy       (E_H)  = 0.543885 Hartree
Exchange-Correlation (E_XC) = -0.809045 Hartree
Total Effective Energy      = -1.911672 Hartree
Solver HOMO Energy          = -1.047270 Hartree
```

To provide a comparison for convergence values in the Kohn-Sham DFT iterations, here are the typical reference values for the energy components of the hydrogen molecule:

- Kinetic Energy ($T$): Approximately 1.2–1.3 Hartree
- External Energy ($E_{\text{ext}}$): between -3.4 and -3.8 Hartree
- Hartree Energy ($E_{\text{H}}$): Usually around 1.2 Hartree
- Exchange-Correlation Energy ($E_{\text{XC}}$): Typically between -0.8 and -0.9 Hartree
- Total Energy: around -1.19 Hartree
- The HOMO energy (highest occupied molecular orbital) is generally around -0.6 Hartree.

Increasing $L$ and $N$ helps to mitigate edge effects, especially for long-range interactions like the Hartree term, and provides a more accurate representation of the continuum electron density, ensuring calculations remain closer to these benchmark values.

### 10.3.2 Physical Interpretation Analysis

One can also plot electron density to see if the calculation can successfully reproduce the well-known bonding and anti-bonding picture from the ground and first excited states.

```python
energies, orbitals = eigsh(H, k=2, which='SA')
psi = orbitals[:, 0].reshape((N, N, N))  # Ground state orbital
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx**3)
rho = np.abs(psi)**2  # Two electrons in ground state

# Plot a cross-section of the electron density along the z=0 plane
plt.contourf(rho[:, :, N//2], extent=(-L/2, L/2, -L/2, L/2), origin='lower')
plt.colorbar(label='Electron Density')
plt.xlabel('x (Bohr)')
plt.ylabel('y (Bohr)')
plt.title('Electron Density for H2 Molecule in 3D')
plt.savefig('H2-bond.png')
plt.close()

psi = orbitals[:, 1].reshape((N, N, N))  # Ground state orbital
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx**3)
rho = np.abs(psi)**2  # Two electrons in ground state
rho *= num_electrons / (np.sum(rho) * dx**3)

plt.contourf(rho[:, :, N//2], extent=(-L/2, L/2, -L/2, L/2), origin='lower')
plt.colorbar(label='Electron Density')
plt.xlabel('x (Bohr)')
plt.ylabel('y (Bohr)')
plt.title('Electron Density for H2 Molecule in 3D')
plt.savefig('H2-antibond.png')
plt.close()
```

<p align="center">
  <img src="https://github.com/qzhu2017/AtomisticSimulation/blob/main/Codes/lec_10_bonding.png" alt="Alt text" width="480"/>
  <img src="https://github.com/qzhu2017/AtomisticSimulation/blob/main/Codes/lec_10_antibonding.png" alt="Alt text" width="480"/>
</p>

## 10.4 Further discussions

- **Interpret the Ground State Energy**. The ground state energy of the H<sub>2</sub> molecule is a critical value that represents the lowest energy configuration of the system. It includes contributions from the kinetic energy of the electrons, the external potential energy due to the nuclei, the Hartree energy representing electron-electron repulsion, and the exchange-correlation energy. 

- **Compare the Computed Energy Values**. The computed total energy of the H<sub>2</sub> molecule should be compared with known reference values from literature or high-precision quantum chemistry calculations. Typically, the total energy for the H<sub>2</sub> molecule is around -1.17 Hartree. Any significant deviation from this value could indicate issues with the numerical setup or the approximations used.

- **Discuss the covalent bonding**. The electron density distribution in the H<sub>2</sub> molecule shows how electrons are shared between the two protons, forming a covalent bond. In the ground state, the electron density is highest between the two nuclei, indicating a strong bonding interaction. This distribution can be visualized through contour plots or 3D density plots.

- **Compare the Results Obtained from LDA and GGA**. The Local Density Approximation (LDA) and the Generalized Gradient Approximation (GGA) are two different approaches to approximating the exchange-correlation energy in DFT. LDA considers only the local electron density, while GGA includes the gradient of the electron density. GGA typically provides more accurate results for systems with varying electron densities, such as molecules and surfaces.

- **Extending DFT to More Complicated Systems (Molecules or Crystals)**. Extending DFT to more complicated systems, such as larger molecules or crystals, involves several numerical challenges. These include the increased computational cost due to the larger number of electrons and nuclei, the need for more sophisticated algorithms to solve the Kohn-Sham equations efficiently, and the handling of periodic boundary conditions in crystals. Advanced techniques such as plane-wave basis sets, pseudopotentials, and parallel computing are often employed to address these challenges.


<!----
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
----!>

