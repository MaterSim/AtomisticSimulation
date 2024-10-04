# 9 Introduction to Density Functional Theory

## 9.1 History
In 1920-1930, Schrödinger Equation provides a way to calculate the wavefunction $\Psi$ of a quantum system, encapsulating all its physical properties. However, the exact solutions exist only for simple systems like the hydrogen atom. Solving the Schrödinger equation for systems with more than one electron becomes intractable due to the electron-electron interactions. Hartree and Hartree-Fock Methods were developed to apply the mean-field approximations to simplify the many-body problem but still faced limitations, particularly in accounting for electron correlation. 

In 1964, Pierre Hohenberg and Walter Kohn established Density Functional Theory (DFT) that the ground-state properties of a many-electron system are uniquely determined by its electron density $\rho(\mathbf{r})$. Later in 1965, Walter Kohn and Lu Jeu Sham developed a practical method to implement DFT by introducing a system of non-interacting electrons that reproduce the same electron density as the interacting system. Nowadays, DFT has become an indispensable tool in physics, chemistry, and materials science for studying the electronic structure of matter. 

## 9.2 The Many-Body Problem in Quantum Mechanics 

The time-independent Schrödinger equation for a system of $N$ interacting electrons is:

$$
\hat{H} \Psi(\mathbf{r}_1, \mathbf{r}_2, \cdots, \mathbf{r}_N) = E \Psi(\mathbf{r}_1, \mathbf{r}_2, \cdots, \mathbf{r}_N)
$$

- $\hat{H}$ is the Hamiltonian operator, including kinetic and potential energies.
- $\Psi$ is the many-electron wavefunction, depending on the positions of all electrons.
- $E$ is the total energy of the system.

The Hamiltonian for such a system is:

$$
\hat{H} = \sum_{i=1}^{N} \left( -\frac{\hbar^2}{2m} \nabla_i^2 + V_{\text{ext}}(\mathbf{r}_i) \right) + \sum^{i \lt j} \frac{e^2}{4\pi\epsilon_0 |\mathbf{r}_i - \mathbf{r}_j|}
$$

where 

- The first term represents the kinetic energy of each electron and $V_{\text{ext}}$ (due to external potential).
- The second summation accounts for electron-electron Coulomb interactions.

It is not easy to solve this equation due to the following issues,
- The wavefunction $\Psi$ depends on $3N$ spatial coordinates, making it computationally infeasible for large $N$.
- Computational resources required grow exponentially with the number of electrons.
- Properly accounting for interactions and correlations between electrons is challenging.
- Methods like Hartree-Fock simplify the problem but neglect electron correlation, leading to inaccuracies.

## 9.3 The Hohenberg-Kohn Theorems

Kohn and Hohenberg proofed two Theorems
1. The ground-state electron density $\rho_0(\mathbf{r})$ uniquely determines the external potential $V_{\text{ext}}(\mathbf{r})$, up to a constant, and hence all properties of the system. 

2. There exists a universal functional $E[\rho]$ such that the ground-state energy $E_0$ can be obtained variationally:

$$
E_0 = \min_{\rho} \left[ E[\rho] = F[\rho] + \int V_{\text{ext}}(\mathbf{r}) \rho(\mathbf{r}) d\mathbf{r} \right]
$$

where $F[\rho]$ is a universal functional comprising the kinetic energy and electron-electron interactions.

From these two theorems, we can find that
- All observables are functionals of the electron density $\rho(\mathbf{r})$.
- One can reduce the problem from dealing with a many-body wavefunction to a function of three spatial variables.
- The ground-state density $\rho_0(\mathbf{r})$ minimizes the energy functional $E[\rho]$.

## 9.4 The Kohn-Sham Formalism
In 1965, Walter Kohn and Lu Jeu Sham applied Variational Principle to develope a method to make DFT practical by introducing a system of non-interacting electrons that produce the same ground-state density as the interacting system.

Following Kohn-Sham Formalism, the total Energy Functional is expressed as:

$$
E[\rho] = T_s[\rho] + E_{\text{ext}}[\rho] + E_H[\rho] + E_{\text{XC}}[\rho]
$$

- $T_s[\rho]$: Kinetic energy of non-interacting electrons.
- $E_{\text{ext}}[\rho]$: Interaction with external potential.
- $E_H[\rho]$: Hartree energy (classical electron-electron repulsion).
- $E_{\text{XC}}[\rho]$: Exchange-correlation energy, encompassing all many-body effects.

Thus, the Schrödinger equation (now called **Kohn-Sham Equations**) becomes:

$$
\left[ -\frac{\hbar^2}{2m} \nabla^2 + V_{\text{eff}}(\mathbf{r}) \right] \phi_i(\mathbf{r}) = \epsilon_i \phi_i(\mathbf{r})
$$

- $\phi_i(\mathbf{r})$: Kohn-Sham orbitals.
- $\epsilon_i$: Orbital energies.
- $V_{\text{eff}}(\mathbf{r})$: Effective potential.

**Effective Potential** is

$$
V_{\text{eff}}(\mathbf{r}) = V_{\text{ext}}(\mathbf{r}) + V_H(\mathbf{r}) + V_{\text{XC}}(\mathbf{r})
$$

- $V_H(\mathbf{r})$ : Hartree potential.
- $V_{\text{XC}}(\mathbf{r}) = \frac{\delta E_{\text{XC}}[\rho]}{\delta \rho(\mathbf{r})}$ .

The **Electron Density** derived from the solution of KS equation:

$$
\rho(\mathbf{r}) = \sum_{i}^{\text{occ}} |\phi_i(\mathbf{r})|^2
$$

To solve the KS equation, **Exchange-Correlation Functional** is a tricky term that represents the difference between the true kinetic and electron-electron interaction energies and those of the non-interacting reference system. The exact form is unknown. Hence approximations are necessary. Some commond choices are

- ``Local Density Approximation (LDA)`` assumes the exchange-correlation energy at a point depends only on the local density.

$$
E_{\text{XC}}^{\text{LDA}}[\rho] = \int \rho(\mathbf{r}) \varepsilon_{\text{XC}}(\rho(\mathbf{r})) d\mathbf{r}
$$


- ``Generalized Gradient Approximation (GGA)`` is a more advanced version by including density gradients to account for inhomogeneities. This can provide better accuracy for molecular and surface systems.

$$
E_{\text{XC}}^{\text{GGA}}[\rho] = \int f(\rho(\mathbf{r}), \nabla \rho(\mathbf{r})) d\mathbf{r}
$$

Since the electron density becomes the only concern of interest, one can start with a initial guess of the wavefunction and then solve the KS equations iteratively. 

## 9.5 Practical solutions of KS equations explained 

Now, we proceed to apply DFT to solve an One-Dimensional Harmonic Oscillator. Here we start by define the grid.
We aim to find the wavefunction of an One-Dimensional Harmonic Oscillator. 

1. **Define a grid** to describe the wavefunction spanning between at $x_{\text{min}}$ and $x_{\text{max}}$.

```python
import numpy as np
import matplotlib.pyplot as plt

# Spatial grid parameters
x_min, x_max = -5.0, 5.0
N = 1000  # Number of grid points
x = np.linspace(x_min, x_max, N)
dx = x[1] - x[0]
```

2. **Express the Kinetic energy** ($T$) and **External Potential** $V_{\text{ext}}(x)$ is simply a Harmonic potential:
```python
T = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2
V_ext = 0.5 * x**2
```

3. **An initial guess** of Electron Density as a function of $r$,

```python
rho = np.ones(N) * 1e-3  # Small initial density
```

4. **Iterative Solution**
```python
tolerance = 1e-6
max_iterations = 100
for iteration in range(max_iterations):
    # Hartree potential
    V_H = np.zeros(N)
    for i in range(N):
        V_H[i] = np.sum(rho * np.abs(x[i] - x) * dx)

    # Exchange-correlation potential (simplified)
    c = 1  # Constant for exchange
    V_XC = -c / rho

    # Effective potential
    V_eff = V_ext + V_H + V_XC

    # Total Hamiltonian
    H = -0.5 * T + np.diag(V_eff)

    # Solve eigenvalue problem
    energies, wavefunctions = np.linalg.eigh(H)

    # Update electron density
    rho_new = np.abs(wavefunctions[:, 0])**2  # Ground state

    # Check for convergence
    if np.linalg.norm(rho_new - rho) < tolerance:
        print(f'Converged after {iteration+1} iterations')
        break

    rho = rho_new
else:
    print('Did not converge')
```

5. Post-analysis
```python
plt.plot(x, rho)
plt.title('Electron Density')
plt.xlabel('Position x')
plt.ylabel('Density ρ(x)')

plt.plot(x, V_eff, label='Effective Potential')
plt.plot(x, wavefunctions[:, 0], label='Ground State Wavefunction')
plt.legend()

plt.show()
```

### 9.6 Notes on kinetic energy term ($\hat{T}$)

In the Kohn-Sham formalism of DFT, the kinetic energy operator for an electron in one dimension is given by the Laplacian operator, which is the second derivative of the wavefunction with respect to position. The corresponding term in the Hamiltonian is:

$$
\hat{T} = -\frac{\hbar^2}{2m} \frac{d^2}{dx^2}
$$

To solve this equation numerically, we approximate the second derivative using the finite difference method. Consider a function $f(x)$ on a discrete grid of points $x_1, x_2, …, x_N$ with spacing $dx$. The second derivative of the function $f$  at point $x_i$ can be approximated by the finite difference formula:

$$
\frac{d^2 f}{dx^2} \bigg|{x_i} \approx \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{dx^2}
$$

The expression tells us how to approximate the curvature of the function at each point using values of the function at neighboring points.

To apply this approximation to a system with $N$ discrete points, we can represent it using a matrix that operates on the values of the function at all grid points. The matrix will essentially encode the coefficients ``(-2, +1, +1)`` that multiply the values $f(x_i), f(x_{i-1}), f(x_{i+1})$.

- ``np.eye(N)``: an identity matrix of $N \times N$.
- ``np.eye(N, k=1)}``: a matrix of $N \times N$ with ones on the first upper diagonal.
- ``np.eye(N, k=-1)`` : a matrix of $N \times N$  with ones on the first lower diagonal.

Thus, the combination of these matrices gives us a matrix representation of the finite difference approximation of the second derivative operator across the entire grid. The division by $dx^2$ normalizes the matrix to account for the grid spacing.

For a 5-point grid, the matrix $T$ would look something like this:

$$
T = \frac{1}{dx^2}
\begin{bmatrix}
-2 & 1  & 0  & 0  & 0  \\
1  & -2 & 1  & 0  & 0  \\
0  & 1  & -2 & 1  & 0  \\
0  & 0  & 1  & -2 & 1  \\
0  & 0  & 0  & 1  & -2 \\
\end{bmatrix}
$$

This matrix will act on a vector representing the function $f$ at discrete grid points to approximate the second derivative for the entire grid. It physically represents the kinetic energy operator in one dimension (ignoring constants like $\hbar^2/2m$ ). When this matrix acts on the wavefunction, it computes the kinetic energy for the system in a discretized space.
Extensions


## 9.7 Extensions to more challenging systems
- Different Potentials: (e.g., a double-well potential)
- Modify the code to include additional occupied orbitals.
- Implement GGA or other exchange-correlation approximations.

