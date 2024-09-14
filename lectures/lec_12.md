# Week 12. Phonon calculation from both classical force field and DFT

## 12.1 A simple spring

In a simple spring, Hooke’s Law states that the force $F$ exerted by a spring is proportional to the displacement $x$ from its equilibrium position:

$$
F = -k x
$$

According to Newton’s second law, the net force $F$ acting on a mass $m$ causes an acceleration $a$:

$$
F = m a = m \frac{d^2 x}{dt^2}
$$

Therefore, we obtain the following differential equation to describe the motion,

$$
m \frac{d^2 x}{dt^2} = -k x
$$

Obviously, the general solution is a sin function.

$$
x(t) = A \cos(\omega t + \phi)
$$

And $\omega = \sqrt{\frac{k}{m}}$ is the angular frequency.

Thus, the vibration frequency $f$ is related to the force constant $k$ and the mass $m$  as:

$$
f = \frac{1}{2\pi} \sqrt{\frac{k}{m}}
$$


## 12.2 An 1D Infinite Chain of Identical Atoms
Now let's consider a relatively more complex system with a infinite number of identical spring. This model helps in understanding phonons, the quantized vibrations in a lattice, which are critical in thermal and electrical properties of solids.

We’ll consider a chain of identical atoms, each with mass $m$, connected by springs with force constant $k$. In this model, the atoms are spaced a distance $a$ apart. Each atom oscillates around its equilibrium position.

## 12.2.1 Problem setup
Let $u_n(t)$ represent the displacement of the $n$-th atom from its equilibrium position at time $t$. For the $n$-th atom, the total force exerted on it comes from its interactions with its nearest neighbors, i.e., the ($n$+1)-th atom and the ($n$-1)-th atom.

According to Hooke’s Law, the force on the $n$-th atom due to its neighbors is:

$$
F_n = -k \left( u_n(t) - u_{n-1}(t) \right) - k \left( u_n(t) - u_{n+1}(t) \right)
$$

$$
F_n = -k \left( 2 u_n(t) - u_{n+1}(t) - u_{n-1}(t) \right)
$$

Using $F_n = m \frac{d^2 u_n}{dt^2}$, the equation of motion for the $n$-th atom becomes:

$$
m \frac{d^2 u_n}{dt^2} = -k \left( 2 u_n(t) - u_{n+1}(t) - u_{n-1}(t) \right)
$$


### 12.2.2 Solution
Intuitively, we can see that the general solution should behave like some waves

$$
u_n(t) = A e^{i (nqa - \omega t)}
$$

Where:

- $A$ is the amplitude of the wave.
- $q$  is the wavevector, related to the wavelength of the vibration.
- $\omega$  is the angular frequency of the oscillation.
- $a$  is the distance between neighboring atoms.

Plugging this solution into the equation of motion, we get:

$$
m \frac{d^2 u_n}{dt^2} = -m \omega^2 u_n(t)
$$

Substitute this into the equation of motion:

$$
m \omega^2 A e^{i(nqa - \omega t)} = -k \left( 2 A e^{i(nqa - \omega t)} - A e^{i((n+1)qa - \omega t)} - A e^{i((n-1)qa - \omega t)} \right)
$$

Simplifying the right-hand side:

$$
k A \left( 2 - e^{iqa} - e^{-iqa} \right) e^{i(nqa - \omega t)}
$$

Using the identity  e^{iqa} + e^{-iqa} = 2 \cos(qa) , we obtain:

$$
m \omega^2 A e^{i(nqa - \omega t)} = -k A \left( 2 - 2 \cos(qa) \right) e^{i(nqa - \omega t)}
$$

Canceling common terms, we are left with:

$$
m \omega^2 = 2k \left( 1 - \cos(qa) \right)
$$

Hence, the angular frequency $\omega$ is given by:

$$
\omega(q) = \sqrt{\frac{2k}{m} \left( 1 - \cos(qa) \right)}
$$

Using the trigonometric identity $1 - \cos(qa) = 2 \sin^2\left( \frac{qa}{2} \right)$ , we can rewrite the angular frequency as:

$$
\omega(q) = 2 \sqrt{\frac{k}{m}} \left| \sin\left( \frac{qa}{2} \right) \right|
$$

Using $f(q) = \frac{\omega(q)}{2\pi}$, we get $f(q)$

$$
f(q) = \frac{1}{\pi} \sqrt{\frac{k}{m}} \left| \sin\left( \frac{qa}{2} \right) \right|
$$


### 12.2.3 Physical Insights
- Dispersion relation. we find the vibrational frequency depends on the wavevector. 
- Long Wavelength Limit: For small $q$, $\sin(q a / 2) \approx q a / 2$, so the frequency becomes approximately linear in $q$, i.e.,  $\omega(q) \propto q$, which is typical for acoustic phonons.
- At the edge of the Brillouin zone, where $q = \frac{\pi}{a}$, the frequency reaches its maximum value:

```python
import numpy as np
import matplotlib.pyplot as plt

# Parameters
k = 1.0  # Force constant (N/m)
m = 1.0  # Mass of each atom (kg)
a = 1.0  # Lattice constant (m)

# Function to compute angular frequency omega(q)
def omega_q(q, k, m, a):
    return 2 * np.sqrt(k / m) * np.abs(np.sin(q * a / 2))

# Generate q-values in the first Brillouin zone (-pi/a to pi/a)
q_values = np.linspace(-np.pi/a, np.pi/a, 100)

# Compute the corresponding frequencies
omega_values = omega_q(q_values, k, m, a)

# Plot the vibrational frequency vs. wavevector q
plt.figure(figsize=(8, 6))
plt.plot(q_values, omega_values, label=r'$\omega(q)$')
plt.xlabel(r'Wavevector $q$ (1/m)')
plt.ylabel(r'Angular frequency $\omega(q)$ (rad/s)')
plt.title('Vibrational Frequency vs. Wavevector in a 1D Chain')
plt.axhline(0, color='black', lw=0.5)
plt.axvline(0, color='black', lw=0.5)
plt.legend()
plt.grid(True)
plt.show()
```

## 12.3 An 1D diatomic chain model
In the 1D diatomic chain model, we consider two different types of atoms (A and B) with masses $m_A$ and $m_B$, alternating along the chain. The atoms are connected by springs with a spring constant $k$.


## 12.3.1 Equation of Motions
The displacement of atom A in the $n$-th unit cell is denoted by $u_n^A(t)$, and the displacement of atom B is  $u_n^B(t)$. And the forces on each atom come from the interactions with neighboring atoms. Using Hooke’s law and Newton’s second law, the equations of motion for atoms A and B are:

$$
m_A \frac{d^2 u_n^A}{dt^2} = -k \left( u_n^A - u_n^B \right) - k \left( u_n^A - u_{n-1}^B \right)
$$

$$
m_B \frac{d^2 u_n^B}{dt^2} = -k \left( u_n^B - u_n^A \right) - k \left( u_n^B - u_{n+1}^A \right)
$$

## 12.3.2 Solutions
We assume wave-like solutions for the displacements of atoms A and B:

$$
u_n^A(t) = A e^{i(nqa - \omega t)}
$$

$$
u_n^B(t) = B e^{i(nqa - \omega t)}
$$

Substituting these into the equations of motion, we get two coupled equations:

$$
-m_A \omega^2 A = -k A \left( 2 - e^{-iqa} - e^{iqa} \right) + k B \left( 1 + e^{iqa} \right)
$$

$$
-m_B \omega^2 B = -k B \left( 2 - e^{-iqa} - e^{iqa} \right) + k A \left( 1 + e^{-iqa} \right)
$$

Using the identity  e^{iqa} + e^{-iqa} = 2 \cos(qa) , we simplify these to:

$$
-m_A \omega^2 A = -2k A \left( 1 - \cos(qa) \right) + k B \left( 2 \cos(qa) \right)
$$

$$
-m_B \omega^2 B = -2k B \left( 1 - \cos(qa) \right) + k A \left( 2 \cos(qa) \right)
$$

Substituting these into the equations of motion, we get a system of equations that can be written as a matrix equation:

$$
\begin{pmatrix}
m_A \omega^2 & 0 \\
0 & m_B \omega^2
\end{pmatrix}
\begin{pmatrix}
A \\
B
\end{pmatrix}
=k
\begin{pmatrix}
2 - 2 \cos(qa) & -2 \cos(qa) \\
-2 \cos(qa) & 2 - 2 \cos(qa)
\end{pmatrix}
\begin{pmatrix}
A \\
B
\end{pmatrix}
$$

Moving the mass terms from left to right, we define the dynamical matrix $D(q)$, which describes the interaction between atoms A and B in the lattice.

$$
D(q) =
\frac{k}{m_A m_B}
\begin{pmatrix}
m_B(2 - 2 \cos(qa)) & -2m_B \cos(qa) \\
-2m_A \cos(qa) & m_A(2 - 2 \cos(qa))
\end{pmatrix}
$$

To find the phonon frequencies, we solve the **eigenvalue problem**:

$$
D(q) \begin{pmatrix} A \\ B \end{pmatrix} = \omega^2 \begin{pmatrix} A \\ B \end{pmatrix}
$$

Solving this determinant equation gives the dispersion relations for the acoustic and optical phonon modes, 

- For the acoustic mode, the atoms move in phase, and the frequency $\omega(q)$ is small at small $q$. The dispersion relation is approximately linear for small  q :

$$
\omega_{\text{acoustic}}(q) \approx \frac{c_s}{a} |q|
$$

Where  $c_s$  is the speed of sound in the material.

- For the optical mode, the atoms move out of phase, and the frequency  $\omega(q)$  is larger. The optical phonon mode has a higher frequency because the atoms are oscillating against each other.


```python
import numpy as np
import matplotlib.pyplot as plt

# Parameters
k = 1.0  # Force constant (N/m)
m_A = 1.0  # Mass of atom A (kg)
m_B = 2.0  # Mass of atom B (kg)
a = 1.0   # Lattice constant (m)

# Dispersion relation for the acoustic and optical phonon modes
def phonon_dispersion(q, k, m_A, m_B, a):
    term = 2 * k * (1 - np.cos(q * a))
    omega_acoustic = np.sqrt(term / (m_A + m_B))
    omega_optical = np.sqrt(term * (m_A + m_B) / (m_A * m_B))
    return omega_acoustic, omega_optical

# Generate q-values in the first Brillouin zone (-pi/a to pi/a)
q_values = np.linspace(-np.pi/a, np.pi/a, 100)

# Compute the corresponding phonon frequencies
omega_acoustic_values = []
omega_optical_values = []
for q in q_values:
    omega_acoustic, omega_optical = phonon_dispersion(q, k, m_A, m_B, a)
    omega_acoustic_values.append(omega_acoustic)
    omega_optical_values.append(omega_optical)

# Plot the acoustic and optical phonon dispersion relations
plt.figure(figsize=(8, 6))
plt.plot(q_values, omega_acoustic_values, label='Acoustic Mode')
plt.plot(q_values, omega_optical_values, label='Optical Mode')
plt.xlabel(r'Wavevector $q$ (1/m)')
plt.ylabel(r'Angular frequency $\omega(q)$ (rad/s)')
plt.title('Phonon Dispersion Relations in 1D Diatomic Chain')
plt.axhline(0, color='black', lw=0.5)
plt.axvline(0, color='black', lw=0.5)
plt.legend()
plt.grid(True)
plt.show()
```


## 12.3.2 Physical Insights
The 1D diatomic chain model gives rise to both acoustic and optical phonon modes. In the acoustic mode, atoms oscillate in phase, leading to low-frequency vibrations, while in the optical mode, atoms oscillate out of phase, leading to higher-frequency vibrations. These two phonon branches are a direct consequence of having two different types of atoms in the unit cell, and they play an essential role in the thermal and optical properties of materials.

The dynamical matrix is a key concept in lattice dynamics and phonon theory. It arises in the context of solving the equations of motion for atoms in a crystal lattice, especially when dealing with harmonic vibrations in the lattice, such as in the 1D chain model with two types of atoms (diatomic chain). The dynamical matrix relates the forces on atoms to their displacements and encapsulates the vibrational properties of the system.

## 12.4 Extension to realistic systems

To compute phonon dispersions for a realistic 3D crystal, we need to extend the concepts from the 1D chain to a 3D lattice, which involves considering all the interactions between atoms in the crystal unit cell. In a 3D crystal, the phonon dispersion relation tells us how the phonon frequencies (or energies) depend on the wavevector $\mathbf{q}$ in different directions of the Brillouin zone.

If there are $N$ atoms in the unit cell, each atom has three degrees of freedom (x, y, z). Therefore, the system has $3N$ degrees of freedom, leading to $3N$ phonon modes at each wavevector $\mathbf{q}$. These include:

- 3 acoustic phonon modes: Low-frequency vibrations where the entire unit cell moves in phase.
- 3$N$ - 3 optical phonon modes: Higher-frequency modes where the atoms in the unit cell vibrate relative to each other.

The dynamical matrix in 3D describes the interaction between all the atoms in the unit cell, taking into account their positions and the interatomic forces. For a crystal with $N$ atoms per unit cell, the dynamical matrix  $D(\mathbf{q})$  is a  $3N \times 3N$  matrix.

The dynamical matrix is computed from the force constants $\Phi$ between atoms, which describe the second derivatives of the potential energy with respect to atomic displacements. The matrix elements of the dynamical matrix are:

$$
D_{\alpha\beta}^{ij}(\mathbf{q}) = \frac{1}{\sqrt{m_i m_j}} \sum_{\mathbf{R}} \Phi_{\alpha\beta}^{ij}(\mathbf{R}) e^{i \mathbf{q} \cdot \mathbf{R}}
$$

Where:

- $i$, $j$  index the atoms in the unit cell.
- $\alpha$, $\beta$  represent the Cartesian components $(x, y, z)$.
- $m_i$, $m_j$  are the masses of atoms $i$ and $j$.
- $\Phi_{\alpha\beta}^{ij}(\mathbf{R})$  are the force constants between atoms $i$ and $j$ in unit cells separated by the vector $\mathbf{R}$.
- $\mathbf{q}$ is the wavevector.

This matrix must be diagonalized to obtain the phonon frequencies  \omega(\mathbf{q})  for each mode at each wavevector.

Below are the steps to Compute Phonon in a 3D Crystal.


1. Get Atomic Positions: Obtain the positions of atoms in the unit cell and lattice vectors that define the crystal structure.

2. Compute Force constants: These describe the interaction between atoms and can either be obtained from ab initio calculations (using DFT) or from experimental data. The force constants can be represented as a matrix that relates atomic displacements to forces.

3. Construct the Dynamical Matrix.  For each wavevector $\mathbf{q}$ in the Brillouin zone, construct the dynamical matrix  $D(\mathbf{q})$  using the force constants. The dynamical matrix will be a $3N \times 3N$  matrix, where $N$ is the number of atoms in the unit cell.

4. Diagonalize the Dynamical Matrix. For each wavevector $\mathbf{q}$, diagonalize the dynamical matrix $D(\mathbf{q})$. The eigenvalues of the matrix give the squared frequencies $\omega^2(\mathbf{q})$, and the eigenvectors describe the polarization vectors (atomic displacements) for the phonon modes.

5. Compute the Phonon Dispersion Relation and density of states. Once the phonon frequencies  $\omega(\mathbf{q})$  are computed, plot them as a function of $\mathbf{q}$ along high-symmetry directions in the Brillouin zone.

6. Optionally, compute the phonon density of states. To compute the phonon DOS, you must sample the Brillouin zone using a k-point grid or a q-point grid. The finer the grid, the more accurate your calculation will be. The total number of wavevectors sampled is denoted as $N_q$.

$$
g(\omega) = \frac{1}{N_q} \sum_{\mathbf{q}, j} \delta(\omega - \omega_j(\mathbf{q}))
$$

Where:

- $\omega_j(\mathbf{q})$  are the phonon frequencies for the  $j$-th mode at wavevector $\mathbf{q}$.
- $N_q$  is the total number of sampled wavevectors in the Brillouin zone.
- The delta function  $\delta(\omega - \omega_j(\mathbf{q}))$  ensures that only the phonon states with frequencies near  $\omega$  contribute to the DOS at $\omega$.
