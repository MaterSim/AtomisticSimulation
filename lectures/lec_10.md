# Week 10. Band Structure Analysis

## 10.1 Formation of Energy Bands: From Molecule to Crystal
The electronic band structure describes the energy levels that electrons can occupy in a solid. Unlike in isolated atoms, where electrons occupy discrete energy levels, in a crystal the periodic potential causes these levels to broaden into bands.

- **Valence Band**: The highest energy band that is completely filled with electrons at 0 K.
- **Conduction Band**: The lowest energy band that is partially filled or completely empty at 0 K.
- **Band Gap**: The energy difference between the top of the valence band and the bottom of the conduction band.

A conductor has overlapping bands (or a partially filled conduction band), while a semiconductor or insulator has a band gap separating the valence and conduction bands.

How are the bands created from the very beginning. You can think about [the process of carbon atoms being brought together to form a diamond crystal](https://en.wikipedia.org/wiki/Band_gap). 
- The right graph shows the energy levels as a function of the spacing between atoms. When far apart (right side of graph) all the atoms have discrete valence orbitals $p$ and $s$ with the same energies.
- When the atoms come closer (left side), their electron orbitals begin to spatially overlap and hybridize into $N$ molecular orbitals each with a different energy, where $N$ is the number of atoms in the crystal. Since $N$ is such a large number, adjacent orbitals are extremely close together in energy so the orbitals can be considered a continuous energy band.
- When the carbon atoms get closer and closer to form a diamond crystal cell, multiple bands are formed, called the valence and conduction bands. In the case of diamond, the highest valence band and the lowest conduction band are separated by a 5.5 eV band gap. The Pauli exclusion principle limits the number of electrons in a single orbital to two, and the bands are filled beginning with the lowest energy.

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/ef/Solid_state_electronic_band_structure.svg/880px-Solid_state_electronic_band_structure.svg.png" alt="Alt text" width="600"/>
</p>

## 10.2 A Quantitative Tight-Binding model 
Although the above explaination provides an intuitive picture about the formation of bands. In this lecture, we plan to be more analytical about this process. Specifically, we will use the tight-binding model to compute the band structure of graphene step by step in Python.

The tight-binding model is a simple yet powerful method for understanding the electronic band structure of materials. It’s particularly useful for systems like graphene, where the electrons are tightly bound to atoms but can still hop between neighboring atomic sites. In graphene, the tight-binding model provides a good approximation for describing the π-bands, which arise from the $p_z$ orbitals.

## 10.3 Application of TB model on Graphene
Graphene’s band structure can be derived from a simple nearest-neighbor tight-binding model. We will focus on the π-bands, which are formed by the $p_z$ orbitals of the carbon atoms.

### 10.3.1 Graphene Lattice and Hamiltonian Setup

- Lattice vectors: Graphene’s lattice is hexagonal with two basis atoms (A and B) per unit cell.
- Nearest-neighbor hopping: Electrons can hop between neighboring A and B atoms with a hopping parameter $t$.

The tight-binding Hamiltonian for graphene can be written as:

$$
H = -t \sum_{\langle i,j \rangle} \left( a_i^\dagger b_j + b_j^\dagger a_i \right)
$$

where $t$ is the hopping energy between neighboring sites (typically around 2.7 eV for graphene), and $a_i^\dagger$ and $b_j^\dagger$ are the creation operators for sublattice A and B, respectively.

We begin by defining the graphene lattice, with lattice vectors $\mathbf{a}_1$  and  $\mathbf{a}_2$, and the three nearest-neighbor vectors  $\delta_1$,  $\delta_2$, and  $\delta_3$.

```python
import numpy as np
import matplotlib.pyplot as plt

# Define lattice vectors for graphene
a = 1.42  # Carbon-carbon bond length in Angstroms

# Lattice vectors
a1 = np.array([np.sqrt(3) * a, 0])
a2 = np.array([np.sqrt(3) / 2 * a, 3 * a / 2])

# Nearest-neighbor vectors (bond vectors)
delta1 = np.array([0, a])
delta2 = np.array([np.sqrt(3)/2 * a, -a / 2])
delta3 = np.array([-np.sqrt(3)/2 * a, -a / 2])

# All nearest-neighbor vectors
deltas = [delta1, delta2, delta3]

# Plot the lattice
plt.figure(figsize=(6, 6))
for i in range(-2, 3):
    for j in range(-2, 3):
        R = i * a1 + j * a2
        plt.plot(R[0], R[1], 'ko')  # Carbon atoms
        for delta in deltas:
            plt.plot([R[0], R[0] + delta[0]], [R[1], R[1] + delta[1]], 'k-')
plt.title('Graphene lattice')
plt.gca().set_aspect('equal')
plt.show()
```

This code defines the graphene lattice and plots it. Each carbon atom is connected to three neighboring carbon atoms, forming the characteristic honeycomb structure.


### 10.3.2 Tight-Binding Hamiltonian

We now construct the tight-binding Hamiltonian for graphene, taking into account the nearest-neighbor hopping between sublattices A and B.

The Hamiltonian in k-space is written as:

$$
H(\mathbf{k}) =
\begin{pmatrix}
0 & f(\mathbf{k}) \\
f^*(\mathbf{k}) & 0
\end{pmatrix}
$$

where $f(\mathbf{k}) = -t \left( e^{i \mathbf{k} \cdot \delta_1} + e^{i \mathbf{k} \cdot \delta_2} + e^{i \mathbf{k} \cdot \delta_3} \right)$.

The eigenvalues of this matrix give us the band energies for the π-bands.

3. Define the Hamiltonian in k-Space

Let’s define the function $f(\mathbf{k})$ that represents the hopping terms in the Hamiltonian.

```python
def f_k(kx, ky, a, t):
    """Function f(k) representing hopping terms in k-space."""
    delta1 = np.array([0, a])
    delta2 = np.array([np.sqrt(3)/2 * a, -a / 2])
    delta3 = np.array([-np.sqrt(3)/2 * a, -a / 2])
    return -t * (np.exp(1j * (kx * delta1[0] + ky * delta1[1])) +
                 np.exp(1j * (kx * delta2[0] + ky * delta2[1])) +
                 np.exp(1j * (kx * delta3[0] + ky * delta3[1])))

# Define the band structure computation
def graphene_band_structure(kx, ky, a, t):
    """Compute the band structure of graphene using the tight-binding model."""
    f = f_k(kx, ky, a, t)
    H_k = np.array([[0, f], [np.conj(f), 0]])  # Tight-binding Hamiltonian in k-space
    eigenvalues = np.linalg.eigvalsh(H_k)
    return eigenvalues
```

4. Compute Band Structure Along High-Symmetry Directions

Now we compute the band structure along high-symmetry points in the Brillouin zone. The typical path is $\Gamma \rightarrow K \rightarrow M \rightarrow \Gamma$.

```python
# High-symmetry points in k-space
K_point = [4*np.pi/(3*np.sqrt(3)*a), 0]
M_point = [np.pi/(np.sqrt(3)*a), np.pi/a]
Gamma_point = [0, 0]

# Interpolation between high-symmetry points
def interpolate_points(p1, p2, n):
    return np.linspace(p1, p2, n)

# Parameters for graphene
t = 2.7  # Hopping energy in eV

# Compute the band structure along high-symmetry directions
n_points = 100
kx_values = []
ky_values = []
kx_values += interpolate_points(Gamma_point[0], K_point[0], n_points).tolist()
ky_values += interpolate_points(Gamma_point[1], K_point[1], n_points).tolist()
kx_values += interpolate_points(K_point[0], M_point[0], n_points).tolist()
ky_values += interpolate_points(K_point[1], M_point[1], n_points).tolist()
kx_values += interpolate_points(M_point[0], Gamma_point[0], n_points).tolist()
ky_values += interpolate_points(M_point[1], Gamma_point[1], n_points).tolist()

# Calculate the energy bands
bands = []
for kx, ky in zip(kx_values, ky_values):
    bands.append(graphene_band_structure(kx, ky, a, t))

bands = np.array(bands)

# Plot the band structure
plt.figure(figsize=(8, 6))
plt.plot(bands[:, 0], label='Lower band')
plt.plot(bands[:, 1], label='Upper band')
plt.xticks(ticks=[0, n_points, 2*n_points, 3*n_points],
           labels=['Gamma', 'K', 'M', 'Gamma'])
plt.ylabel('Energy (eV)')
plt.title('Band Structure of Graphene')
plt.legend()
plt.grid(True)
plt.show()
```

## 10.4 Interpretation of the Band Structure

In the plotted band structure of graphene:

- The Dirac points at the  K -point are where the valence and conduction bands meet, and the energy gap is zero.
- Near the Dirac points, the bands show a linear dispersion, which indicates the massless Dirac fermions that are responsible for graphene’s unique electronic properties.
- The bands are symmetric about the zero energy level, representing the bonding and anti-bonding states from the tight-binding model.

Some notable properties include
1. **Linear Dispersion**: The tight-binding model shows a linear energy dispersion near the Dirac points (located at the $K$ and $K{\prime}$ points in the Brillouin zone). This linear relationship between energy and momentum is characteristic of massless Dirac fermions, which are responsible for the high mobility of charge carriers in graphene.
2. **Zero Band Gap**: At the Dirac points, the valence and conduction bands touch, resulting in zero band gap. This makes graphene a semimetal (or more precisely, a zero-gap semiconductor), where the conduction and valence bands meet at the Fermi level.
3. **Bonding and Anti-Bonding States**: The two bands that arise from the tight-binding Hamiltonian correspond to the bonding (lower) and anti-bonding (upper) states. These bands arise due to the interaction between the p_z orbitals on neighboring carbon atoms in the honeycomb lattice.
4. **High Electron Mobility**: The linear dispersion near the Dirac points leads to high electron mobility in graphene. The electrons behave like massless Dirac fermions, which is why graphene exhibits exceptional electrical conductivity.
5. **Optical Properties**: The unique band structure of graphene also influences its optical properties. The zero band gap allows graphene to absorb light across a wide range of frequencies, making it useful for optoelectronic applications.

This lecture could form a foundation for further studies on more complex materials and advanced band structure concepts using tight binding. Similarly, the calculation of band structure can be done in DFT as well, but with a heavier computational cost.
