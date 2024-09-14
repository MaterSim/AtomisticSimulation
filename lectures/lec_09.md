# Week 9. DFT Simulation of the Periodic System

DFT has proven to be one of the most reliable methods for predicting material properties such as band structure, electronic density, and ground-state energy for crystals.

However, it is not straightforward to extend DFT from molecule to crystals. Unlike molecules, crystals are periodic and extend infinitely. We need methods that respect this periodicity and allow us to work with finite-sized models.

Silicon (Si) is a semiconductor and an important material for electronics. DFT can help us study its electronic band structure, density of states, and structural properties like lattice constants.

## 9.1 Bloch’s Theorem:
Bloch’s theorem simplifies the problem by showing that the wavefunctions in periodic potentials can be written as plane waves modulated by a periodic function. This allows us to work in k-space and reduce the complexity.

According to Bloch’s Theorem, the wavefunctions of an electron in a periodic potential take the form:

$$
\psi_{\mathbf{k}}(\mathbf{r}) = e^{i \mathbf{k} \cdot \mathbf{r}} u_{\mathbf{k}}(\mathbf{r})
$$

Where $\mathbf{k}$ is the wavevector and $u_{\mathbf{k}}(\mathbf{r})$ is a periodic function with the periodicity of the lattice.

- The periodic nature of $u_{\mathbf{k}}(\mathbf{r})$ means we only need to solve the problem within a unit cell.
- Reciprocal Lattice:
- For crystals, it is useful to define the reciprocal lattice, which is the Fourier transform of the real-space lattice. The reciprocal lattice defines k-points that represent wavevectors in reciprocal space.
- The Brillouin zone is the primitive cell in reciprocal space. In the case of silicon, the first Brillouin zone is particularly important for computing the electronic structure.

## 9.2 KS DFT for Crystals

The Kohn-Sham equations for periodic solids are solved in reciprocal space (k-space). The Hamiltonian in k-space is periodic and can be solved using plane-wave basis sets.

In periodic solids, it is convenient to expand the wavefunction  \psi_{\mathbf{k}}(\mathbf{r})  as a sum of plane waves:

$$
\psi_{\mathbf{k}}(\mathbf{r}) = \sum_{\mathbf{G}} C_{\mathbf{k}, \mathbf{G}} e^{i (\mathbf{k} + \mathbf{G}) \cdot \mathbf{r}}
$$

where $\mathbf{G}$ are the reciprocal lattice vectors and $C_{\mathbf{k}, \mathbf{G}}$ are the expansion coefficients.

The accuracy of the calculation depends on the plane-wave energy cutoff, which determines how many plane waves are included in the basis set.

## 9.3 Brillouin Zone Sampling

When studying crystalline solids in DFT, we need to account for the periodic nature of the crystal lattice. The behavior of electrons in a crystal is described in terms of Bloch states, which means that their wavefunctions depend on the wavevector $\mathbf{k}$ in reciprocal space (also known as k-space). The reciprocal space of a crystal is divided into regions known as Brillouin zones, and solving the DFT equations over the entire Brillouin zone is crucial to accurately capture the electronic properties of the system.

In periodic solids, the electronic properties (such as energy levels and densities) depend on the wavevector  \mathbf{k} . Since the wavevector can take continuous values within the Brillouin zone, we cannot solve the Kohn-Sham equations for every possible  \mathbf{k} -point. Instead, we need to sample the Brillouin zone at a discrete set of points and integrate over the zone to compute quantities like the total energy, charge density, and density of states.

Brillouin zone sampling is typically done by choosing a grid of k-points that represent the possible electronic states within the zone. The accuracy of the DFT calculation depends on how well the Brillouin zone is sampled. The more k-points you use, the more accurate your results will be, but it will also increase the computational cost.

```python
def generate_k_points(N_k):
    """Generate a Monkhorst-Pack grid for Brillouin zone sampling."""
    k_points = []
    for i in range(N_k):
        for j in range(N_k):
            for k in range(N_k):
                kx = (i - (N_k // 2)) / N_k
                ky = (j - (N_k // 2)) / N_k
                kz = (k - (N_k // 2)) / N_k
                k_points.append([kx, ky, kz])
    return np.array(k_points)

k_points = generate_k_points(N_k=4)  # 4x4x4 k-point grid
```

## 9.4 Pseudopotentials 
In crystalline materials, there are many electrons, but only the valence electrons (those in the outer shells) play a significant role in chemical bonding and material properties. The core electrons (those closer to the nucleus) are tightly bound and do not change much between different environments. To simplify the problem, DFT often uses pseudopotentials to represent the effect of the core electrons and the nucleus, so that only the valence electrons need to be explicitly treated.

A pseudopotential is a simplified model that replaces the full Coulomb potential of the nucleus and the core electrons with a smoother, effective potential. The pseudopotential is constructed so that:

- The valence electrons feel the correct effective potential due to the nucleus and core electrons.
- The wavefunctions of the valence electrons are smooth and can be described by a plane-wave basis.

In this way, the pseudopotential approximates the effects of the core electrons without explicitly solving for them. This drastically reduces the computational cost, as fewer electrons need to be handled, and the smoother wavefunctions of the valence electrons require fewer plane waves for an accurate representation.

Pseudopotentials are created by:

- Solving the all-electron Schrödinger equation for an isolated atom.
- Removing the core electron states and keeping only the valence states.
- Replacing the Coulomb potential with a smoother potential that reproduces the valence wavefunctions in the regions of interest.

The result is a pseudopotential that accurately mimics the interaction between the valence electrons and the atomic nucleus + core electrons without having to model the core electrons explicitly.

## 9.5 Self-Consistent Field (SCF) Loop

We implement the SCF loop to solve the Kohn-Sham equations iteratively. For each k-point:
1.	Solve the Kohn-Sham equations using a plane-wave basis.
2.	Update the electron density.
3.	Recompute the effective potential.
4.	Repeat until convergence is achieved.

