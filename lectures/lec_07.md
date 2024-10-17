# 7 Representation of Local Atomic Environment

In atomic and molecular systems, the local environment around a particle plays a crucial role in determining its physical properties. Understanding how atoms or molecules are arranged in space provides insights into the material’s structural characteristics, phase transitions, and dynamic behaviors. The local atomic environment can be described through various mathematical tools. In this lecture, we will birefly start with the previously mentioned radial distribution function. While RDF is effective at describing radial distributions, it falls short when it comes to capturing angular information, which is critical in systems with orientational order. Hence, we attempt to explore the descriptors that can deal with both radial and angular information with mathematical rigor.

## 7.1 Radial Distribution Function and its Limiation

Radial distribution function $g(r)$, is a fundamental tool for describing the local structure in a system of particles. It gives the probability of finding a particle at a distance $r$ from a reference particle, normalized by the average particle density,

$$
g(r) = \frac{1}{\rho N} \left\langle \sum_{i} \sum_{j \neq i} \delta(r - r_{ij}) \right\rangle
$$

Where:

- $\rho$ is the number density of the system.
- $r_{ij}$ is the distance between particle $i$ and particle $j$.
- $\delta$ is the Dirac delta function, ensuring that contributions are only made when the distance between particles matches $r$.

While the RDF is a powerful tool for capturing the radial distribution of particles, it is unable to describe angular correlations between particles. As a result, systems with orientational order, such as liquid crystals or crystals with complex angular symmetries, cannot be fully characterized by the pair distribution function alone.

## 7.2 Orientational Order Parameter

To capture the missing angular information in systems with orientational order, we need to introduce additional descriptors, such as the orientational order parameter. This parameter quantifies the degree of angular ordering among neighboring particles. Unlike RDF which focuses solely on the radial distances, orientational order parameters provide insights into how the bonds between particles are aligned in space. These parameters are particularly useful in distinguishing between phases with different degrees of symmetry, such as solid, liquid, or nematic phases.

### 7.2.1 Bond Orientational Order Parameter in a 2D system
Let’s first consider a two-dimensional system, where the angular relationships between a particle and its nearest neighbors can provide critical information about the system’s structure. To quantify these angular relationships, we define the bond orientational order parameter $\psi_m$, which measures how the bonds between a particle and its neighbors are aligned with respect to a reference axis. The order parameter $\psi_m$ is given by the following equation:

$$
\psi_m = \frac{1}{N} \sum_{j=1}^{N} e^{i m \theta_j}
$$

Where:
- $N$ is the number of neighbors around a given particle.
- $\theta_j$ is the angle formed between the bond connecting a particle to its neighbor $j$ and some reference axis.
- $m$ is the symmetry index. For example, $m$ = 6 is used for systems with hexagonal symmetry.

The value of $\psi_m$ is a **complex number**, and it will depend on the degree of angular ordering in the system. 

- **Magnitude**: The magnitude $|\psi_m(i)|$ measures how well the local arrangement of neighbors conforms to a m-fold symmetric structure. If the neighbors are perfectly arranged in a hexagonal pattern, $|\psi_6(i)|$ will be close to 1. In disordered regions, the value of $|\psi_6(i)|$ will be closer to 0.
- **Phase**: The phase (angle of $\psi_m(i)$ ) indicates the orientation of the bond order relative to the reference direction.

Thus, by examining how  $\psi_m$ evolves with temperature, pressure, or density, researchers can gain insights into the structural transformations occurring in the system.

### 7.2.2 Extension to 3D: Neighbor density function

In real-world scenarios, we are primarily dealing with 3D systems. How can we extend the approach we discussed for 2D to 3D? This extension involves additional complexity because, with the introduction of an extra dimension, the alignment of atoms can no longer be described by a single variable, as in 2D. To capture this alignment in 3D, we need to be more rigorous with our mathematical description.

When discussing the 2D Bond Orientational Order Parameter, we derived a formula that aimed to capture the angular relationships between particles. Essentially, we were measuring how the particles around a reference atom were arranged, particularly focusing on their angular orientations. In 3D, this concept becomes more involved. Specifically, we want to describe how atoms are arranged around a central reference atom, taking into account both their radial and angular distributions.

To begin with, we define the **Atomic Neighbor Density Function** to describe the spatial distribution of atoms around a reference atom within a cutoff radius  r_c . The atomic neighbor density function is expressed as:

$$
    \rho(\mathbf{r}) = \sum_i^{r_i \leq rc} \delta(\mathbf{r}-\mathbf{r_i})
$$

Here, $\delta(\mathbf{r} - \mathbf{r_i})$ is the Dirac delta function, which ensures that the function only contributes when a neighboring atom is located at  $\mathbf{r_i}$, and the summation runs over all neighboring atoms within the cutoff radius $rc$.

To capture the angular distribution of neighboring atoms, we can transform the spatial neighbor density function into another domain, similar to how the Fourier transform converts a time-domain signal into its frequency components. In this case, we are interested in projecting the atomic density distribution onto the unit sphere to study the angular arrangement of atoms.

The popular choice for the basis functions on the unit sphere is [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics), which are functions defined on the surface of a sphere. Therefore, we can expand the neighbor density function $\rho(\mathbf{r})$ as a series of spherical harmonics on the 2-sphere:

$$
    \rho(\mathbf{r}) = \sum_{l=0}^{+\infty} \sum_{m=-l}^{+l} c_{lm} Y_{lm} (\mathbf{\hat{r}})
$$

In this expression:

- $Y_{lm}(\mathbf{\hat{r}})$  are the spherical harmonics, which form a complete orthonormal basis on the sphere.
- $\mathbf{\hat{r}}$  is the normalized radial vector, with a unit length of 1.
- $c_{lm}$  are the expansion coefficients, which describe the contribution of each spherical harmonic mode to the overall distribution.

Similar to the Fourier transform, the expansion coefficients $c_{lm}$ are complex numbers that capture the angular characteristics of the neighbor density. These coefficients can be computed by projecting the neighbor density function onto the spherical harmonics:

$$
    c_{lm} = \left< Y_{lm}(\mathbf{\hat{r}})|\rho(\mathbf{r}) \right> = \sum_i^{r_i \leq r_c}Y_{lm}(\mathbf{\hat{r}_i}).
$$

Here, the indices $l$ and $m$ correspond to the angular frequency components, with $l$ denoting the total angular momentum and $m$ representing its projection along a chosen axis. Unlike the 1D Fourier transform, where a single frequency index suffices, spherical harmonics require two indices $l$ and $m$  to describe the angular frequencies in 3D.

These expansion coefficients, $c_{lm}$, contain detailed information about the neighbor density around the reference atom. However, a key issue is that these coefficients are sensitive to rotations of the system. If the system is rotated, the values of $c_{lm}$ will change, which is undesirable when trying to describe the local atomic environment in a way that is independent of orientation.

In practical applications, we aim to find a representation of the local atomic environment that is both real-valued and invariant under translations and rotations of the system.


### 7.2.3 3D Bond Order Parameters

To address this challenge, Steinhardt introduced [bond order parameters](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784) in 1983, which use second- and third-order combinations of the expansion coefficients $c_{lm}$ to quantify the order in liquids and glasses. These bond order parameters are rotationally invariant, making them useful for characterizing local atomic environments without being affected by the orientation of the system.

The bond order parameter $p_l$ is defined as:

$$   
p_l = \sum_{m=-l}^{+l}c_{lm}c^*_{lm}
$$

Note that this was called $Q_l$ in the original paper. However, it was later found that bond order parameter is closely related to the power spectrum. Hence, we will call it $p_l$ from now on. In signal processing, the power spectrum describes how the power of a signal is distributed across different frequency components. Similarly, in this context, $p_l$ measures the **power** of the neighbor density when projected onto the angular frequency components represented by $l$.

In general, the power spectrum $p_l$ is the Fourier transform of the autocorrelation function, and it provides a frequency-domain representation of the dependencies captured by the autocorrelation.

When analyzing the $p_l$ series for different structures, Steinhardt found that $p_4$ and $p_6$ were particularly useful in distinguishing between different crystal structures such as body-centered cubic (bcc), face-centered cubic (fcc), hexagonal close-packed (hcp), and icosahedral arrangements. These parameters have proven to be very useful for analyzing molecular dynamics (MD) simulations, particularly when identifying structural differences between solids, liquids, and glasses. For instance, LAMMPS allows the computation of [bond-orientational order parameters](https://docs.lammps.org/compute_orientorder_atom.html) in several kinds of styles.

**Limitations** of $p_l$ 

- The bond order parameter $p_l$  does not capture any radial information; it is purely an angular descriptor.
- It assumes a neighbor density in the form of a Dirac delta function. This may not be good for the purpose of measuring the similarities between two environments.

Despite these limitations, the idea of using spherical harmonics and power spectra to describe the local atomic environment has inspired many subsequent works and is still widely used in computational materials science today.

## 7.3 Manybody descriptors

In modern computational materials science and atomistic simulations, understanding the local atomic environment goes beyond simple pairwise interactions. Describing how groups of atoms are collectively arranged, including both radial and angular components, is critical for capturing the structural complexity of systems such as liquids, glasses, and complex crystals. Manybody descriptors provide a mathematical framework to represent these arrangements, offering a more detailed picture of the atomic environment than traditional pair distribution functions or angular descriptors alone.

Manybody descriptors, such as the power spectrum and bispectrum, help quantify the relative positions of multiple atoms in a way that is invariant to rotation and translation. 

### 7.3.1 Power Spectrum with the Explicit Radial Component

In practical applications, we often care about how atoms are spatially arranged in both radial and angular space. While previous approaches focused primarily on the angular distribution, the introduction of radial information allows for a more comprehensive description of the local atomic environment.

In their 2012 paper, Bartók et al. introduced an improved manybody descriptor that explicitly incorporates both radial and angular components. This approach overcomes the limitation of describing neighbor density with a Dirac delta function by replacing the delta function with a Gaussian function of limited width $\alpha$. This smoothing allows for a more realistic representation of how atoms are distributed around a reference atom.

The modified neighbor density function is given by:

$$
\begin{equation*}
    \rho'(\mathbf{r}) = \sum_i^{r_i \leq rc} e^{(-\alpha|\mathbf{r}-\mathbf{r_i}|^2)}
                      = \sum_i^{r_i \leq rc} e^{-\alpha(r^2+r_i^2)} e^{2\alpha \mathbf{r} \cdot \mathbf{r_i}} 
\end{equation*}
$$

This expression can be further expanded as:

$$
\rho'(\mathbf{r}) = \sum_i^{r_i \leq rc} \sum_{lm}  4\pi e^{-\alpha(r^2+r_i^2)} I_l(2\alpha r r_i) Y_{lm}^*(\mathbf{\hat{r_i}}) Y_{lm}(\mathbf{\hat{r_i}}),
$$

where $I_l$  is the modified spherical Bessel function of the first kind, which provides the radial dependence.

Bartók also introduced a set of polynomials, $g_n(r)$, which help describe the radial component in a more refined way:

$$
    \phi_\alpha(r) = (rc - r)^{\alpha +2}/N_\alpha
$$

where $N_\alpha$  is a normalization factor given by:

$$
        N_\alpha = \sqrt{\int_0^{rc} r^2(rc-r)^{2(\alpha+2)}dr}
$$


These polynomials are orthonormalized to ensure that the radial functions $g_n(r)$ form a basis. The orthonormalization process is performed through linear combinations of $\phi_\alpha(r)$, and the coefficients are obtained from the overlap matrix $\mathbf{S}$:

$$
    g_n(r) = \sum_{\alpha=1}^{n_{\textrm{max}}}W_{n\alpha}\phi_\alpha(r)
$$

Where $\mathbf{W}$ is constructed from the inverse square root of the overlap matrix $\mathbf{S}$, defined as:

$$
\begin{equation*}
    \begin{split}
       S_{\alpha\beta} &= \int_0^{rc}r^2\phi_\alpha(r)\phi_\beta(r)dr \\ 
       &= \frac{\sqrt{(2\alpha+5)(2\alpha+6)(2\alpha+7)(2\beta+5)(2\beta+6)(2\beta+7)}}{(5+\alpha+\beta)(6+\alpha+\beta)(7+\alpha+\beta)}
    \end{split}
\end{equation*}
$$

The overlap matrix describes how different radial functions overlap with each other and ensures that the final radial basis functions $g_n(r)$ are orthonormal.

The neighbor density function $\rho{\prime}(\mathbf{r})$ can then be expanded in terms of both the radial basis $g_n(r)$ and the spherical harmonics:

$$
            c_{nlm} = \left< g_n(r)Y_{lm}(\mathbf{\hat{r}})|\rho'(\mathbf{r}) \right> 
                    = 4\pi\sum_i^{r_i \leq rc} e^{-\alpha r_i^2} Y_{lm}^*(\mathbf{\hat{r}_i}) 
                     \int_0^{rc} r^2 g_n(r) e^{-\alpha r^2} I_l(2\alpha r r_i) dr
$$


Finally, the **rotation-invariant power spectrum** is obtained by combining these expansion coefficients:

$$
    p_{n_1 n_2 l} = \sum_{m=-l}^{+l}c_{n_1 l m} c^*_{n_2 l m}
$$

This rotation-invariant descriptor provides a comprehensive measure of the local atomic environment by accounting for both radial and angular information, making it a powerful tool for analyzing atomic structures in simulations and experiments.


## 7.3.2 Bispectrum on 4D space

An alternative approach to capturing manybody interactions involves mapping the neighbor density function onto the surface of a 4D hypersphere. This method allows for a richer representation of the local atomic environment by incorporating angular information in 4D.

In this formalism, the coordinates $(x, y, z, r)$ on the 4D hypersphere are given by:

$$
\begin{equation*}
    \begin{split}
        s_1 &= r_0\cos\omega \\
        s_2 &= r_0\sin\omega \cos\theta \\
        s_3 &= r_0\sin\omega \sin\theta \cos\phi \\
        s_4 &= r_0\sin\omega \sin\theta \sin\phi,
    \end{split}
\end{equation*}
$$

$$
\begin{equation*}
    \begin{split}
        r_0 & \geq rc\\
        \theta &= \arccos\left(z/r\right)\\
        \phi &= \arctan\left(y/x\right)\\
        \omega &= \pi r/r_0
    \end{split}
\end{equation*}
$$


Where:

- $r_0$  is a characteristic radius (related to the cutoff radius $rc$).
- $\omega$, $\theta$, and $\phi$ are the spherical coordinates.

The atomic neighbor density function is expressed as:

$$
    \rho(\mathbf{r}) = \delta(\mathbf{r}) + \sum_i{f_c(r)\delta(\mathbf{r}-\mathbf{r_i})}
$$


The first term ensures that the density function remains well-behaved with respect to variations in $\omega$, and $f_c$ is a smooth function to ensure that the it gradually decays to 0 when $r \ge rc$.

By expanding the atomic neighbor density function in terms of the [**Wigner-D matrix**](https://en.wikipedia.org/wiki/Wigner_D-matrix) elements, which represent rotations in the angular coordinates, we obtain:

$$
    \rho(\mathbf{r}) = \sum_{j=0}^{+\infty}\sum_{m',m = -j}^{+j} c^j_{m',m} D^j_{m',m} (2\omega;\theta,\phi)
$$

The coefficients $c_{m{\prime},m}^{j}$ are determined by projecting the neighbor density function onto the Wigner-D matrix elements:

$$
c_{m',m}^{j} = \left\langle D_{m',m}^{j} \middle| \rho \right\rangle
$$

$$
= \int_0^{\pi} d\omega \, \sin^2 \omega \, \omega \int_0^{\pi} d\theta \, \sin \theta \int_0^{2\pi} d\phi \, D_{m',m}^{*j}(2\omega; \theta, \phi) \, \rho(r)
$$

$$
= D_{m',m}^{*j}(0) + \sum_{i} f_{\text{cut}}(r_i) D_{m',m}^{*j}(r_i)
$$
  
Finally, the **bispectrum components**, which capture three-body correlations, can be computed using the triple correlation of the expansion coefficients:

$$
B_{j_1, j_2, j} = \sum_{m', m = -j}^{+j} c_{m',m}^{*j} 
\sum_{m_1', m_1 = -j_1}^{+j_1} c_{m_1',m_1}^{j_1}
\sum_{m_2', m_2 = -j_2}^{+j_2} c_{m_2',m_2}^{j_2} 
C_{m_1 m_2 m}^{j_1 j_2 j} C_{m_1' m_2' m'}^{j_1 j_2 j}
$$

Here, $$C_{m_1 m_2 m}^{j_1 j_2 j}$$  are the **Clebsch-Gordan coefficients**, which ensure proper angular momentum coupling in the bispectrum calculation.

## 7.4 Applications

Manybody descriptors, such as the power spectrum and bispectrum, have found widespread applications in various fields of materials science, condensed matter physics, and machine learning, particularly in the analysis of atomic-scale structures. Their ability to describe both angular and radial components of atomic environments has made them essential tools for understanding complex materials and phenomena.

**Interatomic Potentials and Machine Learning.** One of the most prominent applications of manybody descriptors is in the development of machine-learning-based interatomic potentials. Methods such as the Gaussian Approximation Potential (GAP) and Moment Tensor Potentials (MTP) use manybody descriptors to model the potential energy surface (PES) of atomic systems. These models require descriptors that are invariant to translations, rotations, and permutations of atoms. By providing a compact and invariant representation of the local atomic environment, manybody descriptors allow machine learning models to predict atomic forces and energies with high accuracy, without the need for empirical fitting. This has revolutionized the simulation of large-scale systems, such as materials under extreme conditions or complex chemical reactions.

**Materials Characterization.** Bond order parameters are widely used in molecular dynamics (MD) simulations to distinguish between different crystal structures (e.g., fcc, bcc, hcp) and to identify phase transitions between solid, liquid, and amorphous states. These descriptors allow researchers to quantify the degree of local order or disorder in a material and monitor how this order evolves over time. This is particularly useful in the study of glasses, liquids, and amorphous materials, where traditional descriptors like the pair distribution function fail to capture the full complexity of the atomic arrangement. 

**Phase Transitions.** These descriptors are also widely used to investigate phase transitions in materials. For example, during the melting of a crystalline solid, descriptors like the power spectrum $p_l$ can be employed to monitor changes in the local atomic environment, allowing researchers to pinpoint the onset of disorder as a solid transitions to a liquid. In addition, these descriptors can serve as collective variables in enhanced sampling technqiues to study phase transitions. 

**Physical Property Prediction.** By representing atomic structures in a form that is both compact and invariant, these descriptors enable the construction of high-throughput screening models to predict material properties such as hardness, conductivity, and thermal stability. The use of descriptors like the bispectrum in machine learning pipelines has enabled researchers to explore vast chemical and structural spaces and identify novel materials with desired properties.

**Nanostructures and Catalysis.** Manybody descriptors are also employed to study nanostructures and catalytic surfaces, where the arrangement of atoms plays a key role in determining reactivity and stability. For example, in nanoparticle simulations, these descriptors can be used to characterize the atomic coordination around active sites, providing insight into catalytic behavior. In nanostructured materials, descriptors such as the bispectrum can capture the subtle variations in atomic arrangements that lead to enhanced mechanical or electronic properties.

## 7.5 Conclusion and Further Discussions

The development of manybody descriptors, such as the power spectrum and bispectrum, has provided researchers with powerful tools to describe complex atomic environments in a rigorous and invariant manner. These descriptors have addressed the limitations of simpler pairwise and angular metrics, enabling more accurate characterization of local atomic arrangements. Whether in the context of interatomic potentials, structural analysis, phase transitions, or materials discovery, manybody descriptors have significantly advanced our understanding of atomic-scale phenomena.


