# 7 Representation of Local Atomic Environment

## 7.1 Pair Distribution Function

## 7.2 Orientational Order Parameter

### 7.2.1 Bond Orientational Order Parameter:
Let's first consider a two-dimensional systems, we can quantify the angular relationships between nearest neighbors as follows,

$$
\psi_m = \frac{1}{N} \sum_{j=1}^{N} e^{i m \theta_j}
$$

Where:
- $N$ is the number of neighbors around a given particle.
- $\theta_j$ is the angle formed between the bond connecting a particle to its neighbor $j$ and some reference axis.
- $m$ is the symmetry index. For example, $m$ = 6 is used for systems with hexagonal symmetry.

$\psi_m$  will take values near 1 for highly ordered systems (where neighbors align in a regular, symmetrical pattern), and it will be closer to 0 for disordered systems.

### 7.2.2 Spherical Harmonics-Based Order Parameter:

The **atomic neighbor density function** is the spatial distribution of atoms in a structure within a cutoff radius ($r_{\textrm{cut}}$),

$$
    \rho(\bm{r}) = \sum_i^{r_i \leq r_{\textrm{cut}}} \delta(\bm{r}-\bm{r_i})
$$

Expanding $\rho(\bm{r})$ as a series on the 2-sphere using spherical harmonics.

$$
    \rho(\bm{r}) = \sum_{l=0}^{+\infty}\sum_{m=-l}^{+l}c_{lm}Y_{lm}(\bm{\hat{r}}),
$$
where the expansion coefficients $c_{lm}$ are given by:

$$
    c_{lm} = \left<Y_{lm}(\bm{\hat{r}})|\rho(\bm{r})\right> = \sum_i^{r_i \leq r_{\textrm{cut}}}Y_{lm}(\bm{\hat{r}_i}).
$$


Hence, the expansion coefficients can reflect the neighbor density information. However, these coefficients would change if one rotates the system, which is not 

Steinhardt constructed his [**bond order parameters**](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784) using second and third order combinations of the expansion coefficients to quantify order in liquids and glasses.

$$   
p_l = \sum_{m=-l}^{+l}c_{lm}c^*_{lm}
$$

Note that the derived $p_l$ actually 

The limitations of $p_{l}$
-  no any radial information.
-  only describes a neighbor density in the form of $\delta$ function.

Despite these limitation, the ideal of using spherical harmnics and power spectrum has inspired later works.

## 7.3 Manybody descriptors
In real application, we often care about how atoms are spatially arranged in both radial and angular space. In the *** paper, 

### 7.3.1 Power Spectrum with the Explicit Radial Component
To address the limitation of descrbing neighbor density based on the $\delta$ function, we can turn it to a Gaussian with a limited width ($\alpha$). As such, 

$$
    \rho'(\bm{r}) = \sum_i^{r_i \leq r_{\textrm{cut}}} \exp(-\alpha|\bm{r}-\bm{r_i}|^2),
$$

Expanding it on the 2-sphere through a spherical harmonic transform of $e^{2\alpha\bm{r}\cdot \bm{r_i}}$.

$$
\begin{equation*}
\begin{split}
     \rho'(\bm{r}) &= \sum_{r_i \leq r_{\textrm{cut}}} e^{-\alpha(r^2+r_i^2)} e^{2\alpha\bm{r}\cdot \bm{r_i}}\\
     & =\sum_{r_i \leq r_{\textrm{cut}}} \sum_{lm}  4\pi e^{-\alpha(r^2+r_i^2)} I_l(2\alpha r r_i) Y^*_{lm}(\bm{\hat{r}_i})Y_{lm}(\bm{\hat{r}}),
\end{split}
\end{equation*}
$$

Radial information can be explicitly added.  
Bartok introduced a set of polynomials, $g_n(r)$, 

$$
    \phi_\alpha(r) = (r_{\textrm{cut}} - r)^{\alpha +2}/N_\alpha
$$

where
$$

        N_\alpha = \sqrt{\int_0^{r_{\textrm{cut}}} r^2(r_{\textrm{cut}}-r)^{2(\alpha+2)}dr}
$$

$I_l$ is a modified spherical Bessel function of the first kind. 

Orthonormalizing linear combinations of $\phi_\alpha$.

$$
    g_n(r) = \sum_{\alpha=1}^{n_{\textrm{max}}}W_{n\alpha}\phi_\alpha(r)
$$

$\bm{W}$ is constructed from the overlap matrix $\bm{S}$ by the relation $\bm{W}=\bm{S^{-1/2}}$.  

The overlap matrix is given by the inner product 

$$
\begin{equation*}
    \begin{split}
       S_{\alpha\beta} &= \int_0^{r_{\textrm{cut}}}r^2\phi_\alpha(r)\phi_\beta(r)dr \\ 
       &= \frac{\sqrt{(2\alpha+5)(2\alpha+6)(2\alpha+7)(2\beta+5)(2\beta+6)(2\beta+7)}}{(5+\alpha+\beta)(6+\alpha+\beta)(7+\alpha+\beta)}
    \end{split}
\end{equation*}
$$

Expanding $\rho'(\bm{r})$ on the 2-sphere and radial basis $g(r)$,

$$
\begin{equation*}
    \begin{split}
            c_{nlm} &= \left<g_n(r)Y_{lm}(\bm{\hat{r}})|\rho'(\bm{r})\right> \\ 
            &= 4\pi\sum_i^{r_i \leq r_{\textrm{cut}}} e^{-\alpha r_i^2}Y^*_{lm}(\bm{\hat{r}_i})\int_0^{r_{\textrm{cut}}}r^2g_n(r)e^{-\alpha r^2}I_l(2\alpha r r_i)dr
    \end{split}
\end{equation*}
$$

Finally, the rotation invairant power spectrum can be obtained

$$
    p_{n_1 n_2 l} = \sum_{m=-l}^{+l}c_{n_1 l m} c^*_{n_2 l m}
$$

## 7.3.2 Bispectrum on 4D space

Alternatively, one can map the neighbor density function within a cutoff radius onto the  surface  of 4D hyper-sphere,

$$
\begin{equation*}
    \begin{split}
        s_1 &= r_0\cos\omega \\
        s_2 &= r_0\sin\omega\cos\theta \\
        s_3 &= r_0\sin\omega\sin\theta\cos\phi \\
        s_4 &= r_0\sin\omega\sin\theta\sin\phi,
    \end{split}
\end{equation*}
$$

$$
\begin{equation*}
    \begin{split}
        r_0 & \geq r_{\textrm{cut}}\\
        \theta &= \arccos\left(z/r\right)\\
        \phi &= \arctan\left(y/x\right)\\
        \omega &= \pi r/r_0
    \end{split}
\end{equation*}
$$


The atomic neighbor density function is expressed as  

$$
    \rho(\bm{r}) = \delta(\bm{r}) + \sum_i{f_{\textrm{cut}}(r)\delta(\bm{r}-\bm{r_i})}
$$


The first term is to avoid variance with respect to $\omega$.


Now, the atomic neighbor density function can be represented in an expansion of Wigner-$D$ matrix elements in the angle-axis representation, where $2\omega$ is the rotation angle and $\theta,\phi$ define the axis.

$$
    \rho(\bm{r}) = \sum_{j=0}^{+\infty}\sum_{m',m = -j}^{+j}{c^j_{m',m}D^j_{m',m}\left(2\omega;\theta,\phi\right)}
$$

$$    
\begin{equation*}
    \begin{split}
    c^j_{m',m} &= \left<D^j_{m',m}|\rho\right> \\&= \int_0^\pi{d\omega \sin^2\omega}\int_0^\pi{d\theta\sin\theta}\int_0^{2\pi}d\phi D^{*j}_{m',m}\left(2\omega;\theta,\phi\right)\rho(\bm{r}) \\&=
    D^{*j}_{m',m}(\bm{0}) + \sum_i f_{\textrm{cut}}(r_i)D^{*j}_{m',m}(\bm{r_i})
    \end{split}
\end{equation*}    
$$

  
To obtain the **bispectrum components**, the **triple-correlation** of $c^j_{m',m}$ is used

$$
    B_{j_1,j_2,j} = \sum_{m',m = -j}^{+j}c^{*j}_{m',m}\sum_{m_1',m_1 = -j_1}^{+j_1}c^{j_1}_{m_1',m_1}\times  \sum_{m_2',m_2 = -j_2}^{+j_2}c^{j_2}_{m_2',m_2}C^{jj_1j_2}_{mm_1m_2}C^{jj_1j_2}_{m'm_1'm_2'}
$$

where $C$ is a Clebsch-Gordan coefficient.

Note that the **power spectrum** can be considered the **auto-correlation** of the expansion coefficients in the language of signal process.


## 7.4 Applications

