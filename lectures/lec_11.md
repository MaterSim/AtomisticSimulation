# Week 11. Introduction to VASP

For practical application, it is not recommended to run DFT calculation with your own code. In the past, researchers have developped many excellent choice of DFT-based codes. Among them,  the Vienna Ab initio Simulation Package (VASP) is one of the most widely used computational tools for performing first-principles calculations based on Density Functional Theory (DFT), Hartree-Fock (HF), and hybrid functionals. 
It employs plane-wave basis sets and pseudopotentials or PAW (Projector Augmented-Wave) potentials to represent the wavefunctions and the core-electron interaction, respectively.
VASP is highly efficient in calculating the electronic structure, total energies, and properties of materials.

## 11.1 Key Concepts in VASP

1. **Plane-Wave Basis Sets**. In VASP, the electronic wavefunctions are expanded in terms of plane waves. The accuracy of this expansion is controlled by the energy cutoff  E_{\text{cutoff}} , which limits the maximum kinetic energy of the plane waves included in the calculation.

2. **Pseudopotentials and PAW**. To speed up calculation, Pseudopotentials are used to replace the core electrons with a smoother potential, allowing only the valence electrons to be explicitly treated. In VASP, the so called PAW Potentials were used to reconstruct the all-electron wavefunction from the pseudopotential wavefunction.

3. **Exchange-Correlation Functionals**. VASP support a variety of exchange-correlation functionals to account for the interactions between electrons, including LDA, GGA or more advanced hybrid functionals.

4. **K-Point Sampling**. The Brillouin zone of a crystal is sampled using k-points. VASP uses a Monkhorst-Pack grid or Gamma-centered grid to sample the k-points. Denser k-point grids yield more accurate results, especially for band structure calculations.

5. **Self-Consistent Field (SCF) Iteration**. In VASP, the SCF process iteratively solves the Kohn-Sham equations to obtain the electronic ground state. This iterative procedure continues until convergence is achieved based on the total energy or electronic density.


## 11.2 Workflow for a VASP Calculation

A typical VASP calculation consists of several input files and a few critical steps:

### 11.2.1 Input Files

- POSCAR: Defines the crystal structure (atomic positions, lattice vectors, and unit cell).
- INCAR: Contains the calculation parameters (energy cutoff, type of calculation, convergence criteria).
- POTCAR: Contains the pseudopotentials for the elements in the system.
- KPOINTS: Defines the k-point mesh used to sample the Brillouin zone.

### 11.2.2 Output Files

- OUTCAR: Provides detailed information about the entire calculation, including the total energy and convergence status.
- CONTCAR: Contains the relaxed atomic positions if ionic relaxation was performed.
- EIGENVAL: Contains the eigenvalues, which are used for band structure and DOS calculations.
- DOSCAR: Contains the density of states data.
Many other files 

## 11.3 Example: Band Structure Calculation for Silicon

Let’s walk through an example to compute the band structure of silicon (Si) using VASP.

## 11.3.1 Preparing the Input Files

POSCAR: Silicon Crystal Structure
```
Si
1.0
   0.000000  2.715000  2.715000
   2.715000  0.000000  2.715000
   2.715000  2.715000  0.000000
Si
2
Direct
  0.000000  0.000000  0.000000
  0.250000  0.250000  0.250000
```

INCAR: General Input Parameters
```
SYSTEM = Silicon
PREC = Accurate
ENCUT = 400
ISMEAR = 0
SIGMA = 0.05
IBRION = -1
ISIF = 2
NSW = 0
LWAVE = .FALSE.
LCHARG = .FALSE.
EDIFF = 1E-6
```

POTCAR: Pseudopotential File (from each VASP distribution)
```
Si_PBE
```
KPOINTS: K-Point Mesh for SCF Calculation
```
Automatic mesh
0
Monkhorst-Pack
8 8 8
0 0 0
```


Run the SCF Calculation

Run VASP using the prepared input files to get the ground-state electronic density. This is required before computing the band structure.

After the SCF run, you’ll get files like OUTCAR, EIGENVAL, and CHGCAR.


Plotting the Band Structure

To complete

