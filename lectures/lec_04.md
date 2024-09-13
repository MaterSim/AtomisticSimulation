# Week 4: Realistic MD simulation with LAMMPS

## 4.1 Introduction to LAMMPS
LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is an open-source software tool designed for performing classical molecular dynamics (MD) simulations. It can be used to model an array of particle interactions, ranging from simple atomic systems to complex materials and biomolecular systems. As one of the most population materials simulation package, LAMMPS is specifically optimized for large-scale simulations, which involve millions to billions of particles, making it suitable for high-performance computing (HPC) environments. Its versatility allows for simulations of a variety of systems such as metals, polymers, proteins and membranes. Some typical applications include

- Crystal Defects and Deformation: (e.g., dislocation motion, grain boundary evolution, and fracture in materials).
- Phase Transitions: Simulating phase changes in metals and alloys, such as melting, solidification, or the formation of microstructures.
- Transport properties of materials via Green-Kubo or direct methods.
- Mechanical Properties such as stress-strain relationships, elasticity, and plasticity at the atomic scale.
- Drug Interactions: Simulating how drug molecules interact with proteins or other biological targets.
- Investigating lipid bilayers, ion channels, or membrane proteins.

## 4.2 Why is LAMMPS efficient and popular?
LAMMPS is designed to perform large-scale simulations by taking advantage of parallel computing architectures, including multi-core processors and high-performance computing (HPC) clusters. In particular, LAMMPS uses the domain decomposition technique to divide the simulation space into subdomains. Each processor or core is assigned a subdomain, and they work together by communicating boundary conditions and interacting forces. LAMMPS also uses MPI to handle communication between processors, ensuring minimal overhead and efficient data transfer during the simulation. Thanks to these considerations, LAMMPS has demonstrated excellent scalability across thousands of processors, which makes it suitable for simulating systems with millions to billions of particles over long time scales.

- Support for Multiple Interatomic Potentials (LJ, Embedded Atom Method, REAXX.etc)
- Flexible Input Script: LAMMPS uses a flexible scripting language for setting up simulations. It allowing users to define new potential types, set conditions, and run specific simulation protocols.

## 4.3 Input and Output in LAMMPS

After you compiled the LAMMPS code into an executable (often called `lmp_serial` or `lmp_mpi` by default). 

```
path_to_lmp_mpi < lmp.in > lmp.out
```

This command involves the preparation of input and output that would be discussed as follows.

### 4.3.1 LAMMPS Input
LAMMPS simulations are controlled by input scripts, which consist of a series of commands written in a simple text format. These scripts define the simulation parameters, system setup, and specific instructions for running the simulation.

A typical LAMMPS input script is organized into several key sections:

**Initialization Section** includes global settings, such as the units, boundary conditions, and atom styles).
- ``units real``: Defines the units for physical quantities (e.g., distance in angstroms, energy in kcal/mol).
- ``boundary p p p``: Sets periodic boundary conditions in all three spatial dimensions.
- ``atom_style atomic``: Defines how atoms are represented (e.g., atomic, charge, molecular).

**Atom Definition Section** includes Atomic coordinates, and initial velocities.
- ``read_data data.file``: Reads the atomic configuration and other system properties from a file .
- ``velocity all create 300.0 12345``: Assigns random initial velocities to atoms at a temperature of 300 K with a random seed.

**Force Field Definition Section**
- ``pair_style lj/cut 2.5``: Defines a Lennard-Jones potential with a cutoff distance of 2.5.
- ``pair_coeff * * 0.1 3.0``: Sets the Lennard-Jones coefficients (epsilon and sigma) for the interactions between atom types.

**Simulation Parameters Section**
- ``timestep 1.0``: Sets the time step for integration to 1.0 (in the units defined by the units command).
- ``fix 1 all nve``: Applies a constant energy (NVE) ensemble to all atoms.
- ``fix 2 all temp/rescale 100 300.0 300.0 0.02 1.0``: Rescales the temperature every 100 timesteps to maintain a temperature of 300 K.

**Output Control Section** specifies the frequency and format of the simulation output.
- ``thermo 100``: Prints thermodynamic data (e.g., temperature, pressure, energy) every 100 timesteps.
- ``dump 1 all atom 1000 dump.atom``: Outputs the atomic positions to a file (dump.atom) every 1000 timesteps.

**Run Section**
- ``run 10000``: Runs the simulation for 10,000 timesteps.
- ``minimize 1.0e-4 1.0e-6 1000 10000``: Performs energy minimization on the system with specified tolerance and iteration limits.

### 4.3.2 LAMMPS Output
A typical calculation in LAMMPS may generate the following output files.

- **The log file** is generated automatically for every LAMMPS run and records all the commands executed in the input script. It also contains thermodynamic data (e.g., energy, pressure, temperature) at intervals specified by the thermo command.
```
Step Temp E_pair E_mol TotEng Press Volume
0 300.0 -143.53 0 -123.46 1.35 1000.0
100 310.2 -142.11 0 -122.45 1.40 1000.0
```

- **Dump Files** store detailed trajectory information about the system’s atomic coordinates, velocities, and forces. These files are typically used for post-processing to analyze system configurations, create visualizations, or calculate structural properties.

```
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
1000
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
3 2 0.5 0.5 0.5
```

- **Restart Files** store the entire state of a LAMMPS simulation, allowing users to pause and later continue a simulation from where it left off. These files contain information about the atom positions, velocities, forces, and other system properties.


## 4.4 Simulation Process

LAMMPS begins with setting up the system and defining parameters. This involves specifying the atomic system geometry, simulation box, interatomic potentials, and initial conditions like atomic velocities and temperature.

Below is a simulation of Argon atoms using the Lennard-Jones potential as we discussed in the previous lectures:
```
units real
atom_style atomic
read_data argon.data
pair_style lj/cut 2.5
pair_coeff * * 0.238 3.4   # Argon-specific parameters for Lennard-Jones potential
velocity all create 300.0 12345

fix 1 all nvt temp 300.0 300.0 100.0   # NVT ensemble to control temperature
timestep 1.0                           # Time step of 1.0 fs
run 50000                              # Run simulation for 50,000 timesteps
```

## 4.5 Post-Processing

After a simulation, the results need to be visualized and analyzed. LAMMPS produces several types of output files, which contain thermodynamic data, atom positions, velocities, and forces. VMD (Visual Molecular Dynamics) and OVITO (Open Visualization Tool) are popular tools for visualizing molecular dynamics simulations.

To add a few ovito figures and explain ovito.

## 4.6 Running LAMMPS on HPC
For most research projects, running LAMMPS on a high-performance computing (HPC) environment is essential for large-scale simulations that require significant computational resources. Most modern supercomputers use job schedulers like SLURM to manage computational tasks.

```
#!/bin/bash
#SBATCH --job-name=lammps_job               # Job name
#SBATCH --nodes=4                           # Number of nodes
#SBATCH --ntasks-per-node=32                # Number of tasks (processes) per node
#SBATCH --time=24:00:00                     # Max time limit (HH:MM:SS)
#SBATCH --partition=compute                 # Partition or queue to submit to
#SBATCH --output=job_output.log             # Output log file

module load lammps/3Mar2020                 # Load LAMMPS module
mpirun -np 128 lmp_mpi -in input_file.in    # Run LAMMPS in parallel across 128 processes
```


- `SBATCH --options` are used to specify the number of nodes, tasks, job name, and time limit.
- `mpirun -np 128` launches LAMMPS across 128 processes in parallel, ensuring that the simulation scales across multiple cores.
- `lmp_mpi` is the parallel version of LAMMPS used for multi-node execution.


