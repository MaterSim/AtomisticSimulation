# MEGR 7090/8090: Atomistic Simulation in Materials Modeling
 
## Course Introduction
This is a 3-credit course requires three hours of classroom or direct faculty instruction and six hours of out-of-class student work for the equivalent of approximately 15 weeks. Out-of-class work may include but is not limited to: required reading; homework; studying for quizzes and exams; research; written assignments; and project design, simulation, testing and demonstration.

## Instructor: Qiang Zhu (Battcave 114, qzhu8@uncc.edu)

## Text book, title, author, and year: 
- *Understanding molecular simulation from algorithms to applications*, By Daan Frankel and Berend Smit, 3rd Edition
- *Electronic structure*, By Richard. Martin, 2nd Edition

The lecture notes were made based on these two excellent books. However, the hard copies of textbooks are not strictly required. We will also keep updating this lecture notes and provide more open access video or text tutorials throughout the course.

## Course Description
This course aims to use the atomistic computer simulation to model and understand the properties of real materials and their accompanying process and phenomena. It will primarily focus on two approaches: molecular dynamics and electronic structure calculation based on density functional theory. Some typical examples, codes, analytical tools will be also covered in this course. 

The expected outcomes include: 
- Understand the fundmental of Molecular dynamics simulation and its connection with statistical physics
- Apply the molecular dynamics simulation technique to model the physical process in real materials
- Understand the concept of electronic structure simulation based on density functional theory 
- Use the available softwares LAMMPS and VASP to compute material’s properties

## Tenative schedules

### I: Molecular dynamics simulation
- Week 1: Motivating example 1: Numerical simulation of gas under the NVE ensemble
- Week 2: Motivating example 2: Liquid-gas phase transition under the NVT ensemble
- Week 3: Motivating example 3: Simulation of solids under the NPT ensemble 
- Week 4: Introduction to the LAMMPS package
- Week 5: MD Analysis I: structural characterization (RDF), degree of order
- Week 6: MD Analysis II: transport processes
- Week 7: Selected presentations from students

### II: Electronic structure calculation
- Week 8: Gentle introduction to Density functional theory 
- Week 9: Motivating example 4: DFT treatment of H2 molecule
- Week 10: From molecule to the periodic system
- Week 11: Introduction to VASP
- Week 12: Band structure analysis
- Week 13: Phonon calculation from both classical force field and DFT
- Week 14: Selected presentations from students 


# Week 1: Gas Simulation under the NVE ensemble

## Prehistory of Computer Simulation:
* The Los Alamos MANIAC (Mathematical Analyzer, Numerical Integrator, and Computer) became operational in 1952. This event marks a significant milestone in the history of computing. Nicholas Metropolis, as the most notable early user, developed the Monte Carlo method, a statistical technique that utilizes random sampling to solve complex mathematical problems. 

* The launch of computer also openned the door for the study of many fundamental problems. Most of the material systems consist of many atoms or molecules, how can we infer the properties of such systems? In the past, people have to to do it either analytically (e.g., thermodynamics and stastical mechanics have been developed to study some classical systems such as ideal gas, Ising Model, ferromagentic phase transition and alloys. Some analytic solutions can be derived). They are very intelligent but lacks the atomic detail. An alterative approach is directly model the system (straightforward but very time consuming and tedious). Notable examples include. 
1. [Buffon's needle experiment to compute $\pi$](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem), 
2. [Bernal's ball-bearing model](https://iopscience.iop.org/article/10.1088/0953-8984/26/46/463102), 
3. [Kitaigorodskii’s structure seeker](https://pubs.acs.org/doi/10.1021/acs.cgd.8b00972).

* Molecular dynamics Simulation is generally a technical to directly study the atomic evolution of atoms or molecules in the material system based on the simple law of Newtonian Dynamics. Immediately after the invention of computers, MD simulations have been quickly applied to several systems
1. First ever, 1956, Alder and Wainwright (Livermore), [Phase transition of hard spheres](https://gibbs.ccny.cuny.edu/teaching/s2021/labs/HardDiskSimulation/Alders&Wainwright1957.pdf)
2. Real material, 1959, Gibson, Goland, Milgram, Vineyard (Brookhaven), [Dynamics of radiation damage](https://journals.aps.org/pr/abstract/10.1103/PhysRev.120.1229)
3. Liquid, 1964, Rahman, [Correlations in Liquid Argon](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405)

## Why do we need such kind of direct simulation method like MD
* Go beyond the experimental capability
* Gain some physical insights

Homework: Read the [Rahman](https://en.wikipedia.org/wiki/Aneesur_Rahman)'s 1964 paper. We will try to reproduce some results from this work in this week. 

## A first MD simulation of liquid argon under NVE ensemble.

### MD workflow

### Interatomic potential
Lennard Jones Potential is a popular choice for modelling weakly interacted systems.

* Energy
* Force


### Initialization

### Integrator (updating rule)
* Leapfrog
* Verlet

### Full code

Hopefully, you are able to write a basic code for MD simulation. 

Of course, there are many excellent open-source MD codes with more functional support. For productive research project, you would probably use those codes. In this course, we recommend the use of LAMMPS, one of the most popular code for materials modelling. A short example is provided as follows.
   







