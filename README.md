# README 

`Fermi` is a code made for studying the physical behavour of nuclear reactors and to perform a design of them.

It solves an eigenvalue problem corresponding to the steady state multigroup diffusion equation, which is an aproximation of the Boltzman Transport Equation. This aproximation is only valid in materials were the relation between absortion and scattering collision of neutrons is small.

Macroscopic cross sections should be provided in order to solve the equations in the domain. They may become from a lattice cell homogenization, and possibly, a condensation to few groups in order to decrease the computational effort.

To perform calculations and solve the problem, `Fermi` take advantage of the finite element method to discretizise the equations. With this technique, the code is capable of solving the problem over unstructured 3d meshes.

The design of `Fermi` was aim to be a very simple and easy understanding code. The input is a 1d, 2d or 3d mesh that should be generated with `gmsh` code and an ASCII file with a particular format which contains information about the cross sections of each material that exist in the domain. The output is a `VTK` file which contains the solution of the problem (the scalar flux) and the a file containing the eigenvalue of multiplication factor "keff" of the problem.

## Instalation

For reading this text in a `pdf` format do:

```bash
pandoc README.md -V geometry:margin=.5in --latex-engine=xelatex -o README.pdf
```

###`PETSC` library

Download it from [www.mcs.anl.gov/petsc](www.mcs.anl.gov/petsc) and do:

###`SLEPc` library:

Download it from [http://slepc.upv.es](http://slepc.upv.es) and do:

###`Compilation`:

```bash
   make
```

###`Input Structure`:

```bash
$Mesh 
    mesh_file   fuel_1.msh
$EndMesh 

$Mode
  timedep QSTATIC 
  p0 2.5e6  
  t0 0.0
  tf 0.15
  dt 0.05
$EndMode

$Xs
  egn 1
  #          F D   XA  nXF  eXF CHI
  "FUEL_Z1"  0 1.5 0.2   0.5  0.4 1.0
  "FLUID_Z1" 0 1.5 0.005 0.0  0.0 1.0
$EndXs

$Boundary
  "TOP_SURF"    1 0
  "BOT_SURF"    2 0
  "LAT_SURF"    3 1
$EndBoundary

$Output
  # Power in physical entities
  kind 2
  file pow_phys.dat
  nphy 2
  "FUEL_Z1" "FLUID_Z1"
$EndOutput

$Communication
  kind   1
  friend control
  nphy   2
  "FUEL_Z1" "FLUID_Z1"
$EndCommunication
```

## The future  

Paralelization and performance evaluation 

* Benchmarking

* Documentation

Guido Giuntoli - [giuntoli1991@gmail.com]
