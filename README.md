# README 

`Fermi` is a code made for studying the physical behavour of nuclear reactors and to perform a design of them.

>It solves an eigenvalue problem corresponding to the steady state multigroup diffusion equation, which is an aproximation of the Boltzman Transport Equation. This aproximation is only valid in materials were the relation between absortion and scattering collision of neutrons is small.
>
>Macroscopic cross sections should be provided in order to solve the equations in the domain. They may become from a lattice cell homogenization, and possibly, a condensation to few groups in order to decrease the computational effort.
>
>To perform calculations and solve the problem, `Fermi` take advantage of the finite element method to discretizise the equations. With this technique, the code is capable of solving the problem over unstructured 3d meshes.
>
>The design of `Fermi` was aim to be a very simple and easy understanding code. The input is a 1d, 2d or 3d mesh that should be generated with `gmsh` code and an ASCII file with a particular format which contains information about the cross sections of each material that exist in the domain. The output is a `VTK` file which contains the solution of the problem (the scalar flux) and the a file containing the eigenvalue of multiplication factor "keff" of the problem.

## Instalation

For reading this text in a `pdf` format do:

```bash
pandoc README.md -V geometry:margin=.5in --latex-engine=xelatex -o README.pdf
```

<!--###`OpenMPI` installation: 

Download package from [www.open-mpi.org](www.open-mpi.org) and do:

```bash
gunzip -c openmpi-1.10.3.tar.gz | tar xf -
cd openmpi-${version}
```

and with *root* privileges

```bash
./configure --prefix=/usr/local
make all install
```

to check the instalation do:

```bash
   locate mpicc mpirun 
   mpicc --version
   mpirun --version  
```-->

###`BLAS & LAPACK` libraries

```bash
apt-get install libblas-dev liblapack-dev
```  

###`PETSC` library

Download it from [www.mcs.anl.gov/petsc](www.mcs.anl.gov/petsc) and do:

```bash    
   export PETSC_DIR=/path-to-petsc-installation-directory
   export PETSC_ARCH=arch-linux2-c-opt    
   LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-petsc-directory
```  
we recommend to put them in *.bashrc* or in some start up file

```bash
   tar -xvzf petsc-{version}
   ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-mpich
   make all 
   make test
```

###`SLEPc` library:

Download it from [http://slepc.upv.es](http://slepc.upv.es) and do:

```bash   
   export SLEPC_DIR=/path-to-slepc-installation-directory 
```

we recommend to put them in `.bashrc` or in some start up file

```bash
   tar -xvzf slepc-{version} 
   ./configure
   make all
   make test
```
###`Fermi`:

```bash
   make
```

##Running `Fermi` on multiple processors

```bash
   mpirun -np 2  fermi -mesh <mesh_file> -mat <mat_file> -vtk <vtk_file>
```

##Debbuging `Fermi` on multiple processors

```bash
   mpirun -np 2  xterm -e gdb --args fermi \
                 -mesh  <mesh_file>        \ 
                 -mat   <mat_file>         \
                 -vtk   <vtk_file>         \
                 -epmv  <epmv_file>        \
                 -npmv  <npmv_file>        \
                 -epms  <epms_file>        \
                 -npms  <npms_file> 
```   

## The future  

Paralelization and performance evaluation 

* Scripts to test the code

* Benchmarking

* Documentation

Guido Giuntoli - [giuntoli1991@gmail.com]
