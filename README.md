# Meshfree quadrature rule
This repository contains the code for the paper 

* [An asymptotically compatible meshfree quadrature rule for nonlocal problems with applications to peridynamics](https://www.sciencedirect.com/science/article/pii/S004578251830402X)

There are three test problems: static bond-based peridynamics, dynamic nonlocal diffusion, and Kalthoff-Winkler experiment. All the problems are in two-dimensional space.

# Requirements
* gcc 7.5.0 or newer
* Install BLAS and LAPACK packages, you may need to change Makefile to direct to BLAS and LAPACK libraries.

# Files
* PMB_2Dweight.cpp is the code for the static bond-based peridynamics with Dirichlet boundary condition.
* nonlocaldiff.cpp is the code for the dynamic nonlocal diffusion problem with Dirichlet boundary condition and backward Euler scheme.
* KW_2Dweight_dynamic.cpp is the code for the Kalthoff-Winkler experiment simulation.

# Quick start
## static bond-based peridynamics

compile the run code with the following command:

1. make PD
2. ./PMB_2d.ex \<number of particles\>\<ratio of delta/h\>

Example:

./PMB_2d.ex 20 3.5 0.0

## Nonlocal diffusion problems

compile the run code with the following command:

1. make Nldiff
2. ./nldiff.ex 

## Kalthoff-Winkler experiment

compile the run code with the following command:

1. make KW
2. ./KW.ex \<number of particles\>\<ratio of delta/h\>

Example:

./PMB_2d.ex 32 3.0


# Additional information 

If you find this code or method useful for your project, please cite

@article{trask2019asymptotically,

  title={An asymptotically compatible meshfree quadrature rule for nonlocal problems with applications to peridynamics},
  author={Trask, Nathaniel and You, Huaiqian and Yu, Yue and Parks, Michael L},
  journal={Computer Methods in Applied Mechanics and Engineering},
  volume={343},
  pages={151--165},
  year={2019},
  publisher={Elsevier}
}

