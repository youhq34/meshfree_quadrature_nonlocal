# Meshfree quadrature rule
This repository contains the code for the paper 

* [An asymptotically compatible meshfree quadrature rule for nonlocal problems with applications to peridynamics](https://www.sciencedirect.com/science/article/pii/S004578251830402X)

There are three test problems: static bond-based peridynamics, dynamic nonlocal diffusion, and Kalthoff-Winkler experiment simulation. All the problems are in two-dimensional space.

# Requirements
* gcc 7.5.0 or newer
* Install BLAS and LAPACK packages

# Quick start
## static bond-based peridynamics

compile the run code with the following command 

1. make PD
2. ./PMB_2d.ex <number of particles> <ratio of delta/h> <random perturbation coeffcient>





