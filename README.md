# 1D Finite Element Method

### Author: Carel Soney || eid: cs63674

This Github Repository contains: 

1) README.md
2) project_2_galerkin.m
3) FEM and BEM Plots.pdf
4) weak form.pdf 

1 - README.md 

This markdown file contains information regarding the other files located in this repository. 

2 - project_2_galerkin.m 

This Matlab file contains the code that computes the finite element method and applies it to the heat transfer problem. The backward and forward euler method code is also located in this file. 

3 - FEM and BEM Plots. pdf

This pdf file contains the different graphs for the problem computed using both the Forward Euler Methods (FEM) and Backward Euler Method (BEM). 

4 - weak form.pdf

This pdf file answers the first part of this project which would be deriving the weak form of the given equation. 

## Part A: Galerkin Weak Form 

The weak form transforms a differential equation by relaxing point-wise constraints through the multiplication of a test function, v(x), and integrating over a given domain. The Galerkin method allows the problem to be discretized and represents the test functions are a linear combination defined on a mesh. This leads to a system of algebraic equations that can then be solved numerically. 

## Part B: Forward Euler Method

The stability of the FEM is closely tied to the chosen time step, Δt, or dt as in the code. This numerical method is conditionally stable, meaning there are values of dt in which the solution becomes unstable, leading to inaccurate results. To find the instability, increase the time step and plot the graph to see how the function behaves. 

## Part C: Backward Euler Method 

Unlike the FEM, the BEM is unconditionally stable, meaning it is a reliable numerical method regardless of the value of dt. When dt is equal to the spatial step size, it helps the implicit BEM nullify stability constraints that can be seen in the FEM. If dt is greater than the spatial step size, BEM will still remain stable due to its implicit nature. 
