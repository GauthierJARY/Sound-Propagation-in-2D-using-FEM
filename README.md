### Introduction
This repository contains MATLAB code implementing the Finite Element Method (FEM) for solving partial differential equations (PDEs). The code is organized into different sections corresponding to specific problem formulations and numerical techniques. This README provides an overview of the contents and usage of the code.
For more information, please refer yourself to the written technical report.
### Table of Contents
1. [Stationary Problem with Variable Coefficients and Dirichlet Conditions](#stationary-problem)
2. [Wave Equation with Exact Mass Matrix Integration](#wave-equation)
3. [Wave Equation using Mass Condensation Method](#wave-equation-condensation)

### Stationary Problem with Variable Coefficients and Dirichlet Conditions
- **Formulation and Well-Posedness**: The problem is formulated variably with Dirichlet conditions.
- **Discretization**: The variational formulation is discretized.
- **Matrix Representation of the Linear Problem**: Representation of the linear problem in matrix form.
- **Pseudo-Elimination Technique**: Application of pseudo-elimination technique for solving linear systems.
- **Simulation Domain**: Description of the simulation domain.
- **Mapping to Reference Triangle**: Transformation to the reference triangle.
- **Calculation of Basis Functions on the Reference Triangle**: Method for calculating basis functions.
- **Elementary Mass Matrix**: Exact expression for the elementary mass matrix.
- **Elementary Stiffness Matrix**: Approximation using Gauss-Legendre quadrature.
- **Finite Element Method Framework**: Implementation of the finite element method.
- **Algorithm Testing and Results Analysis**: Testing the algorithm and analyzing the results.

### Wave Equation with Exact Mass Matrix Integration
- **Variational Formulation and Physical Link to Energy**: Formulation of the wave equation with exact mass matrix integration and its relation to energy.
- **Discretization of the Variational Formulation**: Discretization of the variational formulation of the wave equation.
- **Temporal Discretization of the Problem**: Techniques for discretizing the temporal aspect of the problem.
- **Numerical Implementation**: Details of the numerical implementation.
- **CFL Condition**: Explanation of the Courant-Friedrichs-Lewy (CFL) condition.
- **Propagation Function**: Function describing the propagation of waves.
- **Results Analysis**: Analysis of the simulation results.

### Wave Equation using Mass Condensation Method
- **Mass Matrix Condensation**: Application of mass matrix condensation technique.
- **Second Application Case of the Method**: Description of a second case where the mass condensation method is applied.
- **Homogeneous Case**: Considerations for the homogeneous case.

For detailed instructions on running the code and understanding each section, please refer to the corresponding MATLAB files and comments within the code. If you have any questions or issues, feel free to contact me.

