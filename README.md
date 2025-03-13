
# Heat Equation Simulation

This repository contains Python code for simulating the heat equation using various numerical methods.

## Methods Implemented

*   **Explicit Finite Difference Method**: This method uses a forward difference in time and a central difference in space to approximate the solution of the heat equation.
*   **Implicit Finite Difference Method**: This method uses a backward difference in time and a central difference in space. It requires solving a system of linear equations at each time step.
*   **Theta Scheme**: This method is a generalization of the explicit and implicit methods, using a weighted average of the two.

## Files

*   `Explicit Scheme.gif`: Animation of the explicit scheme.
*   `Implicit Scheme.gif`: Animation of the implicit scheme.
*   `theta_ Scheme_animation.gif`: Animation of the theta scheme.
*   `Explicit Finite Difference Method for Heat Equation`: Python script for simulating the heat equation using the explicit finite difference method.
*   `Implicit Scheme: Heat Equation Simulation`: Python script for simulating the heat equation using the implicit finite difference method.
*   `Theta Scheme Simulation of the Heat Equation`: Python script for simulating the heat equation using the theta scheme.

## Parameters

The scripts use the following parameters:

*   `L`: Length of the spatial domain.
*   `n`: Number of spatial intervals.
*   `x`: Spatial discretization.
*   `h`: Spatial step size.
*   `T`: Total time.
*   `tau`: Time step.
*   `m`: Number of time steps.
*   `gam`:  Constant, defined as \$\\tau / h^2\$.
*   `theta`: Weighting parameter for the Theta Scheme.

## Initial and Boundary Conditions

*   **Initial Condition**: The initial condition is set to zero.
*   **Boundary Conditions**: Dirichlet boundary conditions are applied, setting the solution to zero at both ends of the spatial domain.

## Dependencies

The following Python libraries are required:

*   numpy
*   matplotlib
*   scipy

You can install these dependencies using pip:

```bash
pip install numpy matplotlib scipy
```

## Usage

To run the simulations, execute the Python scripts. For example:

```bash
python "Explicit Scheme.py"
python "Implicit Scheme.py"
python "Theta Scheme.py"
```

## True Solution and Error Calculation

The code calculates the absolute error between the numerical solution and the true solution, which is given by  \$u(t, x) = t \\cdot sin(\\pi x)\$. The maximum absolute error is printed for each time step.

## Animations

The simulations generate animations (GIF files) showing the numerical and true solutions over time.

## Notes

*   The explicit scheme is subject to stability constraints, which may require adjusting the time step (`tau`) to ensure stable results.
*   The implicit scheme is unconditionally stable but requires solving a system of linear equations at each time step.
*   The theta scheme allows you to interpolate between the explicit and implicit methods.
```
