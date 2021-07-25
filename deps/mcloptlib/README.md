# mcloptlib

By Matt Overby  
[http://www.mattoverby.net](http://www.mattoverby.net)

mcloptlib is a header-only optimization library for C++ using Eigen and is geared towards lower-dimension graphics problems.
Originally a fork of [Patrick Wieschollek's CppOptimizationLibrary](https://github.com/PatWie/CppNumericalSolvers), but has diverged considerably.

## Contents:

Optimization algorithms:
- Newton's
- Non-linear conjugate gradient
- L-BFGS
- Trust Region with
  - Cauchy Point
  - Dog Leg

Linesearch methods:
- Backtracking (Armijo)
- Backtracking with cubic interpolation
- Bisection
- MoreThuente

## To-do:

- Option of std::function for value/gradient instead of Problem class
- Sparse Hessians
- Auto-diff
