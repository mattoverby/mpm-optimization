# mpm-optimiziation

***2021 Update***: I've copied the original blog post to this readme, updated code formatting, added rendering with libigl, and switched to [mcloptlib](https://github.com/mattoverby/mcloptlib)'s LBFGS.
I wrote this in 2016 as an exercise, so please don't hold it against me :smile:.

Optimization-based implicit material point method with APIC,
implemented from: [SIGGRAPH course notes](http://web.cs.ucla.edu/~cffjiang/mpmcourse/mpm_course_nodes.pdf),
[Optimization Integrator for Large Time Steps](https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf),
and [The Affine Particle-In-Cell Method](https://disney-animation.s3.amazonaws.com/uploads/production/publication_asset/104/asset/apic-aselle-final.pdf).

## Optimization-Based Material Point Method

The Material Point Method (MPM) is a powerful technique for physics simulation.
MPM is great for continuum substances like fluids and soft bodies, as well as for dealing with complex problems like fracture and self collisions.
It has received a lot of attention in computer graphics, spurred by the snow simulation and phase change/heat transfer papers.

Some time around 2016 I coded up MPM to test out some research ideas, and wrote a little text to go along with it.
My ideas didn't work out, but I put this code up on github instead of letting it collect metaphorical dust on my hard drive.
It is not efficient nor bug free.
The demo is nothing more than a big neo-Hookean blob dropped on a ground plane.

### Overview

The basic idea is that particles move around the domain and carry mass.
When it's time to do time integration, the particle mass/momentum is mapped to a grid.
New velocities are computed on the grid, then mapped back to the particles, which are then advected.
The grid essentially acts as a scratch pad for simulation that is cleared and reset each time step.
For a (better) background on MPM, read the introduction chapter from the [SIGGRAPH course notes](http://alexey.stomakhin.com/research/siggraph2016_mpm.pdf).

### Time integration

The typical starting point for computing new velocities is explicit time integration.
Though explicit integration is fast to compute and easy to understand, it has problems with stability and overshooting.
If we take too large of a time step, the particle goes past the physically-correct trajectory.
Intermediate techniques like midpoint method can alleviate such problems, but the silver bullet is implicit integration.

With implicit integration we figure out what the velocity needs to be to make the particles to end up where they should.
The hard part is that we now have both unknowns in the velocity calculation, which requires solving a system of nonlinear equations.
This is computationally expensive and time consuming, especially for large and complex systems.

Fortunately, many important forces (like elasticity) are conservative.
Mathematically, that means we can represent the forces as the negative gradient of potential energy.
By representing our forces as energy potentials, we can do implicit time integration by minimizing the total system energy.
The concept itself is not new, but was recently applied to Finite Elements/MPM in a paper by [Gast et al.](https://www.math.ucla.edu/%7Ejteran/papers/GSSJT15.pdf).
Section 6 of the paper covers the algorithm for optimization-based MPM.
One thing to keep in mind: in MPM the velocity calculations happen on grid nodes, not grid boundaries.

### PIC/FLIP versus APIC

Mapping to/from the grid is an important part of the material point method.
MPM was originally an extension of Fluid Implicit Particle (FLIP), which in turn was an extension of the Particle in Cell (PIC) method.
PIC and FLIP are the standard in fluid simulation, and were used in the above Gast paper.
Mike Seymour wrote [a good article on fxguide](https://www.fxguide.com/fxfeatured/the-science-of-fluid-sims/) that covers what PIC/FLIP is and why they are useful.
However, they suffer from a few well known limitations.
Namely, PIC causes things to "clump" together, and FLIP can be unstable.

A technique called the [affine particle in cell (APIC)](https://dl.acm.org/doi/10.1145/2766996) was recently introduced to combat these limitations.
It's more stable and does a better job at conserving momentum.
In short, it's worth using. To do so, you only need to make a few subtle changes to the MPM algorithm.
The SIGGRAPH course notes explain what to change, but the Gast paper does not.

## Implementation

Everything you need to know to code up MPM can be found in the SIGGRAPH course notes.
If you've ever done fluid simulation with PIC/FLIP, it should be straight forward. If not, here are some things I found useful:

1. It seemed to help with stability to have a cell size large enough to fit 8-16 particles.
2. Keep a fully-allocated grid as the domain boundary, but do minimization on a separate "active grid" list made up of pointers to grid nodes.
3. A "make grid loop" function for particle/grid mapping: given a particle, return the indices of grid nodes within the kernel radius, as well as kernel smoothing values/derivatives. I felt this was the trickiest part of the whole method. Much of this can be preallocated if you have a lot of memory, but I found (in 3D) allocating/deallocating on the spot was helpful and easier to debug.