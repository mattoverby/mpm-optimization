// Copyright (c) 2016 Matt Overby
//
// MPM-OPTIMIZATION Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of
// conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list
// of conditions and the following disclaimer in the documentation and/or other materials
// provided with the distribution.
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef MPM_BASE_HPP
#define MPM_BASE_HPP 1

#include "Particle.hpp"

namespace mpm
{

// Grid node base class
class GridNode
{
  public:
    GridNode(int global_index)
        : v(0, 0, 0)
        , m(0)
        , global_idx(global_index)
        , active_idx(0)
    {
    }
    virtual ~GridNode() {}

    double          m;           // mass
    Eigen::Vector3d v;           // velocity
    unsigned int    global_idx;  // index into the global buffer
    unsigned int active_idx;  // index into the active buffer (set by MPM::p2g_mass)

    std::vector<Particle*> particles;  // particles within interp radius (set by MPM::p2g_mass)
    std::vector<double> wip;  // corresponding particle weight (set by MPM::p2g_mass)
    std::vector<Eigen::Vector3d> dwip;  // corresponding particle deriv. weight (set by MPM::p2g_mass)
};


// MPM is a static class used by different integrators. It has basic functions
// for steps 1,2,7,8 of the SIGGRAPH course notes.
class MPM
{
  public:
    static Eigen::Vector3d cellsize;
    static Eigen::Vector3i gridsize;
    static double          timestep_s;

    //
    //	Particle to Grid functions
    //

    static void clear_grid(std::vector<GridNode*>& active_grid);
    static void p2g_mass(std::vector<Particle*>& particles,
                         std::vector<GridNode*>& grid,
                         std::vector<GridNode*>& active_grid);
    static void init_volumes(std::vector<Particle*>&       particles,
                             const std::vector<GridNode*>& grid);
    static void p2g_velocity(const std::vector<Particle*>& particles,
                             std::vector<GridNode*>&       active_grid);

    //
    //	Grid to Particle functions
    //

    static void g2p_deformation(std::vector<Particle*>&       particles,
                                const std::vector<GridNode*>& grid);
    static void g2p_velocity(std::vector<Particle*>&       particles,
                             const std::vector<GridNode*>& grid);
    static void particle_advection(std::vector<Particle*>& particles);

    //
    //	Collision
    //

    static void grid_collision(std::vector<GridNode*>& active_grid);

    //
    //	Interpolation functions (see Interp.hpp)
    //

    // Pass in (normalized) grid index and returns its index in the grid vector.
    static int grid_index(const Eigen::Vector3i& pos);

    // Inverse of grid_index (position in grid, NOT world space!)
    static Eigen::Vector3i grid_pos(int grid_idx);

    // Makes an easy-to-iterate list of grid nodes that are within interpolation range of a particle.
    static void make_grid_loop(const Particle*               p,
                               std::vector<int>&             grid_idx,
                               std::vector<double>&          Wip,
                               std::vector<Eigen::Vector3d>& dWip);

    // Computes grid position minus particle position (world space) for APIC
    static Eigen::Vector3d xi_xp(const int grid_idx, const Particle* p);

};  // end class MPM

}  // end namespace mpm


#endif
