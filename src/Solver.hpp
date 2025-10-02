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

#ifndef MPM_SOLVER_HPP
#define MPM_SOLVER_HPP 1

#include "MPM.hpp"
#include "MCL/LBFGS.hpp"

namespace mpm
{

class Solver
{
  public:
    Solver()
        : elapsed_s(0)
    {
    }
    ~Solver();

    //
    //  Data
    //

    std::vector<GridNode*> m_grid;
    std::vector<GridNode*> active_grid;  // resized each time step
    std::vector<Particle*> m_particles;
    void                   get_vertices(Eigen::MatrixXd& X) const;
    mcl::optlib::LBFGS<double, Eigen::Dynamic> optimizer;

    //
    //  Settings
    //

    double                elapsed_s;
    const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -9.8, 0.0);

    //
    //  Called by the gui
    //

    // Initialize the system
    bool initialize();

    // Run a timestep.
    // Returns true on success.
    bool step(float screen_dt);

  private:
    double compute_timestep(double screen_dt);

    //
    //	Integration
    //

    void explicit_solve();
    void implicit_solve();

};  // end class system

}  // end namespace mpm

#endif
