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
// ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef MPM_SOLVER_HPP
#define MPM_SOLVER_HPP 1

#include "MPM.hpp"

namespace mpm
{

class Solver
{
  public:
    Solver()
        : gravity(0.0, -9.8, 0.0)
        , elapsed_s(0)
    {
    }
    ~Solver();

    std::vector<GridNode*> grid;
    std::vector<Particle*> particles;
    std::vector<GridNode*> active_grid;  // resized each time step
    MPM::Settings          mpm_settings;
    const Eigen::Vector3d  gravity;
    double                 elapsed_s;

    // Initialize the system
    bool initialize();

    // Run a timestep.
    void step(float screen_dt);

    // Gradient of the objective.
    double gradient(const Eigen::VectorXd& v, Eigen::VectorXd& grad);

    // Retuns current particle positions
    void get_vertices(Eigen::MatrixXd& X) const;

  private:
    double compute_timestep(double screen_dt);
    void   explicit_solve();
    void   implicit_solve();

};  // end class solver

}  // end namespace mpm

#endif
