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

#include "Solver.hpp"
#include "LBFGS.hpp"
#include <igl/parallel_for.h>
#include <iostream>

namespace mpm
{

bool Solver::initialize()
{
    // Fill a sphere with particles
    Eigen::Vector3d sphere_center(3.f, 3.f, 3.f);
    double          sphere_rad = 1.0;
    double          step       = 0.2;
    double          p_mass     = 0.01;
    for(double x = sphere_center[0] - sphere_rad; x < sphere_center[0] + sphere_rad; x += step)
    {
        for(double y = sphere_center[1] - sphere_rad; y < sphere_center[1] + sphere_rad; y += step)
        {
            for(double z = sphere_center[2] - sphere_rad; z < sphere_center[2] + sphere_rad; z += step)
            {
                Eigen::Vector3d pos(x, y, z);
                if((pos - sphere_center).norm() < sphere_rad)
                {
                    Particle* p = new pNeoHookean();
                    p->x        = pos;
                    p->m        = p_mass;
                    particles.push_back(p);
                }
            }
        }
    }

    // The boundaries are determined by these variables, and assumed positive.
    // I.e. the boundary is 0 to gridsize*cellsize.
    mpm_settings.cellsize = Eigen::Vector3d(.5, .5, .5);
    mpm_settings.gridsize = Eigen::Vector3i(16, 16, 16);

    for(int i = 0; i < mpm_settings.gridsize[0] * mpm_settings.gridsize[1] * mpm_settings.gridsize[2]; ++i)
    {
        grid.push_back(new GridNode(i));
    }

    std::cout << "Gridsize: " << mpm_settings.gridsize[0] << ' ' << mpm_settings.gridsize[1] << ' '
              << mpm_settings.gridsize[2] << std::endl;
    std::cout << "Cellsize: " << mpm_settings.cellsize[0] << ' ' << mpm_settings.cellsize[1] << ' '
              << mpm_settings.cellsize[2] << std::endl;
    std::cout << "Num Particles: " << particles.size() << std::endl;
    std::cout << "Num Nodes: " << grid.size() << std::endl;


    return true;
}


Solver::~Solver()
{
    for(size_t i = 0; i < particles.size(); ++i)
    {
        delete particles[i];
    }
    particles.clear();

    for(size_t i = 0; i < grid.size(); ++i)
    {
        delete grid[i];
    }
    grid.clear();

    active_grid.clear();
}

void Solver::get_vertices(Eigen::MatrixXd& X) const
{
    int np = particles.size();
    if(np == 0)
    {
        return;
    }

    if(X.rows() < np)
    {
        X.resize(np, 3);
    }

    for(int i = 0; i < np; ++i)
    {
        X.row(i) = particles[i]->x;
    }
}

void Solver::step(float screen_dt)
{
    // Compute timestep
    mpm_settings.timestep_s = compute_timestep(screen_dt);
    MPM::clear_grid(active_grid);

    // Particle to grid
    MPM::p2g_mass(mpm_settings, particles, grid, active_grid);  // step 1
    if(elapsed_s <= 0.f)
    {
        MPM::init_volumes(mpm_settings, particles, grid);
    }
    MPM::p2g_velocity(mpm_settings, particles, active_grid);  // step 1

    // Steps 4, 5 in course notes
    //	explicit_solve();
    implicit_solve();
    MPM::grid_collision(mpm_settings, active_grid);

    // Grid to particle, and advection
    // Steps 6, 7, 8 in course notes
    MPM::g2p_deformation(mpm_settings, particles, grid);
    MPM::g2p_velocity(mpm_settings, particles, grid);
    MPM::particle_advection(mpm_settings, particles);

    elapsed_s += mpm_settings.timestep_s;

}  // end step


void Solver::explicit_solve()
{
    igl::parallel_for(
        particles.size(),
        [&](int i)
        {
            Particle*       p   = particles[i];
            Eigen::Matrix3d p_F = p->get_deform_grad();
            p->tempP            = p->vol * p->get_piola_stress(p_F) * p->Fe.transpose();
        },
        MPM_MIN_PARALLEL);


    // Velocity update
    igl::parallel_for(
        active_grid.size(),
        [&](int i)
        {
            GridNode*       g = active_grid[i];
            Eigen::Vector3d f(0, 0, 0);
            for(size_t j = 0; j < g->particles.size(); ++j)
            {
                Particle*       p   = g->particles[j];
                Eigen::Matrix3d p_F = p->get_deform_grad();
                f -= p->tempP * g->dwip[j];
            }

            g->v += (mpm_settings.timestep_s * f) / g->m + mpm_settings.timestep_s * gravity;
        },
        MPM_MIN_PARALLEL);

}  // end compute forces


double Solver::gradient(const Eigen::VectorXd& v, Eigen::VectorXd& grad)
{
    using namespace Eigen;
    VectorXd energy_sum = VectorXd::Zero(std::max(particles.size(), active_grid.size()));

    //
    //	Compute Particle Deformation Gradients for new grid velocities
    //
    igl::parallel_for(
        particles.size(),
        [&](int i)
        {
            Particle* p   = particles[i];
            Matrix3d  p_F = p->get_deform_grad();

            Matrix3d vel_grad;
            vel_grad.fill(0.f);

            std::vector<int>      grid_nodes;
            std::vector<double>   wip;
            std::vector<Vector3d> dwip;
            MPM::make_grid_loop(mpm_settings, p, grid_nodes, wip, dwip);
            for(size_t j = 0; j < grid_nodes.size(); ++j)
            {
                int      idx = grid[grid_nodes[j]]->active_idx;
                Vector3d curr_v(v[idx * 3 + 0], v[idx * 3 + 1], v[idx * 3 + 2]);
                vel_grad += (mpm_settings.timestep_s * curr_v) * dwip[j].transpose();
            }

            Matrix3d newF = (Matrix3d::Identity() + vel_grad) * p_F;  // eq 193 course notes
            p->tempP      = p->vol * p->get_piola_stress(newF) * p_F.transpose();

            // Update the total energy
            energy_sum[i] = p->get_energy_density(newF) * p->vol;
        },
        MPM_MIN_PARALLEL);


    //
    //	Compute energy gradient
    //
    igl::parallel_for(
        active_grid.size(),
        [&](int i)
        {
            GridNode* node = active_grid[i];
            Vector3d  curr_v(v[i * 3 + 0], v[i * 3 + 1], v[i * 3 + 2]);

            // Only compute gradient if requested (grad is sized)
            if(grad.rows() == v.rows())
            {
                Vector3d momentum_grad = node->m * (curr_v - node->v);
                Vector3d energy_grad(0, 0, 0);
                for(size_t j = 0; j < node->particles.size(); ++j)
                {
                    energy_grad += node->particles[j]->tempP * node->dwip[j] * mpm_settings.timestep_s;
                }

                grad.segment(i * 3, 3) = momentum_grad + energy_grad;
            }
            energy_sum[i] += 0.5 * node->m * (curr_v - node->v).squaredNorm();
        },
        MPM_MIN_PARALLEL);

    return energy_sum.sum();
}


void Solver::implicit_solve()
{
    using namespace Eigen;

    // Initial guess of grid velocities
    VectorXd v = VectorXd::Zero(active_grid.size() * 3);
    igl::parallel_for(
        active_grid.size(),
        [&](int i)
        {
            for(int j = 0; j < 3; ++j)
            {
                v[i * 3 + j] = active_grid[i]->v[j];
            }
        },
        MPM_MIN_PARALLEL);

    // Minimize
    mcl::LBFGS<VectorXd> optimizer;

    optimizer.gradient = [&](const VectorXd& x, VectorXd& g) -> double
    {
        double objective = gradient(x, g);
        return objective;
    };

    optimizer.minimize(v);

    // Copy solver results back to grid
    igl::parallel_for(
        active_grid.size(),
        [&](int i)
        {
            for(int j = 0; j < 3; ++j)
            {
                active_grid[i]->v[j] = v[i * 3 + j];
            }
            active_grid[i]->v += mpm_settings.timestep_s * gravity;  // explicitly add gravity
        },
        MPM_MIN_PARALLEL);

}  // end implicit solve


double Solver::compute_timestep(double screen_dt)
{
    using namespace Eigen;
    double   max_timestep = 0.01;
    VectorXd squared_vel  = VectorXd::Zero(particles.size());

    igl::parallel_for(
        particles.size(), [&](int i) { squared_vel[i] = particles[i]->v.squaredNorm(); }, MPM_MIN_PARALLEL);

    double max_v = std::sqrt(squared_vel.maxCoeff());
    if(max_v > 1e-8)
    {
        double min_cellsize = std::min(mpm_settings.cellsize[0],
                                       std::min(mpm_settings.cellsize[1], mpm_settings.cellsize[2]));
        double ts_adjust    = 0.04;
        double dt           = ts_adjust * min_cellsize / sqrt(max_v);
        if(dt > max_timestep)
        {
            return max_timestep;
        }
        if(dt > screen_dt)
        {
            return screen_dt;
        }
        return dt;
    }

    return max_timestep;

}  // end adaptive timestep

}  // namespace mpm
