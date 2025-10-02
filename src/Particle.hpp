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

#ifndef MPM_PARTICLE_HPP
#define MPM_PARTICLE_HPP 1

#include <iostream>
#include <memory>
#include "Interp.hpp"

namespace mpm
{

// Projection, Singular Values, SVD's U, SVD's V transpose
static inline void projected_svd(const Eigen::Matrix3d& F,
                                 Eigen::Vector3d&       S,
                                 Eigen::Matrix3d&       U,
                                 Eigen::Matrix3d&       Vt)
{
    using namespace Eigen;
    JacobiSVD<Matrix3d> svd(F, ComputeFullU | ComputeFullV);
    S          = svd.singularValues();
    U          = svd.matrixU();
    Vt         = svd.matrixV().transpose();
    Matrix3d J = Matrix3d::Identity();
    J(2, 2)    = -1.0;
    // Check for inversion
    if(U.determinant() < 0.0)
    {
        U = U * J;
        S[2] *= -1.0;
    }
    if(Vt.determinant() < 0.0)
    {
        Vt = J * Vt;
        S[2] *= -1.0;
    }

}  // end projected svd


// Particle (material point) class
class Particle
{
  public:
    Particle()
        : m(0.1)
        , vol(0.0)
        , v(0, 0, 0)
        , x(0, 0, 0)
    {
        Fe    = Eigen::Matrix3d::Identity();
        tempP = Eigen::Matrix3d::Identity();
        B.setZero();
        D.setZero();
    }

    virtual Eigen::Matrix3d get_piola_stress(Eigen::Matrix3d& currF) const  = 0;
    virtual double          get_energy_density(Eigen::Matrix3d& newF) const = 0;
    virtual Eigen::Matrix3d get_deform_grad() { return Fe; }
    virtual void update_deform_grad(Eigen::Matrix3d velocity_grad, double timestep_s)
    {
        Fe = (Eigen::Matrix3d::Identity() + timestep_s * velocity_grad) * Fe;
    }

    double          m;    // mass
    double          vol;  // rest volume, computed at first timestep
    Eigen::Vector3d v;    // velocity
    Eigen::Vector3d x;    // location
    Eigen::Matrix3d B;    // affine matrix (eq 176 course notes)
    Eigen::Matrix3d D;  // helper matrix, set in p_to_g mass (eq 174 course notes)

    Eigen::Matrix3d tempP;  // temporary piola stress tensor computed by solver
    Eigen::Matrix3d Fe;     // elastic deformation gradient
};

// NeoHookean Particle
class pNeoHookean : public Particle
{
  public:
    pNeoHookean()
        : mu(10)
        , lambda(10)
    {
    }
    double mu;
    double lambda;

    inline Eigen::Matrix3d get_piola_stress(Eigen::Matrix3d& currF) const
    {
        // Fix inversions:
        Eigen::Matrix3d U, Vt, Ftemp;
        Eigen::Vector3d S;
        projected_svd(currF, S, U, Vt);
        if(S[2] < 0.0)
        {
            S[2] *= -1.0;
        }
        Ftemp = U * S.asDiagonal() * Vt;

        // Compute Piola stress tensor:
        double J = Ftemp.determinant();
        assert(isreal(J));
        assert(J > 0.0);
        Eigen::Matrix3d Fit = (Ftemp.inverse()).transpose();  // F^(-T)
        Eigen::Matrix3d P   = mu * (Ftemp - Fit) + lambda * (log(J) * Fit);
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                assert(isreal(P(i, j)));
            }
        }
        return P;

    }  // end compute piola stress tensor

    inline double get_energy_density(Eigen::Matrix3d& newF) const
    {
        // Fix inversions:
        Eigen::Matrix3d U, Vt, Ftemp;
        Eigen::Vector3d S;
        projected_svd(newF, S, U, Vt);
        if(S[2] < 0.0)
        {
            S[2] *= -1.0;
        }
        Ftemp = U * S.asDiagonal() * Vt;

        // Compute energy density:
        double J = Ftemp.determinant();
        assert(isreal(J));
        assert(J > 0.0);
        double t1 = 0.5 * mu * ((Ftemp.transpose() * Ftemp).trace() - 3.0);
        double t2 = -mu * log(J);
        double t3 = 0.5 * lambda * (log(J) * log(J));
        assert(isreal(t1));
        assert(isreal(t2));
        assert(isreal(t3));
        return (t1 + t2 + t3);

    }  // end compute energy density

};  // end class neohookean particle

}  // end namespace mpm

#endif
