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

#ifndef INTERP_HPP
#define INTERP_HPP 1

#include <Eigen/Dense>
#include <vector>

namespace mpm
{

// Cubic B-spline
static inline double cspline(double x)
{
    if(x == 0.0)
    {
        x = 1e-12;
    }
    x = fabs(x);
    if(x < 1.0)
    {
        return x * x * (x * 0.5 - 1.0) + 2.0 / 3.0;
    }
    else if(x < 2.0)
    {
        return x * (x * (-x / 6.0 + 1.0) - 2.0) + 4.0 / 3.0;
    }
    return 0.0;
}

// Slope of cubic spline
static inline double d_cspline(double x)
{
    if(x == 0.0)
    {
        x = 1e-12;
    }
    double abs_x = fabs(x);
    if(abs_x < 1.0)
    {
        return 1.5 * x * abs_x - 2.0 * x;
    }
    else if(x < 2.0)
    {
        return -x * abs_x * 0.5 + 2.0 * x - 2.0 * x / abs_x;
    }
    return 0.0;
}

// Quadratic spline
static inline double qspline(double x)
{
    double fx = fabs(x);
    if(fx < 0.5)
    {
        return (0.75 - x * x);
    }
    else if(fx < 1.5)
    {
        return (0.5 * x * x - 1.5 * fx + 9.0 / 8.0);
    }
    return 0.0;
}

// Slope of quadratic spline
static inline double d_qspline(double x)
{
    double fx = fabs(x);
    if(fx < 0.5)
    {
        return (-2.0 * x);
    }
    else if(fx < 1.5)
    {
        if(x < 0.0)
        {
            return fabs(fx - 1.5);
        }
        else
        {
            return (fx - 1.5);
        }
    }
    return 0.0;
}

static inline bool isreal(double n)
{
    return !std::isnan(n) && !std::isinf(n);
}

}  // end namespace mpm

#endif
