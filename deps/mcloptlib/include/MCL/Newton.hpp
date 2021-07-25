// The MIT License (MIT)
// Copyright (c) 2017 Matt Overby
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef MCL_NEWTON_H
#define MCL_NEWTON_H

#include "Minimizer.hpp"

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM>
class Newton : public Minimizer<Scalar,DIM> {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VectorX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatrixX;

public:
	Newton() {
		this->m_settings.max_iters = 20;
	}

	int minimize(Problem<Scalar,DIM> &problem, VectorX &x){

		VectorX grad, delta_x, x_last;
		if( DIM  == Eigen::Dynamic ){
			int dim = x.rows();
			x_last.resize(dim);
			grad.resize(dim);
			delta_x.resize(dim);
		}

		int verbose = this->m_settings.verbose;
		int max_iters = this->m_settings.max_iters;
		int iter = 0;
		for( ; iter < max_iters; ++iter ){

			problem.gradient(x,grad);
			problem.solve_hessian(x,grad,delta_x);

			Scalar rate = this->linesearch(x, delta_x, problem, 1.0);

			if( rate <= 0 ){
				if( verbose > 0 ){ printf("Newton::minimize: Failure in linesearch\n"); }
				return Minimizer<Scalar,DIM>::FAILURE;
			}

			x_last = x;
			x += rate * delta_x;
			if( problem.converged(x_last,x,grad) ){ break; }
		}

		return iter;
	}

};

} // ns optlib
} // ns mcl

#endif
