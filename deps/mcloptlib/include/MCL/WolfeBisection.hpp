// The MIT License (MIT)
// Copyright (c) 2018 Matt Overby
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

#ifndef MCL_WOLFEBISECTION_H
#define MCL_WOLFEBISECTION_H

#include "Problem.hpp"

namespace mcl {
namespace optlib {

// Bisection method for Weak Wolfe conditions
template<typename Scalar, int DIM>
class WolfeBisection {
private:
	// Strong wolfe conditions: the armijo rule and a stronger curvature condition
	// alpha = step length
	// fx_ap = f(x + alpha p)
	// fx = f(x)
	// pT_gx = p^T ( grad f(x) )
	// pT_gx_ap = p^T ( grad f(x + alpha p) )
	static inline bool strong_wolfe( Scalar alpha,
		Scalar fx, Scalar pT_gx, 
		Scalar fx_ap, Scalar pT_gx_ap,
		Scalar wolfe_c1, Scalar wolfe_c2 ){
		if( !(fx_ap <= fx + wolfe_c1 * alpha * pT_gx) ){ return false; } // armijo rule
		return std::abs( pT_gx_ap ) <= wolfe_c2 * std::abs( pT_gx );
	}

public:
	typedef Eigen::Matrix<Scalar,DIM,1> VecX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatX;

	static inline Scalar search(int verbose, int max_iters, const VecX &x, const VecX &p, Problem<Scalar,DIM> &problem, Scalar alpha0) {

		const Scalar t_eps = std::numeric_limits<Scalar>::epsilon();
		const Scalar wolfe_c1 = 0.0001;
		const Scalar wolfe_c2 = 0.8; // should be 0.1 for CG!
		const int dim = x.rows();
		double alpha = alpha0;
		double alpha_min = 1e-8;
		double alpha_max = 1;

		VecX grad0, grad_new;
		if( DIM == Eigen::Dynamic ){
			grad0 = VecX::Zero(dim);
			grad_new = VecX::Zero(dim);
		}
		Scalar fx0 = problem.gradient(x, grad0);
		const Scalar gtp = grad0.dot(p);
		bool min_set = false;

		int iter = 0;
		for( ; iter < max_iters; ++iter ){

			// Should we stop iterating?
			if( std::abs(alpha_max-alpha_min) <= t_eps ){ break; }

			// Step halfway
			alpha = ( alpha_max + alpha_min ) * 0.5;
			grad_new.setZero();
			Scalar fx_ap = problem.gradient(x + alpha*p, grad_new);
			Scalar gt_ap = grad_new.dot( p );

			// Check the wolfe conditions
			bool happy_wolfe = strong_wolfe( alpha, fx0, gtp, fx_ap, gt_ap, wolfe_c1, wolfe_c2 );
			if( happy_wolfe ){
				alpha_min = alpha;
				min_set = true;
			}
			else { alpha_max = alpha; }

		} // end bs iters

		if( iter == max_iters ){
			if( verbose > 0 ){ printf("WolfeBisection::linesearch Error: Reached max_iters\n"); }
			return -1;
		}
		if( !min_set ){
			if( verbose > 0 ){ printf("WolfeBisection::linesearch Error: LS blocked\n"); }
			return -1;
		}

		return alpha_min;

	}
};

}
}

#endif
