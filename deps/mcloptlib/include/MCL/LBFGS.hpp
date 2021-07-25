// The MIT License (MIT)
// Copyright (c) 2017 University of Minnesota
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

#ifndef MCL_LBFGS_H
#define MCL_LBFGS_H

#include "Minimizer.hpp"

namespace mcl {
namespace optlib {

// L-BFGS implementation based on Nocedal & Wright Numerical Optimization book (Section 7.2)
// DIM = dimension of the problem
// M = history window
//
// Original Author: Ioannis Karamouzas
//
template<typename Scalar, int DIM, int M=8>
class LBFGS : public Minimizer<Scalar,DIM> {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VecX;
	typedef Eigen::Matrix<Scalar,DIM,M> MatM;
	typedef Eigen::Matrix<Scalar,M,1> VecM;

public:
	bool show_denom_warning; // Print out warning for zero denominators

	LBFGS() : show_denom_warning(false) {
		this->m_settings.max_iters = 50;
		show_denom_warning = this->m_settings.verbose > 0 ? true : false;
	}

	// Returns number of iterations used
	int minimize(Problem<Scalar,DIM> &problem, VecX &x){

		MatM s, y;
		VecM alpha, rho;
		VecX grad, q, grad_old, x_old, x_last;

		if( DIM==Eigen::Dynamic ){
			int dim = x.rows();
			s = MatM::Zero(dim,M);
			y = MatM::Zero(dim,M);
			alpha = VecM::Zero(M);
			rho = VecM::Zero(M);
			grad = VecX::Zero(dim);
			q = VecX::Zero(dim);
			grad_old = VecX::Zero(dim);
			x_old = VecX::Zero(dim);
			x_last = VecX::Zero(dim);
		}

		problem.gradient(x, grad);

		Scalar gamma_k = 1.0;
		Scalar alpha_init = 1.0;

		int global_iter = 0;
		int max_iters = this->m_settings.max_iters;
		int verbose = this->m_settings.verbose;

		for( int k=0; k<max_iters; ++k ){

			x_old = x;
			grad_old = grad;
			q = grad;
			global_iter++;
	
			// L-BFGS first - loop recursion		
			int iter = std::min(M, k);
			for(int i = iter - 1; i >= 0; --i){
				rho(i) = 1.0 / ((s.col(i)).dot(y.col(i)));
				alpha(i) = rho(i)*(s.col(i)).dot(q);
				q = q - alpha(i)*y.col(i);
			}

			// L-BFGS second - loop recursion			
			q = gamma_k*q;
			for(int i = 0; i < iter; ++i){
				Scalar beta = rho(i)*q.dot(y.col(i));
				q = q + (alpha(i) - beta)*s.col(i);
			}

			// is there a descent
			Scalar dir = q.dot(grad);
			if(dir <= 0 ){
				q = grad;
				max_iters -= k;
				k = 0;
				alpha_init = std::min(1.0, 1.0 / grad.template lpNorm<Eigen::Infinity>() );
			}

			Scalar rate = this->linesearch(x, -q, problem, alpha_init);

			if( rate <= 0 ){
				if( verbose > 0 ){ printf("LBFGS::minimize: Failure in linesearch\n"); }
				return Minimizer<Scalar,DIM>::FAILURE;
			}

			x_last = x;
			x -= rate * q;
			if( problem.converged(x_last,x,grad) ){ break; }

			problem.gradient(x,grad);
			VecX s_temp = x - x_old;
			VecX y_temp = grad - grad_old;

			// update the history
			if(k < M){
				s.col(k) = s_temp;
				y.col(k) = y_temp;
			}
			else {
				s.leftCols(M - 1) = s.rightCols(M - 1).eval();
				s.rightCols(1) = s_temp;
				y.leftCols(M - 1) = y.rightCols(M - 1).eval();
				y.rightCols(1) = y_temp;
			}
		
			Scalar denom = y_temp.dot(y_temp);
			if( std::abs(denom) <= 0 ){
				if( show_denom_warning ){
					printf("LBFGS::minimize Warning: Encountered a zero denominator\n");
				}
				break;
			}
			gamma_k = s_temp.dot(y_temp) / denom;
			alpha_init = 1.0;

		}

		return global_iter;

	} // end minimize
};

}
}

#endif
