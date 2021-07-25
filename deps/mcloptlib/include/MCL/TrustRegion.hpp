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

#ifndef MCL_TRUSTREGION_H
#define MCL_TRUSTREGION_H

#include "Minimizer.hpp"

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM>
class TrustRegion : public Minimizer<Scalar,DIM> {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VecX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatX;

public:
	//
	// The trust region method operates on an approximate hessian (B).
	// To keep the interface simple, Problem::hessian is used to obtain
	// this approximation. It shouldn't be much of a problem, since you
	// can still use approximations for Newton's (the only other method
	// implemented that requires a hessian evaluation).
	//
	// TODO Replace B with possibly sparse representation
	//
	TrustRegion() {
		this->m_settings.max_iters = 100;
	}

	int minimize(Problem<Scalar,DIM> &problem, VecX &x){

		Scalar delta_k = 2.0; // trust region radius
		const Scalar delta_max = 8.0; // max trust region radius
		const Scalar eta = 0.125; // min reduction ratio allowed (0<eta<0.25)
		const TRMethod m = this->m_settings.tr_method;
		int max_iters = this->m_settings.max_iters;
		int verbose = this->m_settings.verbose;

		VecX grad, dx, x_last;
		MatX B; // Approximate hessian
		if( DIM  == Eigen::Dynamic ){
			int dim = x.rows();
			x_last.resize(dim); // last variable
			grad.resize(dim); // grad at k
			dx.resize(dim); // descent direction
			B.resize(dim,dim);
		}

		// Init gradient and hessian
		Scalar fxk = problem.gradient(x,grad); // gradient and objective
		problem.hessian(x,B); // get hessian (or approximation)
		problem.solve_hessian(x,grad,dx); // attempt with newtons

		int iter = 0;
		for( ; iter < max_iters; ++iter ){

			// If it's outside the trust region, pick new descent
			if( dx.norm() > delta_k ){
				eval_subproblem(m, delta_k, grad, B, dx);
			}

			Scalar fxdx = problem.value(x+dx);

			// Compute reduction ratio
			Scalar rho_k = eval_reduction(fxk, fxdx, dx, grad, B);
			Scalar dx_norm = dx.norm();

			// Update trust region radius
			if( rho_k < 0.25 ){ delta_k = 0.25*delta_k; }

			// Full step, good approximation
			else if( rho_k > 0.75 && std::abs(dx_norm-delta_k) <= 0.0 ){
				delta_k = std::min( 2.0*delta_k, delta_max );
			}

			// Take a step, otherwise need to re-eval sub problem
			if( rho_k > eta ){

				x_last = x;
				x = x + dx;

				// Only need to compute gradient and hessian
				// if x has actually changed.
				fxk = problem.gradient(x,grad); // gradient and objective
				if( problem.converged(x_last,x,grad) ){ break; }

				// I think I should improve this as to not call both hessian
				// and solve_hessian, which likely causes redundant computation.
				problem.hessian(x,B); // get hessian (or approximation)
				problem.solve_hessian(x,grad,dx); // attempt with newtons
			}

			if( std::isnan(rho_k) ){
				if( verbose ){ printf("\n**TrustRegion Error: NaN reduction"); }
				return Minimizer<Scalar,DIM>::FAILURE;
			}

		}

		return iter;

	} // end minimize

protected:

	// Assumes coeffs size 3
	static inline Scalar max_roots(Scalar a, Scalar b, Scalar c){

		Scalar d = (b*b) - (4.0*a*c);
		if( d > 0 ){
			Scalar sqrt_d = std::sqrt(d);
			Scalar r1 = (-b + sqrt_d) / (2.0*a);
			Scalar r2 = (-b - sqrt_d) / (2.0*a);
			return std::max(r1,r2);
		}
		// Should I do something with the real/imaginary parts?
		// Scalar real_part = -b/(2.0*a);
		// Scalar imaginary = std::sqrt(-d)/(2.0*a);
		throw std::runtime_error("TrustRegion Error: Problem in quadratic roots");
		return 0;
	}

	// fxk = objective at x_k
	// fxdx = objective at x_k + dx
	// dx = descent direction
	// grad_k = gradient at x_k 
	// B_k = hessian guess.
	static inline Scalar eval_reduction( Scalar fxk, Scalar fxdx, const VecX &dx,
		const VecX &grad_k, const MatX &B_k ){
		// rho = ( f(x) - f(x-dx) ) / ( model(0) - model(dx) )
		// with model = f(x) + dx^T grad + 0.5 dx^T B dx
		Scalar num = fxk - fxdx;
		Scalar denom = fxk - ( fxk + dx.dot(grad_k) + 0.5 * dx.dot( B_k * dx ) );
		return num/denom;
	}

	static inline void eval_subproblem(
		const TRMethod &m, Scalar delta_k,
		const VecX &grad, const MatX &B,
		VecX &dx ){

		Scalar gTBg = grad.dot(B*grad);

		switch( m ){

			// Generally requires a large number of outer solver iterations
			case TRMethod::CauchyPoint: {

				Scalar grad_norm = grad.norm();
				Scalar tau = 1.0;
				if( gTBg > 0.0 ){ tau = std::min( 1.0, std::pow(grad_norm,3) / (delta_k*gTBg) ); }
				dx = ( -tau * delta_k / grad_norm ) * grad;

			} break;

			// Uses steepest descent if possible
			case TRMethod::DogLeg: {

				Scalar gTg = grad.dot(grad);
				VecX dx_U = ( -gTg / gTBg ) * grad;
				Scalar dx_U_norm = dx_U.norm();

				// Use steepest descent
				if( dx_U_norm >= delta_k ){ dx = delta_k/dx_U_norm * dx_U; }
				else{
					// Compute tau and update descent
					VecX dx_C = dx - dx_U;
					Scalar dx_C_norm = dx_C.norm();
					Scalar tau = max_roots( // Ax^2 + Bx + c
						dx_C_norm*dx_C_norm,
						2.0*dx_C.dot(dx_U),
						dx_U_norm*dx_U_norm - delta_k*delta_k
					);
					dx = dx_U + tau*dx_C;
				}

			} break;

		} // end swithc method

	} // end eval sub problem

};

} // ns optlib
} // ns mcl

#endif
