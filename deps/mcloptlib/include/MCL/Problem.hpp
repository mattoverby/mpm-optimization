// The MIT License (MIT)
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

#ifndef MCL_PROBLEM_H
#define MCL_PROBLEM_H

#if MCL_DEBUG == 1
#include <iostream>
#endif

#include <Eigen/Dense>

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM>
class Problem {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VecX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatX;

public:
	// Returns true if the solver has converged
	// x0 is the result of the previous iteration
	// x1 is the result at the current iteration
	// grad is the gradient at the last iteration
	virtual bool converged(const VecX &x0, const VecX &x1, const VecX &grad) = 0;

	// Compute just the value
	virtual Scalar value(const VecX &x) = 0;

	// Compute the objective value and the gradient
	virtual Scalar gradient(const VecX &x, VecX &grad){
		finiteGradient(x, grad);
		return value(x);
	}

	// Compute hessian
	virtual void hessian(const VecX &x, MatX &hessian){
		finiteHessian(x, hessian);
	}

	// Solve dx = H^-1 -g (used by Newton's)
	virtual void solve_hessian(const VecX &x, const VecX &grad, VecX &dx){
		MatX hess;
		if( DIM  == Eigen::Dynamic ){
			int dim = x.rows();
			hess = MatX::Zero(dim,dim);
		}
		hessian(x,hess); // hessian at x_n

		// Going with with high-accurate, low requirements as default factorization for lin-solve
		// Copied from https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
		//	Method			Requirements	Spd (sm)	Spd (lg)	Accuracy
		//	partialPivLu()		Invertible	++		++		+
		//	fullPivLu()		None		-		- -		+++
		//	householderQr()		None		++		++		+
		//	colPivHouseholderQr()	None		++		-		+++
		//	fullPivHouseholderQr()	None		-		- -		+++
		//	llt()			PD		+++		+++		+
		//	ldlt()			P/N SD 		+++		+		++
		if( DIM == Eigen::Dynamic || DIM > 4 ){
			dx = hess.fullPivLu().solve(-grad);
		}
		else {
			dx = -hess.inverse()*grad;
		}
	}

	// Gradient with finite differences
	inline void finiteGradient(const VecX &x, VecX &grad){
		const int accuracy = 0; // accuracy can be 0, 1, 2, 3
		const Scalar eps = 2.2204e-6;
		const std::vector< std::vector<Scalar> > coeff =
		{ {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} };
		const std::vector< std::vector<Scalar> > coeff2 =
		{ {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} };
		const std::vector<Scalar> dd = {2, 12, 60, 840};
		int dim = x.rows();
		if( grad.rows() != dim ){ grad = VecX::Zero(dim); }
		else{ grad.setZero(); }
		for(int d = 0; d < dim; ++d){
			for (int s = 0; s < 2*(accuracy+1); ++s){
				VecX xx = x.eval();
				xx[d] += coeff2[accuracy][s]*eps;
				grad[d] += coeff[accuracy][s]*value(xx);
			}
			grad[d] /= (dd[accuracy]* eps);
		}
	} // end finite grad

	// Hessian with finite differences
	inline void finiteHessian(const VecX &x, MatX &hess){
		const Scalar eps = std::numeric_limits<Scalar>::epsilon()*10e7;
		int dim = x.rows();
		if( hess.rows() != dim || hess.cols() != dim ){ hess = MatX::Zero(dim,dim); }
		for(int i = 0; i < dim; ++i){
			for(int j = 0; j < dim; ++j){
				VecX xx = x;
				Scalar f4 = value(xx);
				xx[i] += eps;
				xx[j] += eps;
				Scalar f1 = value(xx);
				xx[j] -= eps;
				Scalar f2 = value(xx);
				xx[j] += eps;
				xx[i] -= eps;
				Scalar f3 = value(xx);
				hess(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);
			}
		}
	} // end finite hess
};

}
}

#endif
