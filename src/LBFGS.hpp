// Copyright Matt Overby 2021.
// Distributed under the MIT License.
// From: https://github.com/mattoverby/mclgeom

#ifndef MCL_LBFGS_HPP
#define MCL_LBFGS_HPP 1

#include <Eigen/Dense>
#include <vector>
#include <functional>

namespace mcl
{

// L-BFGS implementation based on Nocedal & Wright Numerical Optimization book (Section 7.2)
// adapted from code by Ioannis Karamouzas.
//
// Function pointers are meant to be used with lambdas, e.g.
//   LBFGS<MatrixXd> lbfgs;
//   lbfgs.gradient = [&](const MatrixXd &x, MatrixXd &g)->Scalar { return ... };
//   double obj = lbfgs.minimize(x);
//
// If the g arg is not sized, don't compute gradient. This is to avoid
// extra work when only the objective is needed (line search) and to
// avoid redundant code. Example:
//   lbfgs.gradient = [&](const MatrixXd &x, MatrixXd &g)->Scalar
//   {
//     objective = ...
//     if (g.rows() == x.rows())
//     {
//        g = ... (reuse computation for objective)
//     }
//     return objective;
//   };	
//
template <typename MatrixType>
class LBFGS
{
public:
    typedef typename MatrixType::Scalar Scalar;

	struct Options
	{
		int min_iters;
		int max_iters;
		Scalar abs_tol; // absolute tol if converged(...) not set
		Scalar rel_tol; // relative tol if converged(...) not set
		int M; // history window size
		Scalar gamma; // init Hessian = gamma * I
		Options() :
			min_iters(0),
			max_iters(100),
			abs_tol(1e-5),
			rel_tol(0),
			M(6),
			gamma(1)
			{}
	} options;

	LBFGS();

	// Output from the last call to minimize(...)
	int iters() const { return num_iters; }
	Scalar gamma() const { return gamma_k; }

	// Resizes buffers and sets to zero.
	// Called during minimize(...) ONLY if there is a change in dof or M.
	// Otherwise it's assumed you're picking up where you left off
	// on the previous call to minimize(...)
	void reset(int rows, int cols);

	// Required:
	// computes objective value and gradient
	//   obj = gradient(x, g)
	// If the g arg is not sized, don't compute gradient
	std::function<Scalar(const MatrixType&, MatrixType&)> gradient;

	// Optional:
	// Returns true if the solver should exit, default uses ||g||<abs_tol or ||g||<||x||rel_tol
	//   is_converged = converged(obj, x_prev, x, grad)
	std::function<bool(Scalar obj, const MatrixType&, const MatrixType&, const MatrixType&)> converged;

	// Optional:
	// Linesearch function, default uses bracketing weak wolfe (slow!)
	//   obj_k1 = (x, grad, descent, alpha)
	// Returns new objective value and updates both x AND gradient
	std::function<Scalar(MatrixType&, MatrixType&, const MatrixType&, Scalar&)> linesearch;

	// Optional:
	// Filter descent direction, p = B(p)
	// Otherwise p = gamma_k * p is used.
	std::function<void(MatrixType&)> filter;

	// Calls initialize(x) once and iterate(x) until converged
	Scalar minimize(MatrixType& x);

	// Initialize the solver
	void initialize(MatrixType& x);

	// Take an iteration
	// Returns objective
	Scalar iterate(MatrixType& x);

	// i.e. bisection with weak Wolfe conditions
	Scalar bracketing_weakwolfe(
		MatrixType& x,
		MatrixType& grad,
		const MatrixType& p,
		Scalar &alpha) const;

	// Used if converged not set
	// Returns true if:
	// grad.norm <= abs_tol
	// or
	// grad.norm() <= rel_tol * x.norm()
	bool default_converged(
		Scalar curr_obj,
		const MatrixType& xprev,
		const MatrixType& x,
		const MatrixType& grad) const;

	Scalar inner(
		const MatrixType &a,
		const MatrixType &b) const;

protected:
	bool initialized;
	int num_iters, max_iters, k;
	Scalar gamma_k, obj_0, obj_k;
	std::vector<MatrixType> s;
	std::vector<MatrixType> y;
	Eigen::Matrix<Scalar,Eigen::Dynamic,1> alpha;
	Eigen::Matrix<Scalar,Eigen::Dynamic,1> rho;
	MatrixType grad;
	MatrixType q;
	MatrixType descent;
	MatrixType grad_old;
	MatrixType x_old;
	MatrixType x_last;
	MatrixType s_temp;
	MatrixType y_temp;

}; // end class LBFGS

//
// Implementation
//

template <typename MatrixType>
LBFGS<MatrixType>::LBFGS() :
	initialized(false),
	num_iters(0),
	max_iters(0),
	k(0),
	gamma_k(1),
	obj_0(0),
	obj_k(0)
	{}

template <typename MatrixType>
void LBFGS<MatrixType>::reset(int rows, int cols)
{
    using namespace Eigen;
	num_iters = 0;
	max_iters = 0;
	k = 0;
	gamma_k = 1;
	obj_0 = std::numeric_limits<Scalar>::max();
	obj_k = std::numeric_limits<Scalar>::max();
	int M = options.M;
	s = std::vector<MatrixType>(M, MatrixType::Zero(rows,cols));
	y = std::vector<MatrixType>(M, MatrixType::Zero(rows,cols));
	alpha = VectorXd::Zero(M);
	rho = VectorXd::Zero(M);
	grad = MatrixType::Zero(rows, cols);
	q = MatrixType::Zero(rows, cols); // inv descent
	descent = q;
	grad_old = MatrixType::Zero(rows, cols);
	x_old = MatrixType::Zero(rows, cols);
	x_last = MatrixType::Zero(rows, cols);
	s_temp = MatrixType::Zero(rows, cols);
	y_temp = MatrixType::Zero(rows, cols);
}

// Returns number of iterations used
template <typename MatrixType>
typename LBFGS<MatrixType>::Scalar
LBFGS<MatrixType>::minimize(MatrixType &x)
{
	initialize(x);

	// Did we start at the initializer?
	if (num_iters >= options.min_iters &&
		converged(obj_k,x_last,x,grad))
	{
		num_iters = 1;
		return obj_k;
	}

	for (; k<max_iters; ++k)
	{
		iterate(x);
	} // end loop lbfgs iters

	return obj_k;

} // end minimize

template <typename MatrixType>
void LBFGS<MatrixType>::initialize(MatrixType& x)
{
	initialized = true;

	if (gradient == nullptr) {
	    throw std::runtime_error("no gradient function");
	}

	if (converged == nullptr)
	{
		using namespace std::placeholders;
		converged = std::bind(&LBFGS::default_converged, this, _1, _2, _3, _4);
	}

	if (linesearch == nullptr)
	{
		using namespace std::placeholders;
		linesearch = std::bind(&LBFGS::bracketing_weakwolfe, this, _1, _2, _3, _4);
	}

	// Resize/resize variables?
	if (alpha.rows() != options.M ||
		grad.rows() != x.rows() ||
		grad.cols() != x.cols()) {
		reset(x.rows(), x.cols());
	}

	num_iters = 0;
	max_iters = std::max(options.min_iters, options.max_iters);
	gamma_k = options.gamma;
	obj_0 = gradient(x, grad);
	obj_k = obj_0;
	k = 0;
}

template <typename MatrixType>
typename LBFGS<MatrixType>::Scalar
LBFGS<MatrixType>::iterate(MatrixType& x)
{
	if (!initialized) {
		throw std::runtime_error("not initialized");
	}

	x_old = x;
	grad_old = grad;
	q = grad;
	num_iters++;

	//
	// Two-loop recursion
	//
	{
		// L-BFGS first - loop recursion		
		int iter = std::min(options.M, k);
		for(int i = iter - 1; i >= 0; --i)
		{
			Scalar denom = inner(s[i], y[i]);
			if (std::abs(denom) <= 0.0)
			{
				rho(i) = 0;
				alpha(i) = 0;
				continue;
			}
			rho(i) = 1.0 / denom;
			alpha(i) = rho(i)*inner(s[i], q);
			q -= alpha(i) * y[i];
		}

		if (filter != nullptr) { filter(q); }
		else { q = gamma_k*q; }

		// L-BFGS second - loop recursion
		for(int i = 0; i < iter; ++i)
		{
			Scalar beta = rho(i)*inner(q, y[i]);
			q += (alpha(i) - beta)*s[i];
		}
	}

	//
	// Perform step
	//
	{
		// If our hess approx is bad and we start going in
		// the wrong direction, restart memory
		Scalar step_size = 1.0;
		Scalar dir = inner(q, grad);
		if (dir <= 0)
		{
			q = grad;
			max_iters -= k; // Restart memory
			k = 0;
			step_size = std::min(1.0, 1.0 / grad.template lpNorm<Eigen::Infinity>() );
		}

		// We've hit local minima, we have to exit
		if (q.squaredNorm() <= 0.0) {
			return obj_k;
		}

		descent = -q;
		x_last = x;
		obj_k = linesearch(x, grad, descent, step_size);
		if (num_iters >= options.min_iters && converged(obj_k,x_last,x,grad)) {
			return obj_k;
		}
	}

	//
	// Correction term
	//
	{
		s_temp = x - x_old;
		y_temp = grad - grad_old;

		// update the history
		if (k < options.M)
		{
			s[k] = s_temp;
			y[k] = y_temp;
		}
		else
		{
			for (int i=0; i<options.M - 1; ++i)
			{
				s[i] = s[i+1];
				y[i] = y[i+1];
			}
			s.back() = s_temp;
			y.back() = y_temp;
		}

		Scalar denom = inner(y_temp, y_temp);
		if (std::abs(denom) > 0.0)
		{
			gamma_k = inner(s_temp, y_temp) / denom;
		}
	}

	return obj_k;
}

template <typename MatrixType>
typename LBFGS<MatrixType>::Scalar
LBFGS<MatrixType>::bracketing_weakwolfe(
		MatrixType& x,
		MatrixType& g,
		const MatrixType &drt,
		Scalar &step) const
{
    using namespace Eigen;

	Scalar c1 = 10e-4;
	Scalar c2 = 0.9;
	if (step <= 0.0) { step = 1.0; }
    int x_cols = x.cols();

	g = MatrixType::Zero(x.rows(), x_cols);
	MatrixType xp = x;
	const Scalar fx_init = gradient(x,g);
	Scalar fx = fx_init;
	
	Scalar dg_init = inner(g, drt);
    if (dg_init > 0) {
        throw std::runtime_error("direction increases objective");
    }

	const Scalar test_decr = c1 * dg_init;
	Scalar lower = 0;
	Scalar upper = std::numeric_limits<Scalar>::infinity();
	int maxiter = 200;
	int iter = 0;
	for (iter = 0; iter < maxiter; iter++)
	{
		x = xp + step * drt;
		fx = gradient(x, g);
		if (fx > fx_init + step * test_decr) { // Armijo rule
			upper = step;
		}
		else
		{		
			Scalar dg = inner(g, drt);			
			if(dg < c2 * dg_init){ // Weak wolfe
				lower = step;
			} else {
				break; // both met
			}
		}
		step = std::isinf(upper) ? 2*step : lower/2 + upper/2;
	}
	return fx;

} // end linesearch

template <typename MatrixType>
bool LBFGS<MatrixType>::default_converged(
		Scalar curr_obj,
		const MatrixType& xprev,
		const MatrixType& x,
		const MatrixType& g) const
{
	(void)(curr_obj);
	(void)(xprev);
	Scalar gnorm = g.norm();
	if (gnorm <= options.abs_tol) {
	    return true;
	}
	Scalar xnorm = x.norm();
	if (gnorm < options.rel_tol * xnorm) {
	    return true;
	}
	return false;
}

template <typename MatrixType>
typename LBFGS<MatrixType>::Scalar
LBFGS<MatrixType>::inner(
    const MatrixType &a,
    const MatrixType &b) const
{
    int cols = std::min(a.cols(), b.cols());
    Scalar dot = 0;
    for (int i=0; i<cols; ++i) {
        dot += a.col(i).dot(b.col(i));
    }
    return dot;
} // end inner product

} // end ns mcl

#endif