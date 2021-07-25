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

#ifndef MCL_MINIMIZER_H
#define MCL_MINIMIZER_H

#include "Problem.hpp"
#include "Backtracking.hpp"
#include "MoreThuente.hpp"
#include "WolfeBisection.hpp"
#include <memory>

namespace mcl {
namespace optlib {

// The different line search methods currently implemented
enum class LSMethod {
	None = 0, // use step length = 1, not recommended
	MoreThuente, // TODO test this one for correctness
	Backtracking, // basic backtracking with sufficient decrease
	BacktrackingCubic, // backtracking with cubic interpolation
	WeakWolfeBisection // slow
};

// Trust region subproblem method (see TrustRegion.hpp)
enum class TRMethod {
	CauchyPoint,
	DogLeg
};

//
// Base class for optimization algs
//
template<typename Scalar, int DIM>
class Minimizer {
public:
	typedef Eigen::Matrix<Scalar,DIM,1> VecX;
	static const int FAILURE = -1; // returned by minimize if an error is encountered

	struct Settings {
		int verbose; // higher = more printouts
		int max_iters; // usually changed by derived constructors
		int ls_max_iters; // max line search iters
		Scalar ls_decrease; // sufficient decrease param
		LSMethod ls_method; // see LSMethod (above)
		TRMethod tr_method; // see TRMethod (above)

		Settings() : verbose(0), max_iters(100),
			ls_max_iters(100000), ls_decrease(1e-4),
			ls_method(LSMethod::BacktrackingCubic),
			tr_method(TRMethod::DogLeg)
			{}
	} m_settings;

	//
	// Performs optimization
	//
	virtual int minimize(Problem<Scalar,DIM> &problem, VecX &x) = 0;


protected:

	// Line search method/options can be changed through m_settings.
	Scalar linesearch(const VecX &x, const VecX &p, Problem<Scalar,DIM> &prob, double alpha0) const {
		double alpha = alpha0;
		int mi = m_settings.ls_max_iters;
		int v = m_settings.verbose;
		Scalar sd = m_settings.ls_decrease;
		switch( m_settings.ls_method ){
			default:{
				alpha = Backtracking<Scalar,DIM>::search(v, mi, sd, x, p, prob, alpha0);
			} break;
			case LSMethod::None: { alpha = 1.0; } break;
			case LSMethod::MoreThuente: {
				alpha = MoreThuente<Scalar,DIM>::search(x, p, prob, alpha0);
			} break;
			case LSMethod::Backtracking: {
				alpha = Backtracking<Scalar,DIM>::search(v, mi, sd, x, p, prob, alpha0);
			} break;
			case LSMethod::BacktrackingCubic: {
				alpha = BacktrackingCubic<Scalar,DIM>::search(v, mi, sd, x, p, prob, alpha0);
			} break;
			case LSMethod::WeakWolfeBisection: {
				alpha = WolfeBisection<Scalar,DIM>::search(v, mi, x, p, prob, alpha0);
			} break;
		}
		return alpha;
	} // end do linesearch

}; // class minimizer

} // ns optlib
} // ns mcl

#endif
