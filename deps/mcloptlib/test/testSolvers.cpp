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

#include <iostream>
#include "TestProblem.hpp"
#include "MCL/LBFGS.hpp"
#include "MCL/NonLinearCG.hpp"
#include "MCL/Newton.hpp"
#include "MCL/TrustRegion.hpp"
#include <memory>

using namespace mcl::optlib;
typedef std::shared_ptr< Minimizer<double,2> > MinPtr2; // rb
typedef std::shared_ptr< Minimizer<double,Eigen::Dynamic> > MinPtrD; // linear

bool test_linear( std::vector<MinPtrD> &solvers, std::vector<std::string> &names ){

	std::cout << "\nTest linear:" << std::endl;
	typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
	bool success = true;

	// Using multiple dimensions in the linear case, which changes
	// how Newton's performs a solve.
	std::vector<int> test_dims = { 4, 16, 64 };

	int n_test_dims = test_dims.size();
	for( int i=0; i<n_test_dims; ++ i ){

		int dim = test_dims[i];

		// High/Low dimensions, linear, and Eigen::Dynamic.
		// The solvers should work, since the Vecs/Mats are resized at run time.
		DynProblem cp(dim);

		int n_solvers = solvers.size();
		for( int i=0; i<n_solvers; ++i ){
			bool curr_success = true;

			// Special case for Newtons: should converge in 1 iter because its 2nd order
			if( names[i] == "newton" ){ solvers[i]->m_settings.max_iters = 1; }
			else { solvers[i]->m_settings.max_iters = 100; }
			solvers[i]->m_settings.verbose = 1;
			VecX x = VecX::Zero(dim);
			solvers[i]->minimize( cp, x );
			for( int i=0; i<dim; ++i ){
				if( std::isnan(x[i]) || std::isinf(x[i]) ){
					std::cerr << "(" << names[i] << ") Bad values in x: " << x[i] << std::endl;
					curr_success = false;
				}
			}
			VecX r = cp.A*x - cp.b;
			double rn = r.norm(); // x should minimize |Ax-b|
			if( rn > 1e-3 ){
				std::cerr << "(" << names[i] << ") Failed to minimize: |Ax-b| = " << rn << std::endl;
				curr_success = false;
			}

			if( curr_success ){ std::cout << "(" << names[i] << ") Linear (" << dim << "): Success" << std::endl; }
			else{ success = false; }

		} // end loop solvers

	} // end loop test dims

	return success;

}


bool test_rb( std::vector<MinPtr2> &solvers, std::vector<std::string> &names ){

	std::cout << "\nTest Rosenbrock:" << std::endl;
	Rosenbrock rb; // Dim = 2, also tests finite gradient/hessian
	bool success = true;

	int n_solvers = solvers.size();
	for( int i=0; i<n_solvers; ++i ){
		bool curr_success = true;

		solvers[i]->m_settings.max_iters = 1000;
		solvers[i]->m_settings.verbose = 1;
		Eigen::Vector2d x = Eigen::Vector2d::Zero();
		solvers[i]->minimize( rb, x );

		for( int i=0; i<2; ++i ){
			if( std::isnan(x[i]) || std::isinf(x[i]) ){
				std::cerr << "(" << names[i] << ") Bad values in x: " << x[i] << std::endl;
				curr_success = false;
			}
		}

		double rn = (Eigen::Vector2d(1,1) - x).norm();
		if( rn > 1e-8 ){
			std::cerr << "(" << names[i] << ") Failed to minimize: Rosenbrock = " << rn << std::endl;
			curr_success = false;
		}

		if( curr_success ){ std::cout << "(" << names[i] << ") Rosenbrock: Success" << std::endl; }
		else{ success = false; }
	}

	return success;
}

// Test what happens when the energy is already minimized
bool test_zero( std::vector<MinPtrD> &solvers, std::vector<std::string> &names ){

	std::cout << "\nTest zero energy:" << std::endl;
	typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
	bool success = true;
	int dim = 2;

	// High/Low dimensions, linear, and Eigen::Dynamic.
	// The solvers should work, since the Vecs/Mats are resized at run time.
	DynProblem cp(dim);

	int n_solvers = solvers.size();
	for( int i=0; i<n_solvers; ++i ){
		bool curr_success = true;

		// Set max iterations to 1, should exit right away (at min)
		solvers[i]->m_settings.max_iters = 1;
		solvers[i]->m_settings.verbose = 1;
		VecX x = cp.A.inverse() * cp.b;
		solvers[i]->minimize( cp, x );

		for( int i=0; i<dim; ++i ){
			if( std::isnan(x[i]) || std::isinf(x[i]) ){
				std::cerr << "(" << names[i] << ") Bad values in x: " << x[i] << std::endl;
				curr_success = false;
			}
		}
		VecX r = cp.A*x - cp.b;
		double rn = r.norm(); // x should minimize |Ax-b|
		if( rn > 1e-10 ){
			std::cerr << "(" << names[i] << ") Failed to minimize: |Ax-b| = " << rn << std::endl;
			curr_success = false;
		}

		if( curr_success ){ std::cout << "(" << names[i] << ") Linear (" << dim << "): Success" << std::endl; }
		else{ success = false; }

	} // end loop solvers

	return success;
}


int main(int argc, char *argv[] ){
	srand(100);
	std::vector< std::string > names;
	std::vector< MinPtr2 > min2;
	std::vector< MinPtrD > minD;

	std::string mode = "all";
	if( argc == 2 ){ mode = std::string(argv[1]); }

	if( mode=="lbfgs" || mode=="all" ){
		min2.emplace_back( std::make_shared< LBFGS<double,2> >( LBFGS<double,2>() ) );
		minD.emplace_back( std::make_shared< LBFGS<double,Eigen::Dynamic> >( LBFGS<double,Eigen::Dynamic>() ) );
		names.emplace_back( "lbfgs" );
	}
	if( mode=="cg" || mode=="all" ){
		min2.emplace_back( std::make_shared< NonLinearCG<double,2> >( NonLinearCG<double,2>() ) );
		minD.emplace_back( std::make_shared< NonLinearCG<double,Eigen::Dynamic> >( NonLinearCG<double,Eigen::Dynamic>() ) );
		names.emplace_back( "cg" );
	}
	if( mode=="newton" || mode=="all" ){
		min2.emplace_back( std::make_shared< Newton<double,2> >( Newton<double,2>() ) );
		minD.emplace_back( std::make_shared< Newton<double,Eigen::Dynamic> >( Newton<double,Eigen::Dynamic>() ) );
		names.emplace_back( "newton" );
	}
	if( mode=="trustregion" || mode=="all" ){
		min2.emplace_back( std::make_shared< TrustRegion<double,2> >( TrustRegion<double,2>() ) );
		minD.emplace_back( std::make_shared< TrustRegion<double,Eigen::Dynamic> >( TrustRegion<double,Eigen::Dynamic>() ) );
		names.emplace_back( "trustregion" );
	}

	bool success = true;
	success &= test_linear( minD, names );
	success &= test_rb( min2, names );
	success &= test_zero( minD, names );
	if( success ){
		std::cout << "\nSUCCESS!" << std::endl;
		return EXIT_SUCCESS;
	}
	else{ std::cout << "\n**FAILURE!" << std::endl; }
	return EXIT_FAILURE;
}



