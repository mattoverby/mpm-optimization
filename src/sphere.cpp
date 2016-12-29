// Copyright (c) 2016 University of Minnesota
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
//
// By Matt Overby (http://www.mattoverby.net)

#include "Solver.hpp"
#include <vdb.h>

int main(int argc, char *argv[]){
	using namespace mpm;

	double timestep = 0.01;

	Solver solver;
	solver.initialize();

	std::cout << "Hold enter to step" << std::endl;; 
	int max_steps = 300;
	for( int i=0; i<max_steps; ++i ){
		std::cin.get(); 

		solver.step( timestep );
		vdb_frame();
		for( int i=0; i<solver.m_particles.size(); ++i ){
			Particle *p = solver.m_particles[i];
			vdb_point(p->x[0],p->x[1],p->x[2]);
		}
	}

	return 0;
}

