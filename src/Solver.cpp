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

using namespace mpm;

bool Solver::initialize(){

	// Fill a sphere with particles
	Eigen::Vector3d sphere_center(3.f,3.f,3.f);
	double sphere_rad = 1.0;
	double step = 0.2;
	double p_mass = 0.01;
	for( double x=sphere_center[0]-sphere_rad; x<sphere_center[0]+sphere_rad; x+=step ){
	for( double y=sphere_center[1]-sphere_rad; y<sphere_center[1]+sphere_rad; y+=step ){
	for( double z=sphere_center[2]-sphere_rad; z<sphere_center[2]+sphere_rad; z+=step ){
		Eigen::Vector3d pos( x, y, z );
		if( (pos-sphere_center).norm()<sphere_rad ){			
			Particle *p = new pNeoHookean();
			p->x = pos;
			p->m = p_mass;
			m_particles.push_back(p);
		}
	}}}

	MPM::cellsize = Eigen::Vector3d(.5,.5,.5);
	MPM::gridsize = Eigen::Vector3i(16,16,16);

	for( int i=0; i<MPM::gridsize[0]*MPM::gridsize[1]*MPM::gridsize[2]; ++i ){ m_grid.push_back( new GridNode(i) ); }

	std::cout << "Gridsize: " << MPM::gridsize[0] << ' ' << MPM::gridsize[1] << ' ' << MPM::gridsize[2] << std::endl;
	std::cout << "Cellsize: " << MPM::cellsize[0] << ' ' << MPM::cellsize[1] << ' ' << MPM::cellsize[2] << std::endl;
	std::cout << "Num Particles: " << m_particles.size() << std::endl;
	std::cout << "Num Nodes: " << m_grid.size() << std::endl;

	// Set up the solver
	cppoptlib::Options settings;
	settings.maxIter = 20;
	settings.gradTol = 1e-8;
	solver.settings_ = settings;

	return true;
}


Solver::~Solver(){
	for( int i=0; i<m_particles.size(); ++i ){ delete m_particles[i]; } m_particles.clear();
	for( int i=0; i<m_grid.size(); ++i ){ delete m_grid[i]; } m_grid.clear();
	active_grid.clear();
}


//
//	Run a system step
//
bool Solver::step( float screen_dt ){

	// Compute timestep
	MPM::timestep_s = compute_timestep( screen_dt );
	MPM::clear_grid( active_grid );

	// Particle to grid
	MPM::p2g_mass( m_particles, m_grid, active_grid ); // step 1
	if( elapsed_s <= 0.f ){ MPM::init_volumes( m_particles, m_grid ); }
	MPM::p2g_velocity( m_particles, active_grid ); // step 1

	// Steps 4, 5 in course notes
//	explicit_solve();
	implicit_solve();
	MPM::grid_collision( active_grid );

	// Grid to particle, and advection
	// Steps 6, 7, 8 in course notes
	MPM::g2p_deformation( m_particles, m_grid );
	MPM::g2p_velocity( m_particles, m_grid );
	MPM::particle_advection( m_particles );

	elapsed_s += MPM::timestep_s;

	// Done
	return true;

} // end step


void Solver::explicit_solve(){

	#pragma omp parallel for
	for( int i=0; i<m_particles.size(); ++i ){
		Particle *p = m_particles[i];
		Eigen::Matrix3d p_F = p->get_deform_grad();
		p->tempP = p->vol * p->get_piola_stress( p_F ) * p->Fe.transpose();
	} // end loop particles

	// Velocity update
	#pragma omp parallel for
	for( int i=0; i<active_grid.size(); ++i ){
		GridNode *g = active_grid[i];

		Eigen::Vector3d f(0,0,0);
		for( int j=0; j<g->particles.size(); ++j ){
			Particle *p = g->particles[j];
			Eigen::Matrix3d p_F = p->get_deform_grad();
			f -= p->tempP*g->dwip[j];
		}

		g->v += ( MPM::timestep_s*f )/g->m + MPM::timestep_s*gravity;

	} // end for all active grid

} // end compute forces


double Objective::value_gradient(const cppoptlib::Vector<double> &v, cppoptlib::Vector<double> &grad){
	using namespace Eigen;

	double tot_energy = 0.0;

	//
	//	Compute Particle Deformation Gradients for new grid velocities
	//
	#pragma omp parallel for reduction (+:tot_energy)
	for( int i=0; i<solver->m_particles.size(); ++i ){
		Particle *p = solver->m_particles[i];
		Eigen::Matrix3d p_F = p->get_deform_grad();

		Matrix3d vel_grad; vel_grad.fill(0.f);

		std::vector<int> grid_nodes;
		std::vector<double> wip;
		std::vector<Eigen::Vector3d> dwip;
		MPM::make_grid_loop( p, grid_nodes, wip, dwip );
		for( int j=0; j<grid_nodes.size(); ++j ){
			int idx = solver->m_grid[ grid_nodes[j] ]->active_idx;
			Eigen::Vector3d curr_v( v[idx*3+0], v[idx*3+1], v[idx*3+2] );
			vel_grad += ( MPM::timestep_s*curr_v ) * dwip[j].transpose();
		}

		Eigen::Matrix3d newF = ( Matrix3d::Identity() + vel_grad ) * p_F; // eq 193 course notes
		p->tempP = p->vol * p->get_piola_stress( newF ) * p_F.transpose();

		// Update the total energy
		double e = p->get_energy_density( newF )*p->vol;
		tot_energy = tot_energy + e;

	} // end loop particles


	//
	//	Compute energy gradient
	//
	#pragma omp parallel for reduction (+:tot_energy)
	for( int i=0; i<solver->active_grid.size(); ++i ){

		GridNode *node = solver->active_grid[i];
		Eigen::Vector3d curr_v( v[i*3+0], v[i*3+1], v[i*3+2] );

		Eigen::Vector3d momentum_grad = node->m * ( curr_v - node->v );
		Eigen::Vector3d energy_grad(0,0,0);

		for( int j=0; j<node->particles.size(); ++j ){
			energy_grad += node->particles[j]->tempP * node->dwip[j];
		}

		grad.segment( i*3, 3 ) = momentum_grad + energy_grad;
		tot_energy = tot_energy + 0.5 * node->m * ( curr_v - node->v ).squaredNorm();

	} // end loop grid

	return tot_energy;
}


void Solver::implicit_solve(){

	// Initial guess of grid velocities
	cppoptlib::Vector<double> v( active_grid.size()*3 );
	#pragma omp parallel for
	for( int i=0; i<active_grid.size(); ++i ){
		for( int j=0; j<3; ++j ){
			v[i*3+j] = active_grid[i]->v[j];
		}
	}

	// Minimize
	Objective obj( this );
	solver.minimize( obj, v );

	// Copy solver results back to grid
	#pragma omp parallel for
	for( int i=0; i<active_grid.size(); ++i ){
		for( int j=0; j<3; ++j ){
			active_grid[i]->v[j] = v[i*3+j];
		}
		active_grid[i]->v +=  MPM::timestep_s*gravity; // explicitly add gravity
	}

} // end implicit solve


double Solver::compute_timestep( double screen_dt ){

	double max_v = 0.0;
	double max_timestep = 0.01;

	#pragma omp parallel for reduction(max:max_v)
	for( int i=0; i<m_particles.size(); ++i ){
		double curr_v = m_particles[i]->v.squaredNorm();
		if( curr_v > max_v ){ max_v = curr_v; }
	} max_v = sqrtf(max_v);

	if( max_v > 1e-8 ){ 
		double min_cellsize = std::min( MPM::cellsize[0], std::min( MPM::cellsize[1], MPM::cellsize[2] ) );
		double ts_adjust = 0.04;
		double dt = ts_adjust * min_cellsize / sqrt(max_v);
		if( dt > max_timestep ){ return max_timestep; }
		if( dt > screen_dt ){ return screen_dt; }
		return dt;
	}

	return max_timestep;

} // end adaptive timestep


