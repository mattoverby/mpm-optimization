// Copyright (c) 2016 Matt Overby
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

#include "Solver.hpp"
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char* argv[])
{
    using namespace mpm;
    using namespace Eigen;

    Solver solver;
    solver.initialize();

    MatrixXd V;
    solver.get_vertices(V);

    // Also add ground plane
    MatrixXd gV(4, 3);
    gV.row(0) = RowVector3d(0, 0.75, 0);
    gV.row(1) = RowVector3d(6, 0.75, 0);
    gV.row(2) = RowVector3d(6, 0.75, 6);
    gV.row(3) = RowVector3d(0, 0.75, 6);
    MatrixXi gF(2, 3);
    gF.row(0) = RowVector3i(0, 2, 1);
    gF.row(1) = RowVector3i(0, 3, 2);

    std::cout << "Press A to toggle animation" << std::endl;

    // Create viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.core().is_animating = false;
    viewer.data().set_mesh(gV, gF);
    viewer.data().add_points(V, Eigen::RowVector3d(1, 0, 0));
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool
    {
        if(viewer.core().is_animating)
        {
            solver.step(0.04);
            solver.get_vertices(V);
            viewer.data().clear();
            viewer.data().set_mesh(gV, gF);
            viewer.data().add_points(V, Eigen::RowVector3d(1, 0, 0));
        }
        return false;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
