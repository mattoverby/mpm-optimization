# mpm-optimiziation

Optimization-based implicit material point method with APIC,
implemented from the [SIGGRAPH course notes](http://web.cs.ucla.edu/~cffjiang/mpmcourse/mpm_course_nodes.pdf), [optimization paper](https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf), and [APIC paper](https://disney-animation.s3.amazonaws.com/uploads/production/publication_asset/104/asset/apic-aselle-final.pdf).

By Matt Overby  
[http://www.mattoverby.net](http://www.mattoverby.net)

# notes

- I coded up this project very quickly to test out some research ideas. It is not optimized or bug free.  
- Currently only Neohookean particles are implemented.
- There is no fancy renderer. I visualize the particles with [VDB](https://github.com/zdevito/vdb).

# dependencies

Included is a modified version of [CppOptimizationLibrary](https://github.com/PatWie/CppNumericalSolvers).  
Not included, but required for compilation, is [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page).

# license

Copyright (c) 2016 University of Minnesota

MCLSCENE Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)  
Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:  
1. Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.  
2. Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.  
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
