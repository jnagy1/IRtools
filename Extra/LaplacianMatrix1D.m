function L = LaplacianMatrix1D(n)
%LaplacianMatrix1D  Auxiliary function for IR Tools
%
% L = LaplacianMatrix1D(n)
%
% Builds finite difference approximation matrix for the 1D Laplacian. Here
% we assume zero boundary conditions to get a square, nonsingular matrix L.
%
% Input:  n = number of grid points for x and y
% Output: L = matrix

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

e = ones(n,1);
L = spdiags([-e, 2*e, -e], [-1,0,1], n, n);
