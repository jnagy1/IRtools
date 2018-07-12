function L = LaplacianMatrix2D(n)
%LaplacianMatrix2D  Auxiliary function for IR Tools
%
% L = LaplacianMatrix12(n)
%
% Builds finite difference approximation matrix for the 2D Laplacian.  Here
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

n2 = sqrt(n);
if n2 ~= fix(n2)
    error('When using 2D Laplacian, need size(A,2) = n = N*N')
end
T = LaplacianMatrix1D(n2);
L = kron(T, speye(n2)) + kron(speye(n2), T);