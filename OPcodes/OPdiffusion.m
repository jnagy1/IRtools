function w = OPdiffusion(u,R,S,Tsteps,tflag)
% OPdiffusion The forward computation and its adjoint for PRdiffusion
%
% w = OPdiffusion(u,R,S,Tfinal,Tsteps,tflag)
%
% Performs the forward computation or returns the problem dimensions.
%
% Input: u      - the vector to be operated upon
%        R,S    - two matrices needed for the iterations
%        Tfinal - diffusion time
%        Tsteps - number of time steps
%        tflag  - string that determines the computation:
%                 'notransp', 'transp' or 'size'

% Based on original code by Allan P. Engsig-Karup, DTU Compute.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

n = size(R,1);

if strcmpi(tflag,'size')
    w(1) = n;
    w(2) = n;
elseif strcmpi(tflag,'notransp')
    w = u;
    for i = 1:Tsteps
        w = R\(S*w);
    end
else
    w = u;
    for i = 1:Tsteps
        w = S'*((R')\w);
    end
end