function PRshowb(b,ProbInfo,h)
%PRshowb Show the right-hand side b (the data) in IR Tools
%
% PRshowb(b,ProbInfo)
% PRshowb(b,ProbInfo,h)
%
% This function uses the information stored in ProbInfo to display the
% right-hand side b in the correct way for that particular test problem.
%
% An optional third parameters specifies a figure handle.
%
% See also: PRshowx

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

if nargin < 3
  h = [];
end
if isempty(h)
  h = gcf;
end

switch ProbInfo.bType
    case {'image2D'}
        switch ProbInfo.problemType
            case {'deblurring'}
                figure(h)
                ImageDisplayGui2D(reshape(b, ProbInfo.bSize), h)
            case {'tomography'}
                figure(h)
                p = ProbInfo.bSize(1);
                ntheta = ProbInfo.bSize(2);
                ImageDisplayGui2D(reshape(b,[p,ntheta]), h)
        end
    case {'surf2D'}
        switch ProbInfo.problemType
            case{'nmr'}
                [tau2g, tau1g] = meshgrid(ProbInfo.tau2, ProbInfo.tau1);
                figure(h)
                surf(log10(tau2g), log10(tau1g), reshape(b, ProbInfo.bSize))
                shading interp
                axis tight
            case {'invinterp2'}
                figure(h)
                plot3(ProbInfo.Xr,ProbInfo.Yr,b,'.', 'color', [0 0 0.7])
        end
    case {'fem'}
        trisurf(ProbInfo.elmtab,ProbInfo.elmX,ProbInfo.elmY,b)
        axis([0 1 0 1 -0.1 1])
        shading interp
    otherwise
        error('Can only do image2D, surf2D, and fem')
end