function PRshowx(x, ProbInfo,h)
%PRshowx Show the solution x in IR Tools
%
% PRshowx(x,ProbInfo)
% PRshowx(x,ProbInfo,h)
%
% This function uses the information stored in ProbInfo to display the
% solution x in the correct way for that particular test problem.
%
% An optional third parameters specifies a figure handle.
%
% See also: PRshowb

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

switch ProbInfo.xType
    case {'image2D'}
        switch ProbInfo.problemType
            case {'deblurring','tomography'}
                figure(h)
                ImageDisplayGui2D(reshape(x, ProbInfo.xSize), h)
        end
    case {'surf2D'}
        switch ProbInfo.problemType
            case{'nmr'}
                [T2g, T1g] = meshgrid(ProbInfo.T2, ProbInfo.T1);
                figure(h)
                surf(log10(T2g), log10(T1g), reshape(x, ProbInfo.xSize))
                shading interp
                axis tight
            case {'invinterp2'}
                figure(h)
                n = ProbInfo.xSize(1);
                points = linspace(0, 1, n);
                [Xg,Yg] = meshgrid(points);
                plot3(Xg(:), Yg(:), x, '.', 'color', [0 0 0.7])
        end
    case {'fem'}
        trisurf(ProbInfo.elmtab,ProbInfo.elmX,ProbInfo.elmY,x)
        axis([0 1 0 1 -0.1 1])
        shading interp
    otherwise
        error('Can only do image2D, surf2D, and fem')
end