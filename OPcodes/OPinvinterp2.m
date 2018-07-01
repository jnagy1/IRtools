function w = OPinvinterp2(u,X,Y,Xr,Yr,method,tflag)
% OPinvinterp2 The forward computation and its adjoint for PRinvinterp2
%
% w = OPinvinterp2(u,X,Y,Xr,Yr,method,tflag)
%
% Performs either the forward computation, its ajoint, or returns the
% problem dimensions.
%
% Input: u      - the vector to be operated upon
%        X,Y    - coordinates to the mesh points
%        Xr,Yr  - coordinate sets for the data points
%        method - string that defines the interpolation method:
%                 'nearest', 'linear', 'cubic' or 'spline'
%        tflag  - string that determines the computation:
%                 'notransp', 'transp' or 'size'

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

N = length(X);

if strcmpi(tflag,'size')
    w(1) = length(Xr);
    w(2) = N^2;
elseif strcmpi(tflag,'notransp')
    u = reshape(u,N,N);
    w = interp2(X,Y,u,Xr,Yr,method);
    w = w(:);
else
    if strcmpi(method,'nearest')
        Xr = round((N-1)*Xr)+1;
        Yr = round((N-1)*Yr)+1;
        w = zeros(N,N);
        for k=1:N^2
            w(Yr(k),Xr(k)) = w(Yr(k),Xr(k)) + u(k);
        end
        w = w(:);
    elseif strcmpi(method,'linear')
        dx = rem((N-1)*Xr,1);
        dy = rem((N-1)*Yr,1);
        alpha = dx.*dy;
        beta  = dx.*(1-dy);
        gamma = (1-dx).*dy;
        delta = (1-dx).*(1-dy);
        Xr = floor((N-1)*Xr)+1;
        Yr = floor((N-1)*Yr)+1;
        w = zeros(N,N);
        for k=1:N^2
            I = sub2ind([N,N],[Yr(k) Yr(k)   Yr(k)+1 Yr(k)+1],...
                              [Xr(k) Xr(k)+1 Xr(k)   Xr(k)+1]    );
            w(I) = w(I) + [delta(k),beta(k),gamma(k),alpha(k)]*u(k);
        end
        w = w(:);
    elseif strcmpi(method,'cubic') || strcmpi(method,'spline')
        w = zeros(N^2,1);
        I = 0;
        for j=1:N
            for i=1:N
                I = I + 1;
                eI = zeros(N,N);
                eI(i,j) = 1;
                aI = interp2(X,Y,eI,Xr,Yr,method);
                w(I) = aI(:)'*u;
            end
        end
    else
       error('''transp'' not implemented for this interpolation method')
    end
end