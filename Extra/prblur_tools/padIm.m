function [Ppad, center_pad] = padIm(P, dim, center)
%PADPSF Pad an array P with zeros to make it bigger.
%
%      Ppad = padIm(PSF, dim);
%
%  Pad P with zeros to make it an m-by-n array. 
%
%  Input:
%        P  Array (matrix)
%      dim  Desired dimension of padded array.  
%             If dim is a scalar, then n = m.
%  Optional Input:
%    center integer array of two values indicating location of "center"
%           of P. This is used for PSFs
%
%  Output:
%        P  Padded m-by-n array.

% Reference: See Chapter 4, 
%            "Deblurring Images - Matrices, Spectra, and Filtering"
%            by P. C. Hansen, J. G. Nagy, and D. P. O'Leary,
%            SIAM, Philadelphia, 2006.

%
% Set default parameters.
%
switch nargin
    case 1
        error('Need desired dimension of padded array')
    case 3
        if length(center) == 1
            error('center should be a vector with two integers')
        else
            ci = center(1); cj = center(2);
        end
end
if length(dim) == 1
    m = dim; n = dim;
else
    m = dim(1); n = dim(2);
end
%
% Pad the with zeros.
%
[mp, np] = size(P);
Ppad = zeros(m, n);
putidx = fix((size(Ppad) - size(P))/2) + 1;
Ppad(putidx(1):putidx(1)+mp-1, putidx(2):putidx(2)+np-1) = P;

if nargin == 3
    center_pad = [ci+putidx(1)-1, cj+putidx(2)-1];
else
    center_pad = [];
end

