function [w, PSF, c] = whichSVD(PSF, c)
%
%          [w, PSF, center] = whichSVD(PSF, center)
%
%    Determine whether to use DCT or separable approximation
%    to compute SVD.pwd
%
%  Input: PSF  -- array 
%           c  -- center of PSF
%
%  Output:  w  -- string: 'dct' or 'kron'
%         PSF  -- possibly a symmetric approximation of the PSF
%

% J. Nagy, 01/12/03

%
%  first check symmetric:
%
n = size(PSF);
P = padarray(PSF, max(n-2*c+1,0), 'pre');
P = padarray(P, max(2*c-n-1,0), 'post');
P_sym = (P + flipdim(P,1) + flipdim(P,2) + flipdim(flipdim(P,1),2))/4;
e_sym = norm(P_sym(:) - P(:));
if e_sym <= sqrt(eps)*norm(P(:))
  w = 'dct';
else
  %
  % now check seprable:
  %
  [U, S, V] = svd(PSF);
  P_sep = S(1,1)*U(:,1)*V(:,1)';
  e_sep = norm(P_sep(:) - PSF(:));
  if e_sep < e_sym
    w = 'kron';
  else
    w = 'dct';
    np = size(P_sym);
    k = max(np - n);
    if k > 0
      PSF = P_sym(k+1:end-k,k+1:end-k);
      c = (size(PSF) + 1)/2;
      PSF = padarray(PSF, n - size(PSF), 'post');
    else
      PSF = P_sym;
    end
  end
end


