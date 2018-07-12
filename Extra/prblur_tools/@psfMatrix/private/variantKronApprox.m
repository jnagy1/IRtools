function [A, B] = variantKronApprox(K, m, n, Rtol)
%  
%      [A, B] = variantKronApprox(K, m, n, Rtol);
%
%  Construct a (possibly spatially variant) Kronecker product
%  approximation of a blurring matrix.
%
%  Input:    K -  psfMatrix object
%            m -  size of matrix A (assumed square)
%            n -  size of matrix B (assumed square)
%         Rtol -  relative tolerance (<=1); determines # kron products
%                 in the approximation.  If si/s1 > Rtol, where si is
%                 the ith largest singular values of PSF, then
%                 Ai (x) Bi will be included in the approximation.
%                 IF ONLY ONE TERM IN THE SUM IS DESIRED, USE Rtol = 1.
%  Output:   K -  kronMatrix object, B (x) A, which approximates
%                 the original psfMatrix.
%

%
%  07/07/03, J. Nagy
%            (This is a preliminary version)

PSF = K.psf;
center = PSF.center;
PSF = PSF.image;

%  First determine number of terms based on the PSF in the middle
%  of the image.
%
kk = ceil(size(PSF)/2);
PP = PSF{kk(1),kk(2)};
s = svd(PP);

nKronTerms = sum(s/max(s(:)) >= Rtol);



imsize = [m, n, 1];
nregions = [size(PSF), 1];
rsize = floor(imsize ./ nregions);
[RIidx, RJidx, RKidx] = region_indices(nregions, rsize);
RIidx(RIidx == max(RIidx(:))) = imsize(1);
RJidx(RJidx == max(RJidx(:))) = imsize(2);

k = ceil(nregions/2);

%
%  We only consider 2D case here
%
if imsize(3) > 1
  error('can only handle 2D case')
end

%
%  First compute B (which acts on rows of image)
%
B = cell(nKronTerms,1);
for jj = 1:nKronTerms
  B{jj} = zeros(m);
end
for i = 1:nregions(1)
  [U, S, V] = svd(PSF{i,k(2)});
  %
  %  check to make sure first column looks like a Gaussian, and is
  %  not inverted
  %
  minU = abs(min(min(U(:,1))));
  maxU = max(max(abs(U(:,1))));
  if minU == maxU
    U = -U;
  end
  cc = center{i, k(2)};
  for jj = 1:nKronTerms
    b = sqrt(S(jj,jj))*U(:,jj);
    switch i
      case{1, nregions(1)}
        switch K.boundary
          case 'zero'
            Bt = build_toep(b, cc(1), imsize(1));
          case 'reflexive'
            Bt = build_toep(b, cc(1), imsize(1)) + buildHank(b, cc(1), imsize(1));
          otherwise
            error('BC should be zero or reflexive')
        end
    otherwise
      Bt = build_toep(b, cc(1), imsize(1));
    end
    B{jj}(RIidx(i,1):RIidx(i,2),:) = Bt(RIidx(i,1):RIidx(i,2),:);
  end
end

%
%  Now compute B (which acts on columns of image)
%
A = cell(nKronTerms,1);
for jj = 1:nKronTerms
  A{jj} = zeros(n);
end
for i = 1:nregions(2)
  [U, S, V] = svd(PSF{k(1),i});
  %
  %  check to make sure first column looks like a Gaussian, and is
  %  not inverted
  %
  minV = abs(min(min(V(:,1))));
  maxV = max(max(abs(V(:,1))));
  if minV == maxV
    V = -V;
  end
  cc = center{k(1),i};
  for jj = 1:nKronTerms
    a = sqrt(S(jj,jj))*V(:,jj);
    switch i
      case{1, nregions(2)}
        switch K.boundary
          case 'zero'
            At = build_toep(a, cc(2), imsize(2));
          case 'reflexive'
            At = build_toep(a, cc(2), imsize(2)) + buildHank(a, cc(2), imsize(2));
          otherwise
            error('BC should be zero or reflexive')
        end
      otherwise
        At = build_toep(a, cc(2), imsize(2));
    end
    A{jj}(RJidx(i,1):RJidx(i,2),:) = At(RJidx(i,1):RJidx(i,2),:);
  end
end
