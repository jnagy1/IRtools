function [psfMatData, p] = constructMatrix_new( PSF, center, boundary, imsize )
%
%       psfMatData = constructMatrix_new( PSF, center, boundary, imsize );
%
%  Construct psfMatrix data.
%
%  Given a PSF and the locatio of the corresponding point source,
%  this function sets up the data needed to do efficient matrix-vector
%  multiplication.
%
%  This is meant only for spatially invariant PSFs.
%
%  Input:
%         PSF  -  PSF image array
%      center  -  location of point source
%    boundary  -  desired boundary condition
%      imsize  -  size of restored image
%
%  Output:
%   psfMatData -  array containing the (complex) data needed to
%                 do efficient matrix-vector multiplications.
%
switch boundary
    case 'periodic'
        padSize = imsize - size(PSF);
        PSF = padarray(PSF, padSize, 'post');
        psfMatData = fft2(circshift(PSF, 1-center));
        p = [];
    otherwise
        [U, S, V] = svd(PSF);
        [m, n] = size(PSF);
        imsize_pad = 2.^nextpow2(2*imsize);
        minU = abs(min(min(U(:,1))));
        maxU = max(max(abs(U(:,1))));
       if minU == maxU
           U = -U;
           V = -V;
       end
       %
       %  Find number of Kronecker terms to represent the matrix.
       %
       Rtol = 10*eps;
       MaxTerms = sum(diag(S)/S(1,1) >= Rtol);
       
       %
       %  Get quantities needed to multiply by the Toeplitz pieces:
       %
       bidx_top1 = 1:n-center(2)+1;
       bidx_top2 = center(2):n;
       bidx_bot1 = imsize_pad(2)-center(2)+2:imsize_pad(2);
       bidx_bot2 = 1:center(2)-1;
       cidx_top1 = 1:m-center(1)+1;
       cidx_top2 = center(1):m;
       cidx_bot1 = imsize_pad(1)-center(1)+2:imsize_pad(1);
       cidx_bot2 = 1:center(1)-1;
       
              
       BT = V(:,1:MaxTerms)*sqrt(S(1:MaxTerms,1:MaxTerms));
       CT = U(:,1:MaxTerms)*sqrt(S(1:MaxTerms,1:MaxTerms));
       
       BT_hat = zeros(imsize_pad(2), MaxTerms);
       CT_hat = zeros(imsize_pad(1), MaxTerms);
       
       BT_hat(bidx_top1,:) = BT(bidx_top2,:);
       BT_hat(bidx_bot1,:) = BT(bidx_bot2,:);
       CT_hat(cidx_top1,:) = CT(cidx_top2,:);
       CT_hat(cidx_bot1,:) = CT(cidx_bot2,:);
       lambdaT = fft(CT_hat);
       deltaT = fft(BT_hat);
       
       %
       %  Get quantities needed to multiply by the Hankel pieces:
       %
       bidx_top1 = 1:n-center(2);
       bidx_top2 = center(2)+1:n;
       bidx_bot1 = imsize_pad(2)-center(2)+1:imsize_pad(2)-1;
       bidx_bot2 = 1:center(2)-1;
       cidx_top1 = 1:m-center(1);
       cidx_top2 = center(1)+1:m;
       cidx_bot1 = imsize_pad(1)-center(1)+1:imsize_pad(1)-1;
       cidx_bot2 = 1:center(1)-1;
              
       BH_hat = zeros(imsize_pad(2), MaxTerms);
       CH_hat = zeros(imsize_pad(1), MaxTerms);
       BH_hat(bidx_top1,:) = BT(bidx_top2,:);
       BH_hat(bidx_bot1,:) = BT(bidx_bot2,:);
       CH_hat(cidx_top1,:) = CT(cidx_top2,:);
       CH_hat(cidx_bot1,:) = CT(cidx_bot2,:);
       lambdaH = fft(CH_hat);
       deltaH = fft(BH_hat);
       
       %
       %  Need to permute these with a shift matrix:
       %
       pB = [1;(imsize_pad(2):-1:2)'];
       pC = [1;(imsize_pad(1):-1:2)'];
       
       lambdaH = lambdaH(pC,:);
       deltaH = deltaH(pB,:);
       
       psfMatData = zeros(imsize_pad(1), imsize_pad(2), 4);
       psfMatData(:,:,1) = lambdaT*deltaT';
       psfMatData(:,:,2) = lambdaT*deltaH';
       psfMatData(:,:,3) = lambdaH*deltaT';
       psfMatData(:,:,4) = lambdaH*deltaH';
       p = [pB, pC];
       
end
       
       

