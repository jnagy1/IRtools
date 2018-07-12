function Y = invariantMultiply_new( A, X)
%
%           Y = invariantMultiply( A, X );
%
%  This function computes the multiplication of a spatially invariant
%  point spread function (PSF) times an image 
%                y = A*x
%
%  Here we use a Kronecker product approach to better handle reflexive
%  boundary conditions
%
%  Input:
%            A  -  psfMatrix object
%            X  -  array containing the image to which the psfMatrix
%                  is to be multiplied.
%
%  Output:
%            Y  -  contains the result after PSF multiplication.
%

%  J. Nagy  7/5/2016

imsize = size( X );
psfMatData = A.matdata;

switch A.boundary
    case 'periodic'
        if A.transpose
            Y = real(ifft2(conj(psfMatData).*fft2(X)));
        else
            Y = real(ifft2(psfMatData.*fft2(X)));
        end
    case 'zero'
        p = A.p;
        padSize = [size(psfMatData,1),size(psfMatData,2)] - size(X);
        Xpad = padarray(X, padSize, 'post');
        Xtilde = fft2(Xpad);
        Xtilde = Xtilde(:,p(:,1));
        if A.transpose
            Y = conj(psfMatData(:,:,1)).*Xtilde;
        else
            Y = psfMatData(:,:,1).*Xtilde;
        end
        Y = real(ifft2(Y));
        Y = Y(:, p(:,1));
        Y = Y(1:imsize(1), 1:imsize(2));
    case 'reflexive'
        p = A.p;
        padSize = [size(psfMatData,1),size(psfMatData,2)] - size(X);
        Xpad = padarray(X, padSize, 'post');
        Xtilde = fft2(Xpad);
        Xtilde = Xtilde(:,p(:,1));
        if A.transpose
            term1 = conj(psfMatData(:,:,1)).*Xtilde;
            term2 = conj(psfMatData(:,:,2)).*Xtilde(:,p(:,1));
            term3 = conj(psfMatData(:,:,3)).*Xtilde(p(:,2),:);
            term4 = conj(psfMatData(:,:,4)).*Xtilde(p(:,2),p(:,1));
        else
            term1 = psfMatData(:,:,1).*Xtilde;
            term2 = psfMatData(:,:,2).*Xtilde;
            term2 = term2(:,p(:,1));
            term3 = psfMatData(:,:,3).*Xtilde;
            term3 = term3(p(:,2),:);
            term4 = psfMatData(:,:,4).*Xtilde;
            term4 = term4(p(:,2),p(:,1));
        end
        Y = real(ifft2(term1 + term2 + term3 + term4));
        Y = Y(:,p(:,1));
        Y = Y(1:imsize(1), 1:imsize(2));
    otherwise
        error('still working on this')
end



