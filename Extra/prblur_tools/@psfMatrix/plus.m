function S = plus( arg1, arg2 )
%
% This function overloads the binary + operator for psfMatrix 
% objects A and B
%
% Input: A and B, the two psfMatrix objects to be added
%
%  OUTPUT is a new psfMatrix with fields as follows:
%        .psf is created from the summation of the original psf fields
%        .matdata as created as usual (by constructMatrix)
%        .boundary is the same as A
%        .type is determined by number of images in the psf
%        .transpose = 0
%

%  L. Perrone  4/9/02


% First check that the original psf .image fields have the 
% same number of images
AA = arg1.psf;
Apsf=AA.image;  % the .image field is a cell with one or more psfs in it
sa = size(Apsf);
Acenter = AA.center;
BB = arg2.psf;
Bpsf=BB.image;  % the .image field is a cell with one or more images in it
sb = size(Bpsf);
Bcenter = BB.center;
if (sa ~= sb)
  error('Input psfMatrix objects must have the same number of PSF images')
end

% Make everything look like a 3D image in order to accommodate 
% actual 3D images in the loop (i,j,k)
if length(sa)==1
  sa = [sa 1 1];
  sb = [sb 1 1];
elseif length(sa) == 2
  sa = [sa 1];
  sb = [sb 1];
end

% If the input matrices are transposes of each other,
% send them off to APlusAtranspose to finish the job.
% Otherwise, continue along.
same=zeros(sa);
for i=1:sa(1)
  for j=1:sa(2)
      if all(all(Apsf{i,j}==Bpsf{i,j})) & all(all(Acenter{i,j}==Bcenter{i,j}) )
         same(i,j) = 1;
      else 
         same(i,j) = 0;
      end
  end
end
if all(same) & ~arg1.transpose & arg2.transpose
     S = APlusAtranspose(Apsf, Acenter ,arg1.boundary);
elseif all(same) & arg1.transpose & ~arg2.transpose
     S = APlusAtranspose(Bpsf, Bcenter, arg2.boundary);
else % the entire rest of the code is the "else" statement.
       % In this case, either the two input psfMatrices have
       % different psfs, or the user is adding A+A or A' + A'.
%     disp('You did not enter the transpose zone')

% Next check that the dimensions of the images are the same
% (Both must be 2D or both must be 3D)
if length(Acenter{1}) ~= length(Bcenter{1})
  error('Input psfMatrix objects must have the same dimensions')
end


% Now prepare to add.
newImages=cell(sa);
newCenters=cell(sa);

for k=1:sa(3);
 for j=1:sa(2);
  for i=1:sa(1);

    % If a psfMatrix is transposed, handle it here.
    summand1 = Apsf{i,j,k};
    if (arg1.transpose)
      summand1 = flipdim(summand1,1);
      summand1 = flipdim(summand1,2);
      summand1 = flipdim(summand1,3);
      Acenter{i,j,k} = size(summand1) +1 - Acenter{i,j,k};
    end
    summand2 = Bpsf{i,j,k};
    if(arg2.transpose)
      summand2 = flipdim(summand2,1);
      summand2 = flipdim(summand2,2);
      summand2 = flipdim(summand2,3);
      Bcenter{i,j,k} = size(summand2) +1 - Bcenter{i,j,k};
    end

    % Pad around each image, however necessary,
    % so that the resulting images are equal in size and have
    % the same center coordinates.
    Asize = size(summand1);
    centerA = Acenter{i,j,k};
    Bsize = size(summand2);
    centerB = Bcenter{i,j,k};
    Aprepad = centerB - centerA;
    Aprepadsize = (Aprepad + abs(Aprepad))/2;
    summand1 = padarray(summand1,Aprepadsize,'pre');
    newCenter = centerA + Aprepadsize;
    Bprepad = -Aprepad;
    Bprepadsize = (Bprepad + abs(Bprepad))/2;
    summand2 = padarray(summand2,Bprepadsize,'pre');
    Apost = Asize - centerA;
    Bpost = Bsize - centerB;
    Apostpad = Bpost - Apost;
    Apostpadsize = (Apostpad + abs(Apostpad))/2;
    summand1 = padarray(summand1,Apostpadsize,'post');
    Bpostpad = -Apostpad;
    Bpostpadsize = (Bpostpad + abs(Bpostpad))/2;
    summand2 = padarray(summand2,Bpostpadsize,'post');

    % Normalize each psf image so the sum of its entries is 1,
    % then add the psfs together to get the sum psf.
    % Assign this sum psf to its proper place in the image cell.
    summand1 = summand1/(sum(summand1(:)));
    summand2 = summand2/(sum(summand2(:)));
    newPsf = summand1 + summand2; 
    newImages{i,j,k} = newPsf/sum(newPsf(:));
    newCenters{i,j,k} = newCenter;
  end
 end
end

% Now build the new psfMatrix from the summed psf image
NewPsf = psf(newImages,newCenters);
S = psfMatrix(NewPsf,arg1.boundary);
S = S*2;

end % end the big "else" statement (the one which split A+A' methods
    % from A+B methods)




