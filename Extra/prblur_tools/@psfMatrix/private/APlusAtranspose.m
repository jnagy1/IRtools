
function S = APlusAtranspose( PSF, Acenter, boundary)
%
% This function exists specifically to handle the addition
% A + Atranspose, where A is a psfMatrix.  It's called by
% the overloaded method plus.m whenever the input matrices
% are transposes of each other.
%
%  INPUT: PSF = the psf.image field of the psfMatrix that's being added
%                to its own transpose
%         Acenter = the center coordinates of the input PSF(s)
%         boundary = the boundary condition of A
%
%  OUTPUT is a new psfMatrix with fields as follows:
%        .psf is created from the summation of the original psf fields.
%             Centers are aligned before adding.
%             Care is taken so that the resulting psf is no larger than
%             the original psf (if possible).
%        .matdata as created as usual (by constructMatrix)
%        .boundary is the same as A
%        .type is determined by number of images in the psf
%        .transpose = 0
%

%  L. Perrone  4/9/02

%disp('You are entering the transpose zone')
sa = size(PSF);

% Make everything look like a 3D image in order to accommodate 
% actual 3D images in the loop (i,j,k)
if length(sa)==1
  sa = [sa 1 1];
elseif length(sa) == 2
  sa = [sa 1];
end

newPsfs=cell(sa);
for k=1:sa(3);
 for j=1:sa(2);
  for i=1:sa(1);

    summand1 = PSF{i,j,k};
    summand1 = summand1/sum(summand1(:));
    summand2 = summand1;
    summand2 = flipdim(summand2,1);
    summand2 = flipdim(summand2,2);
    summand2 = flipdim(summand2,3);
    Bcenter{i,j,k} = size(summand1) - Acenter{i,j,k};

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

    % Add the psfs together to get the sum psf.
    newPsf = summand1 + summand2;
    newPsf = newPsf/sum(newPsf(:));

    % Now for the messy stuff.
    % Cut off excess rows &/or columns so the resultant PSF is 
    % no bigger than the input PSFs... but check the size of  
    % what you want to cut off before chopping.  If it's too large
    % a piece, leave it in.
    maxSize = max( Asize, Bsize);
    sumSize = size(newPsf);
    chopSize = sumSize - maxSize;
    if length(chopSize)==1
       chopSize=[chopSize 0 0];
    elseif length(chopSize)==2
       chopSize=[chopSize 0];
    end
    if (chopSize(1) > 0)
      newPsf = circshift(newPsf,[chopSize(1) 0 0]);
      rowNorms = zeros(2*chopSize(1),1);
      for ii=1:2*chopSize(1)
       rowNorms(ii) = norm(squeeze(newPsf(ii,:,:)));
      end
      blooble = zeros(chopSize(1)+1,1);
      for ii=1:chopSize(1)+1
        blooble(ii) = norm(rowNorms(ii:ii+chopSize(1)-1));
      end
      cutStart = find(min(blooble)==blooble);
      if length(cutStart) ~= 1 %% (if more than one min is identified)
         cutStart = cutStart(fix(length(cutStart)/2)+1);
      end
      if blooble(cutStart) < .01*norm(newPsf(:))
         newPsf = cat(1, newPsf(1:cutStart-1,:,:), newPsf(cutStart+chopSize(1):sumSize(1),:,:));
         newPsf = circshift(newPsf, [cutStart-1 0 0]);
     else 
        newPsf = circshift(newPsf,[-chopSize(1) 0 0]);
        disp('Row norms are too large for chopping')
      end
    end %%%   if (chopSize(1) > 0)

    if (chopSize(2) > 0)
      newPsf = circshift(newPsf,[0 chopSize(2) 0]);
      colNorms = zeros(2*chopSize(2));
      for ii=1:2*chopSize(2)
      colNorms(ii) = norm(squeeze(newPsf(:,ii,:)));
      end
      blooble2 = zeros(chopSize(2)+1,1);
      for ii=1:chopSize(2)+1
        blooble2(ii) = norm(colNorms(ii:ii+chopSize(2)-1));
      end
      cutStart2 = find(min(blooble2)==blooble2);
      if length(cutStart2) ~= 1 %% (if more than one min is identified)
         cutStart2 = cutStart2(fix(length(cutStart2)/2)+1);
      end
      if blooble2(cutStart2) < .01*norm(newPsf(:))
         newPsf = cat(2,newPsf(:,1:cutStart2-1,:), newPsf(:,cutStart2+chopSize(2):sumSize(2),:));
         newPsf = circshift(newPsf,[0 cutStart2-1 0]);
     else
        newPsf = circshift(newPsf, [0 chopSize(2) 0]);
        disp('Column norms are too large for chopping')
      end
    end %%%   if (chopSize(2) > 0)

    if (chopSize(3) > 0)
      newPsf = circshift(newPsf,[0 0 chopSize(3)]);
      stackNorms = zeros(2*chopSize(3));
      for ii=1:2*chopSize(3)
      stackNorms(ii) = norm(squeeze(newPsf(:,:,ii)));
      end
      blooble3 = zeros(chopSize(3)+1,1);
      for ii=1:chopSize(3)+1
        blooble3(ii) = norm(stackNorms(ii:ii+chopSize(3)-1));
      end
      cutStart3 = find(min(blooble3)==blooble3);
      if length(cutStart3) ~= 1 %% (if more than one min is identified)
         cutStart3 = cutStart3(fix(length(cutStart3)/2)+1);
      end
      if blooble3(cutStart3) < .01*norm(newPsf(:))
         newPsf = cat(3,newPsf(:,:,1:cutStart3-1), newPsf(:,:,cutStart3+chopSize(3):sumSize(3)));
         [ci,cj,ck]=find(max(newPsf(:))==newPsf);
         newPsf = circshift(newPsf,[0 0 cutStart3-1]);
      else
         newPsf = circshift(newPsf,[0 0 chopSize(3)]);
         disp('Stack norms are too large for chopping')
      end
    end %%%   if (chopSize(3) > 0)

    % Now that the chopping is complete, normalize the sum psf
    % and assign it to its compartment in newPsfs
    newPsfs{i,j,k} = newPsf/sum(newPsf(:));
  end
 end
end

% Now build the new psfMatrix from the summed psf image
NewPSF = psf(newPsfs);
S = psfMatrix(NewPSF,boundary);
S=S*2;


