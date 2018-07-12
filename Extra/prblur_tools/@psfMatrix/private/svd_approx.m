function [U, s, V] = svd_approx(K)
%
%   [U, s, V] = svd_approx(K);
%
%  Given a kronMatrix object, this function computes
%  an "SVD" factorization that is to be used for a
%  psfMatrix object.
%
%  On Entry:
%     K  -  kronMatrix obect
%
%  On Exit:
%     U, V  -  depends -- could be transform matrix objects (fft, dct)
%              or kronMatrix object.
%        s  -  column vector containing the singular values
%              of A.  Note that these are not sorted from largest
%              smallest.
%

%
%  first get Ua, Ub, Va, Vb
%
A1 = K.a{1};
B1 = K.b{1};
[Ua, Sa, Va] = svd(A1);
[Ub, Sb, Vb] = svd(B1);

%
%  Now find the best diagonal for Kron, and a measure of
%  the off diagonal elaments.
%
l = length(K);

%disp(sprintf('     Using %2d terms in the Kronecker product approximation,', l))
s_kron = kron(diag(Sa), diag(Sb));
OffDiagKron = 0;
for i = 2:l
  Da = Ua' * K.a{i} * Va;
  Db = Ub' * K.b{i} * Vb;
  
  sa = diag(Da);
  sb = diag(Db);
  
  OffDiagKron = OffDiagKron + norm(Da-diag(sa),'fro') + norm(Db-diag(sb),'fro');
  
  s_kron = s_kron + kron(sa, sb);
end

s_fft = zeros(size(s_kron));
OffDiagFFT = 0;
for i = 1:l
  n = length(K.a{i});
  scale = 1 / sqrt(n);

  Da = scale*fft( scale*fft(K.a{i})' )';
  Db = scale*fft( scale*fft(K.b{i})' )';
  
  sa = diag(Da);
  sb = diag(Db);
  
  OffDiagFFT = OffDiagFFT + norm(Da-diag(sa),'fro') + norm(Db-diag(sb),'fro');
  
  s_fft = s_fft + kron(sa, sb);
end

s_dct = zeros(size(s_kron));
OffDiagDCT = 0;
for i = 1:l
  
  Da = dct( dct(K.a{i})' )';
  Db = dct( dct(K.b{i})' )';
  
  sa = diag(Da);
  sb = diag(Db);
  
  OffDiagDCT = OffDiagDCT + norm(Da-diag(sa),'fro') + norm(Db-diag(sb),'fro');
  
  s_dct = s_dct + kron(sa, sb);
end

[mx, idx] = min([OffDiagKron, OffDiagFFT, OffDiagKron]);

switch idx
  case 1
    U = kronMatrix(Ua, Ub);
    V = kronMatrix(Va, Vb);
    s = s_kron;
    %disp('     with the SVD basis.')
  case 2
    U = transformMatrix('fft');
    U = U';
    V = transformMatrix('fft');
    V = V';
    s = s_fft;
    %disp('     with the FFT basis.')
  case 3
    U = transformMatrix('dct');
    U = U';
    V = transformMatrix('dct');
    V = V';
    s = s_dct;
    %disp('     with the DCT basis.')
  otherwise
    error('something is wrong here!')
end
