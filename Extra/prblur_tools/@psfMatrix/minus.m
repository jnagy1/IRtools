function A = minus(arg1,arg2)
%
%   Overload Matrix subtraction operations for psfMatrix
%
%   Implements A-B where A and B are psfMatrices
%
%   Result returned is a psfMatrix
%

%  L. Perrone  4/16/01

A=arg1+(-arg2);

