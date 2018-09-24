function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.8 $  $Date: 2002/03/15 15:58:11 $

out = repmat(feval(class, value), size);

