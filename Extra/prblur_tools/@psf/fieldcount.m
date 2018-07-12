function num_fields = fieldcount(psf_obj)
%
%  Determines the number of fields in a psf object.
%  Used by psf child class methods.
%

% J. Nagy 2/11/01

num_fields = length(fieldnames(psf_obj));