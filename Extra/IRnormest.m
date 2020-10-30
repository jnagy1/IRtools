function sigma = IRnormest(A, b)
%  We sometimes need an estimate of the 2-norm of the matrix A.  We want to 
%  allow for user-defined objects and function handles that implement matrix
%  vector multiplication with A.  So we cannot use MATLAB's built-in
%  normest, nor can we use svds.  So we'll use a few iterations of our 
%  Lanczos bidiagonalization as implemented in HyBR to get an estimate of 
%  the largest singular value.
% HyBRoptions = HyBRset;
% HyBRoptions = HyBRset(HyBRoptions, 'InSolv', 'Tikhonov', 'RegPar', 0, 'Iter', 5, 'Reorth', 'on', 'verbosity', 'off');
% [~, HyBRout] = HyBR_for_IRtools(A, b, [], HyBRoptions, 'off', 'off');
% sigma = max(svd(HyBRout.B));

optnrm.RegParam = 0;
optnrm.MaxIter = 5;
optnrm.Reorth = 'on';
optnrm.DecompOut = 'on';
optnrm.IterBar = 'off';
optnrm.verbosity = 'off';
[~, infonrm] = IRhybrid_lsqr(A, b, optnrm);
sigma = max(svd(infonrm.B));