%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR Tools
% April 2018
% 
% This file is part of the IR Tools package and is distributed under the
% 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2018 Silvia Gazzola, University of Bath, Per Christian Hansen,
% Technical University of Denmark and James G. Nagy, Emory University.
% 
% Contact: s.gazzola@bath.ac.uk, pcha@dtu.dk, jnagy@emory.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Installating
%   IRtools_setup       - Set up search paths to IRtools
%
% Example scripts:
%   EXblur_cgls_hybrid  - Example script, speckle deblurring problem
%   EXdiffusion_rrgmres - Example script, inverse diffusion problem
%   EXinvinterp2_cgls   - Example script, inverse interpolation problem
%   EXnmr_cgls_mrnsd    - Example script, 2D NMR relaxometry
%   EXsparsity          - Example script, deblurring without and with sparsity
%
% Iterative methods:
%   IRart              - Algebraic Reconstruction Technique
%   IRcgls             - Conjugate Gradient algorithm for Least Squares problems
%   IRconstr_ls        - Least squares solver with box and energy constraints
%   IRell1             - Least squares solver with 1-norm penalization term
%   IRenrich           - Enriched CGLS method
%   IRfista            - FISTA algorithm for constrained least squares problems
%   IRget              - Get options for IR Tools functions
%   IRhtv              - Least squares solver with heuristic total variation penalization
%   IRhybrid_fgmres    - Hybrid version of FGMRES for enforcing 1-norm penalization
%   IRhybrid_gmres     - Hybrid version of GMRES algorithm for square systems
%   IRhybrid_lsqr      - Hybrid version of LSQR algorithm
%   IRirn              - Least squares solver with 1-norm penalization term
%   IRmrnsd            - Modified Residual Norm Steepest Descent method
%   IRnnfcgls          - Modified flexible CGLS for nonnegatively constrained LS problems
%   IRrestart          - Restarted Krylov subspace methods
%   IRrrgmres          - Range Restricted GMRES for square systems
%   IRset              - Set options for IR Tools functions
%   IRsirt             - Simuletaneous Iterative Reconstruction Technique
%
% Operators (functions) for some test problems:
%   OPdiffusion        - The forward computation and its adjoint for PRdiffusion
%   OPinvinterp2       - The forward computation and its adjoint for PRinvinterp2
%   OPnmr              - The forward computation and its adjoint for PRnmr
%
% Test problems and associated functions:
%   PRblur             - Generates data for use in image deblurring problems
%   PRblurdefocus      - Image deblurring problem with a defocus point spread function
%   PRblurgauss        - Image deblurring problem with a Gaussian point spread function
%   PRblurmotion       - Image deblurring problem with linear motion blur
%   PRblurrotation     - Image deblurring problem with rotational motion blur
%   PRblurshake        - Image deblurring problem with random camera motion (shaking)
%   PRblurspeckle      - Image deblurring problem with a speckle point spread function
%   PRdiffusion        - Generates data for use in an inverse diffusion problems
%   PRget              - Get options for IR Tools test problems
%   PRinvinterp2       - Generates data for use in 2D inverse interpolation problems
%   PRnmr              - Generates data for 2D NMR relaxometry problems
%   PRnoise            - Add noise of a specified distribution and level
%   PRseismic          - Generates data for seismic travel-time tomography problems
%   PRset              - Set options for IR Tools test problems
%   PRshowb            - Show the right-hand side b (the data) in IR Tools
%   PRshowx            - Show the solution x in IR Tools
%   PRspherical        - Generates data for spherical means tomography problems
%   PRtomo             - Generates data for X-ray tomographic reconstruction problems