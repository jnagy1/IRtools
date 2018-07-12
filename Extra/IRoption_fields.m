function allfields = IRoption_fields
%IRoption_fields Define option field names for IR methods 
%
% This function is used by IRset and IRget to create a cell array of
% all possible fields for the input options structure to IR methods.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

allfields = {'x0'; 
    'MaxIter'; 
    'MaxIterIn'; 
    'MaxIterOut'; 
    'x_true';
    'NoStop'; 
    'NoStopIn'; 
    'NoStopOut'; 
    'NE_Rtol'; 
    'IterBar'; 
    'BestEnrm'; 
    'RegMatrix'; 
    'RegParam'; 
    'Reorth'; 
    'restart'; 
    'TotIterMax'; 
    'ktotcount'; 
    'Ktot'; 
    'stopCrit'; 
    'stopGCV'; 
    'resflatTol'; 
    'sirt_method'; 
    'GCVweight'; 
    'GCVflatTol'; 
    'GCVminTol'; 
    'tolX'; 
    'stopOut'; 
    'stabOut';
    'thr0'; 
    'LSQRtols'; 
    'omega'; 
    'stoprule'; 
    'taudelta'; 
    'transf'; 
    'none';
    'enrichment'; 
    'stdCGLS_out'; 
    'noiseType'; 
    'noiseParam'; 
    'MatUpdate';
    'tollalpha'; 
    'verbosity'; 
    'nonneg'; 
    'Ubound'; 
    'relaxParam'; 
    'RR_Rtol';
    'NoiseLevel';
    'eta'; 
    'RegParam0';
    'xMin';
    'xMax';
    'xEnergy';
    'trunc';
    'inSolver'; 
    'adaptConstr';
    'nonnegativity';
    'DecompOut'};


