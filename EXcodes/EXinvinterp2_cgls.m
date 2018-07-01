%EXinvinterp2_cgls  Example script, inverse interpolation problem

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% Clear workspace and command window.
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots').
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width.
MS = 10; % Size of markers on plots.

rng(0);  % Make sure this test is repeatable

% Define the test problem.
n = 32;
NoiseLevel = 0.05;
[A,b,x,ProbInfo] = PRinvinterp2(n);
bn = PRnoise(b, NoiseLevel);

% Compute a standard CGLS solution.
K = 1:200;
options = IRset('x_true', x, 'NoStop', 'on');
[X1, IterInfo1] = IRcgls(A, bn, K, options);

% Compute CGLS with a 'Laplacian2' priorconditioner, with
% zero boundary conditions everywhere.
options.RegMatrix = 'Laplacian2D';
[X2, IterInfo2] = IRcgls(A, bn, K, options);

% Compute CGLS with "our" priorconditioner that enforces the
% "correct" boundary conditions.
L1 = spdiags([ones(n,1),-2*ones(n,1),ones(n,1)],[-1,0,1],n,n);
L1(1,1:2) = [1,0];
L1(n,n-1:n) = [0,1];
L2 = L1;
L2(n,n-1:n) = [-1,1];
L = [ kron(speye(n),L2) ; kron(L1,speye(n)) ];
L = qr(L,0);
options.RegMatrix = L;
[X3, IterInfo3] = IRcgls(A, bn, K, options);

% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends.
if strcmp(dispres, 'subplots')
    figure(1), clf
    subplot(3,3,1)
    PRshowx(X1(:,end), ProbInfo)
    title(['Standard CGLS, $k$ = ',num2str(IterInfo1.StopReg.It)],...
    'interpreter','latex','fontsize',14)
    axis([0 1 0 1 -0.2 1.2])
    %
    subplot(3,3,2)
    PRshowx(X2(:,end), ProbInfo)
    title(['{\tt ''Laplacian2D''}, $k$ = ',num2str(IterInfo2.StopReg.It)],...
    'interpreter','latex','fontsize',14)
    axis([0 1 0 1 -0.2 1.2])
    %
    subplot(3,3,3)
    PRshowx(X3(:,end), ProbInfo)
    title(['Our $L$, $k$ = ',num2str(IterInfo3.StopReg.It)],'interpreter','latex','fontsize',14)
    axis([0 1 0 1 -0.2 1.2])
else
    figure(1), clf
    PRshowx(X1(:,end), ProbInfo)
    axis([0 1 0 1 -0.2 1.2])
    title(['Standard CGLS, $k$ = ',num2str(IterInfo1.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',32)
    %
    figure(2), clf
    PRshowx(X2(:,end), ProbInfo)
    axis([0 1 0 1 -0.2 1.2])
    title(['{\tt ''Laplacian2D''}, $k$ = ',num2str(IterInfo2.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',32)
    %
    figure(3), clf
    PRshowx(X3(:,end), ProbInfo)
    title(['Our $L$, $k$ = ',num2str(IterInfo3.StopReg.It)],'interpreter','latex','fontsize',18)
    axis([0 1 0 1 -0.2 1.2])
    set(gca,'fontsize',32)
end

return

% A number of instructions useful to save the displayed figures follow;
% the defualt is not to execute them. If you wish to save the displayed
% figures in the dedicated 'Results' folder, please comment the above
% return statement
oldcd = cd;
if strcmp(dispres, 'subplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -dpng -r300 EXinviterp2
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -depsc -r300 EXinviterp2_a
    figure(2), print -depsc -r300 EXinviterp2_b
    figure(3), print -depsc -r300 EXinviterp2_c
end
cd(oldcd)

% Uncomment the following return statement if you wish to save the
% displayed figures as MATLAB figures

% return

oldcd = cd;
if strcmp(dispres, 'subplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXinviterp2.fig')
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXinvinterp2_a.fig')
    saveas(figure(2), 'EXinvinterp2_b.fig')
    saveas(figure(3), 'EXinvinterp2_c.fig')
end
cd(oldcd)