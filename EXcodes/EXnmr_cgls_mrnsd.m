%EXnmr_cgls_mrnsd Example script, 2D NMR relaxometry

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% Clear workspace and command window.
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots').
dispres = 'subplots';
% dispres = 'manyplots';

LW = 2;  % Plot line width.
MS = 10; % Size of markers on plots.

rng(0);  % Make sure this test is repeatable.

% Define the test problem.
n = 64;
NoiseLevel = 0.05;
[A, b, x, ProbInfo] = PRnmr(n);
bn = PRnoise(b, NoiseLevel);

% Compute a CGLS solution, specifying the noise level for the discrepancy
% principle, and use true solution to compute error norms.
eta = 1.02;
options = IRset('x_true', x, 'NoiseLevel', NoiseLevel, 'eta', eta, 'NoStop', 'on');
[x_cgls, IterInfo_cgls] = IRcgls(A, bn, 1:500, options);

% Compute MRNSD reconstruction.
K = [1, 100:100:20000];
[x_mrnsd, IterInfo_mrnsd] = IRmrnsd(A, bn, K, options);

% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends.
if strcmp(dispres, 'subplots')
    figure(1), clf
    subplot(3,3,1), PRshowx(x, ProbInfo), colormap hsv
    title('True solution','interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
    %
    subplot(3,3,4), PRshowx(IterInfo_cgls.BestReg.X, ProbInfo), colormap hsv
    title('Best CGLS solution','interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
    %
    subplot(3,3,2), semilogy(IterInfo_mrnsd.Enrm,'linewidth',1.5)
    hold on
    hl = legend('{\tt info.Enrm}');
    % hl = legend('IRmrnsd error');
    set(hl,'interpreter','latex','fontsize',12)
    semilogy(IterInfo_mrnsd.BestReg.It, IterInfo_mrnsd.BestReg.Enrm, 'ro', 'LineWidth', 1.5, 'MarkerSize', 6)
    axis([0 max(K) 0.08 1.2])
    semilogy(IterInfo_mrnsd.StopReg.It, IterInfo_mrnsd.StopReg.Enrm, 'ms', 'LineWidth', 1.5, 'MarkerSize', 6)
    set(gca,'fontsize',12)
    title('Error history','interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
    %
    subplot(3,3,5), PRshowx(IterInfo_mrnsd.BestReg.X, ProbInfo)
    title(['Best MRNSD sol., $k$ = ',num2str(IterInfo_mrnsd.BestReg.It)],...
    'interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
    %
    subplot(3,3,3)
    semilogy(K,IterInfo_mrnsd.Rnrm(K),'-',K,eta*NoiseLevel*ones(size(K)),'--','linewidth',1.5)
    hl = legend('{\tt info.Rnrm}','{\tt eta*NoiseLevel}','location','northwest');
    % hl = legend('IRmrnsd residual','{\tt eta*NoiseLevel}','location','northwest');
    set(hl,'interpreter','latex')
    axis([0 max(K) 0.045 0.08])
    title('Residula history','interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
    %
    subplot(3,3,6), PRshowx(IterInfo_mrnsd.StopReg.X, ProbInfo)
    title(['DP MRNSD sol., $k$ = ',num2str(IterInfo_mrnsd.StopReg.It)],...
    'interpreter','latex','fontsize',24)
    set(gca,'fontsize',10)
elseif strcmp(dispres, 'manyplots')
    figure(1), clf
    PRshowx(x,ProbInfo), colormap hsv
    title('True solution','interpreter','latex','fontsize',24)
    set(gca,'fontsize',24)
    %
    figure(2), clf
    PRshowx(IterInfo_cgls.BestReg.X, ProbInfo), colormap hsv
    title('Best CGLS solution','interpreter','latex','fontsize',24)
    set(gca,'fontsize',24)
    %
    figure(3), clf
    semilogy(IterInfo_mrnsd.Enrm,'linewidth',LW)
    hold on
    semilogy(IterInfo_mrnsd.BestReg.It, IterInfo_mrnsd.BestReg.Enrm, 'ro', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(IterInfo_mrnsd.StopReg.It, IterInfo_mrnsd.StopReg.Enrm, 'ms', 'LineWidth', LW, 'MarkerSize', MS)
    axis([0 max(K) 0.08 1.2])
    % hl = legend('{\tt info.Enrm}','location','North');
    % set(hl,'interpreter','latex','fontsize',24)
    hl = legend('IRmrnsd error','optimal stopping iteration', ...
      'DP stopping iteration','location','North');
    set(hl,'interpreter','latex','fontsize',16)
    title('Error history','interpreter','latex','fontsize',24)
    set(gca,'fontsize',30)
    %
    figure(4), clf
    PRshowx(IterInfo_mrnsd.BestReg.X, ProbInfo), colormap hsv
    title(['Best MRNSD sol., $k$ = ',num2str(IterInfo_mrnsd.BestReg.It)],...
    'interpreter','latex','fontsize',24)
    set(gca,'fontsize',24)
    %
    figure(5), clf
    semilogy(K,IterInfo_mrnsd.Rnrm(K),'-',K,eta*NoiseLevel*ones(size(K)),'--','linewidth',LW)
    % hl = legend('{\tt info.Rnrm}','{\tt eta*NoiseLevel}','location','north');
    hl = legend('IRmrnsd residual','{\tt eta*NoiseLevel}','location','north');
    set(hl,'interpreter','latex','fontsize',24)
    title('Residual history','interpreter','latex','fontsize',24)
    axis([0 max(K) 0.045 0.08])
    set(gca,'fontsize',30)
    %
    figure(6), clf
    PRshowx(IterInfo_mrnsd.StopReg.X,ProbInfo), colormap hsv
    title(['DP MRNSD sol., $k$ = ',num2str(IterInfo_mrnsd.StopReg.It)],...
    'interpreter','latex','fontsize',24)
    set(gca,'fontsize',24)
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
    figure(1), print -dpng -r300 EXnmr
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    % save as eps
    figure(1), print -depsc -r300 EXnmr_a.eps
    figure(2), print -depsc -r300 EXnmr_d.eps
    figure(3), print -depsc -r300 EXnmr_b.eps
    figure(4), print -depsc -r300 EXnmr_e.eps
    figure(5), print -depsc -r300 EXnmr_c.eps
    figure(6), print -depsc -r300 EXnmr_f.eps
    % save as png
    figure(1), print -dpng -r300 EXnmr_a
    figure(2), print -dpng -r300 EXnmr_d
    figure(3), print -dpng -r300 EXnmr_b
    figure(4), print -dpng -r300 EXnmr_e
    figure(5), print -dpng -r300 EXnmr_c
    figure(6), print -dpng -r300 EXnmr_f
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
    figure(1), saveas('EXnmr.fig')
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXnmr_a.fig')
    saveas(figure(2), 'EXnmr_b.fig')
    saveas(figure(3), 'EXnmr_c.fig')
    saveas(figure(4), 'EXnmr_d.fig')
    saveas(figure(5), 'EXnmr_e.fig')
    saveas(figure(6), 'EXnmr_f.fig')
end
cd(oldcd)