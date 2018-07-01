%EXdiffusion_rrgmes  Example script, inverse diffusion problem

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% Clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots').
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width.
MS = 10; % Size of markers on plots.

rng(0);  % Cake sure this test is repeatable.

% Define test problem
n = 64;                                   % Problem size.
NoiseLevel = 0.005;                       % Relative noise level in data.
[A,b,x,ProbInfo] = PRdiffusion(n);        % Get the test problem.
[bn,NoiseInfo] = PRnoise(b, NoiseLevel);  % Add Gaussian noise.

% Compute RRGMRES reconstruction, with these options:
%    x_true     - allows the method to compute error norms between iterates
%                 and true solution
%    NoStop     - set this to 'on' so that the method continues to max
%                 number of iterations, even if stopping rule is satisfied
%    NoiseLevel - needed to use discrepancy principle for stopping rule
%    eta        - safety factor in the discrepancy principle
K = 1:100;   % Iterations.
eta = 1.01;  % Safety factor.
options = IRset('x_true', x, 'NoStop', 'on', 'NoiseLevel', NoiseLevel, 'eta', eta);
% Now run RRGMRES.
[X, IterInfo] = IRrrgmres(A,bn,K,options);

% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends.
if strcmp(dispres, 'subplots')
    figure(1), clf
    subplot(3,3,1)
    PRshowx(x,ProbInfo)
    title('True solution','interpreter','latex','fontsize',16)
    set(gca,'fontsize',12)
    %
    subplot(3,3,4)
    PRshowb(bn,ProbInfo)
    title('Noisy data','interpreter','latex','fontsize',16)
    set(gca,'fontsize',12)
    %
    subplot(3,3,2)
    semilogy(K,IterInfo.Enrm,'linewidth',1.5)
    title('Error history','interpreter','latex','fontsize',16)
    axis([0 max(K) 0.15 100])
    set(gca,'fontsize',12)
    % hl = legend('{\tt info.Enrm}');
    % hl = legend('IRrrgmres error');
    % set(hl,'interpreter','latex')
    hold on
    semilogy(IterInfo.BestReg.It, IterInfo.BestReg.Enrm, 'ro', 'LineWidth', 1.5, 'MarkerSize', 6)
    semilogy(IterInfo.StopReg.It, IterInfo.StopReg.Enrm, 'ms', 'LineWidth', 1.5, 'MarkerSize', 6)
    hl = legend('IRrrgmres error','optimal stopping iteration', ...
      'DP stopping iteration','location','North');
    set(hl,'interpreter','latex')
    %
    subplot(3,3,5)
    PRshowx(IterInfo.BestReg.X,ProbInfo)
    title(['Best sol., $k$ = ',num2str(IterInfo.BestReg.It)],'interpreter','latex','fontsize',16)
    set(gca,'fontsize',12)
    %
    subplot(3,3,3)
    semilogy(K,IterInfo.Rnrm,'-',K,eta*NoiseLevel*ones(size(K)),'--','linewidth',1.5)
    % hl = legend('{\tt info.Rnrm}','{\tt eta*NoiseLevel}','location','northwest');
    hl = legend('IRrrgmres residual','{\tt eta*NoiseLevel}','location','northwest');
    set(hl,'interpreter','latex')
    ht = title('Residual history','interpreter','latex','fontsize',16);
    set(ht,'interpreter','latex','Position',[71,1.0422e-2,0]);
    axis([0 max(K) 0.0045 0.01])
    set(gca,'fontsize',12)
    %
    subplot(3,3,6)
    PRshowx(IterInfo.StopReg.X,ProbInfo)
    title(['DP sol., $k$ = ',num2str(IterInfo.StopReg.It)],'interpreter','latex','fontsize',16)
    set(gca,'fontsize',12)
elseif strcmp(dispres, 'manyplots')
    figure(1), clf
    PRshowx(x,ProbInfo)
    title('True solution','interpreter','latex','fontsize',18)
    set(gca,'fontsize',24)
    %
    figure(2), clf
    PRshowb(bn,ProbInfo)
    title('Noisy data','interpreter','latex','fontsize',18)
    set(gca,'fontsize',24)
    %
    figure(3), clf
    semilogy(K,IterInfo.Enrm,'linewidth',LW)
    hold on
    semilogy(IterInfo.BestReg.It, IterInfo.BestReg.Enrm, 'ro', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(IterInfo.StopReg.It, IterInfo.StopReg.Enrm, 'ms', 'LineWidth', LW, 'MarkerSize', MS)
    hl = legend('IRrrgmres errors', 'optimal IRrrgmres stopping iteration','IRrrgmres DP stopping iteration');
    set(hl,'interpreter','latex','fontsize',28)
    title('Error history','interpreter','latex','fontsize',18)
    axis([0 max(K) 0.15 100])
    set(gca,'fontsize',24)
    %
    figure(4), clf
    PRshowx(IterInfo.BestReg.X,ProbInfo)
    title(['Best sol., $k$ = ',num2str(IterInfo.BestReg.It)],'interpreter','latex','fontsize',18)
    set(gca,'fontsize',24)
    %
    figure(5), clf
    semilogy(K,IterInfo.Rnrm,'-',K,eta*NoiseLevel*ones(size(K)),'--','linewidth',LW)
    hl = legend('{\tt info.Rnrm}','{\tt eta*NoiseLevel}','location','northwest');
    % hl = legend('IRrrgmres residual','{\tt eta*NoiseLevel}','location','northwest');
    set(hl,'interpreter','latex','fontsize',28)
    % ht = title('Residual history','interpreter','latex','fontsize',18);
    % set(ht,'interpreter','latex','Position',[71,1.0422e-2,0]);
    axis([0 max(K) 0.0045 0.01])
    set(gca,'fontsize',24)
    %
    figure(6), clf
    PRshowx(IterInfo.StopReg.X,ProbInfo)
    title(['DP sol., $k$ = ',num2str(IterInfo.StopReg.It)],'interpreter','latex','fontsize',18)
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
    figure(1), print -dpng -r300 EXdiffusion
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -depsc -r300 EXdiffusion_a
    figure(2), print -depsc -r300 EXdiffusion_d
    figure(3), print -depsc -r300 EXdiffusion_b
    figure(4), print -depsc -r300 EXdiffusion_e
    figure(5), print -depsc -r300 EXdiffusion_c
    figure(6), print -depsc -r300 EXdiffusion_f
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
    figure(1), saveas('EXdiffusion.fig')
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXdiffusion_a.fig')
    saveas(figure(2), 'EXdiffusion_b.fig')
    saveas(figure(3), 'EXdiffusion_c.fig')
    saveas(figure(4), 'EXdiffusion_d.fig')
    saveas(figure(5), 'EXdiffusion_e.fig')
    saveas(figure(6), 'EXdiffusion_f.fig')
end
cd(oldcd)