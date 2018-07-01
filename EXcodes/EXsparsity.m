%EXsparsity  Example script, deblurring without and with sparsity

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% clear workspace and window.
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots').
% dispres = 'subplots';
dispres = 'manyplots';

rng(0) % make sure this test is repeatable

% Test problem: Gaussian deblurring of a sparse 'dotk' image (starry night).
n = 256;                        % Image is n-by-n.
NoiseLevel = 0.10;              % Relative noise level in data.
PRoptions.trueImage = 'dotk';   % Sparse test image.
[A, b, x, ProbInfo] = PRblur(n,PRoptions);
bn = PRnoise(b,'gauss',NoiseLevel);

% Compute the CGLS reconstruction.
options = IRset('MaxIter', 80, 'x_true', x, 'NoStop', 'on');
[Xcgls, info_cgls] = IRcgls(A, bn, options);

% Compute ell1 sparse reconstructions with default (GCV) inner stopping rule.
[Xell1_GCV, info_ell1_GCV] = IRell1(A, bn, options);

% Compute ell1 sparse reconstruction with DP as inner stopping rule.
options = IRset(options, 'RegParam', 'discrep', 'NoiseLevel', NoiseLevel,...
                'eta', 1.1);
[Xell1_DP, info_ell1_DP] = IRell1(A, bn, options);

options.NoStopOut = 'on';
K = 80;
[Xirn_DP, info_irn_DP] = IRirn(A, bn, K, options);

% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends.
if strcmp(dispres, 'subplots')
    figure(1), clf
    subplot(2,3,1)
    imagesc(reshape(x,n,n)), axis image off
    title('Exact image','fontsize',16,'interpreter','latex')
    subplot(2,3,2)
    imagesc(reshape(b,n,n)), axis image off
    title('Blurred image','fontsize',16,'interpreter','latex')
    subplot(2,3,3)
    imagesc(reshape(max(0,info_cgls.BestReg.X),n,n)), axis image off
    title(['Best IRcgls after ',num2str(info_cgls.BestReg.It),' iterations'],'fontsize',16,'interpreter','latex')
    subplot(2,3,4)
    imagesc(reshape(max(0,info_ell1_GCV.BestReg.X),n,n)), axis image off
    title(['Best IRell1 w/ GCV after ',num2str(info_ell1_GCV.BestReg.It),...
          ' iterations'],'fontsize',16,'interpreter','latex')
    subplot(2,3,5)
    imagesc(reshape(max(0,info_ell1_DP.BestReg.X),n,n)), axis image off
    title(['Best IRell1 w/ DP ',num2str(info_ell1_DP.BestReg.It),...
          ' iterations'],'fontsize',16,'interpreter','latex')
    subplot(2,3,6)
    imagesc(reshape(max(0,info_irn_DP.BestReg.X),n,n)), axis image off
    title(['Best IRirn w/ DP after', num2str(info_irn_DP.BestReg.It),...
         ' iterations'],'fontsize',16,'interpreter','latex')
elseif strcmp(dispres, 'manyplots')
    figure(1), clf
    XX = reshape(x,n,n);
    imagesc(XX), axis image off
    title('Exact image','fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(XX(227:256,227:256)), axis image, caxis([0 max(XX(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
    %
    figure(2), clf
    BB = reshape((max(b,0)),n,n);
    imagesc(BB), axis image off
    title('Blurred image','fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(BB(227:256,227:256)), axis image, caxis([0 max(BB(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
    %
    figure(3), clf
    XCG = reshape((max(0,info_cgls.BestReg.X)),n,n);
    imagesc(XCG), axis image off
    title(['Best IRcgls after ',num2str(info_cgls.BestReg.It),' iterations'],'fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(XCG(227:256,227:256)), axis image, caxis([0 max(XCG(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
    %
    figure(4), clf
    I = reshape((max(0,info_ell1_DP.BestReg.X)),n,n);
    imagesc(I), axis image off
    title(['Best IRell1 w/ GCV after ',num2str(info_ell1_GCV.BestReg.It),...
          ' iterations'],'fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(I(227:256,227:256)), axis image, caxis([0 max(I(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
    %
    figure(5), clf
    II = reshape((max(0,info_ell1_DP.BestReg.X)),n,n);
    imagesc(II), axis image off
    title(['Best IRell1 w/ DP after ',num2str(info_ell1_DP.BestReg.It),...
          ' iterations'],'fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(II(227:256,227:256)), axis image, caxis([0 max(II(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
    %
    figure(6), clf
    III = reshape((max(0,info_irn_DP.BestReg.X)),n,n);
    imagesc(III), axis image off
    title(['Best IRirn w/ DP after ',num2str(info_irn_DP.BestReg.It),...
          ' iterations'],'fontsize',16,'interpreter','latex')
    axes('units','normalized','position',[0.17 0.1 1.15 0.25]);
    imagesc(III(227:256,227:256)), axis image, caxis([0 max(III(:))])
    set(gca,'Xtick',[],'Ytick',[],'LineWidth',2,'Xcolor','red','Ycolor','red')
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
    figure(1), print -dpng -r300 EXsparsity
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -depsc -r300 EXsparsity_a
    figure(2), print -depsc -r300 EXsparsity_b
    figure(3), print -depsc -r300 EXsparsity_c
    figure(4), print -depsc -r300 EXsparsity_d
    figure(5), print -depsc -r300 EXsparsity_e
    figure(6), print -depsc -r300 EXsparsity_f
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
    figure(1), saveas('EXsparsity.fig')
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXsparsity_a.fig')
    saveas(figure(2), 'EXsparsity_b.fig')
    saveas(figure(3), 'EXsparsity_c.fig')
    saveas(figure(4), 'EXsparsity_d.fig')
    saveas(figure(5), 'EXsparsity_e.fig')
    saveas(figure(6), 'EXsparsity_f.fig')
end
cd(oldcd)