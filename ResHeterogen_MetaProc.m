%%
Calculate = 1;
%Metamodel = 'Grids';
Metamodel = 'Single';
switch Metamodel
    case 'Grids'
    %Resource supply cloud moves along the diagonal (1:1 resource ratio)
    %SimulatIDs = {'Sigma=Mu^1.5', 'Sigma=Mu^2',    'Sigma=Mu^2.5',      'Sigma=Mu^3'};
    %resource supply cloud covers a grid
    SimulatIDs = {'GridSigma=Mu^1.5', 'GridSigma=Mu^2',    'GridSigma=Mu^2.5',      'GridSigma=Mu^3'};
    %SimulatIDs = {'Sigma=Mu^1.5_dyn05', 'Sigma=Mu^2_dyn05', 'Sigma=Mu^2.5_dyn05'};
    %SimulatIDs = {'Sigma=Mu^1.5_dyn50', 'Sigma=Mu^2_dyn50', 'Sigma=Mu^2.5_dyn50', 'Sigma=Mu^3_dyn50'};
    alphas = [1.5, 2, 2.5, 3];  %% sigma^2 = mu^alpha
    
    tResultsTbls = {};
    tParamsTbls = {};
    tFuncIndTbls = {};
    tRes_BEFTbls = {};
    tRstTbls = {};
    
    for i = 1:length(SimulatIDs)
        %try
            i
            [tResultsTbls{i}, tParamsTbls{i}, tFuncIndTbls{i}, tRes_BEFTbls{i}] = NSpecCompAdv_Cluster_Run_Model('ResHeterogen', SimulatIDs{i}, Calculate);
        %catch ME
            display(ME)
        %end
    end
    case 'Single'
           [tResultsTbls{i}, tParamsTbls{i}, tFuncIndTbls{i}, tRes_BEFTbls{i}] = NSpecCompAdv_Cluster_Run_Model('ResHeterogen', 'Single', Calculate);
         
end
%% plot results
close all;
fg1 = f_MakeFigure(1, [10, 10, 800, 800]);
fg11 = f_MakeFigure(11, [100, 10, 800, 800]);
fg2 = f_MakeFigure(2, [10, 10, 700, 800]);
fg3 = f_MakeFigure(3, [10, 10, 1200, 800]);  %%histrgram of N:P 
leg = {};
leg_alpha = {};
leg_muS = {};
clrs = linspecer(length(SimulatIDs));
LineHandles = [];
ResRatioRange_log = linspace(-3, 3, 30);
for it = 1:length(SimulatIDs)
    tResults = tResultsTbls{it};
    ResRatio_Avbl = [];
    ResRatio_Fin  = [];
    tFuncInd = tFuncIndTbls{it};
    ResSuppMean_1 = [];
    ResSuppVar_1 = [];
    BiomMean_1 = [];
    Resavailable_mean = [];
    Resavailable_var = [];
    for i = 1:size(tResults, 1)
        S = tResults.S{i};
        S = S(:, :, 1);
        ResSuppMean_1(end + 1) = mean(S(:));
        ResSuppVar_1(end + 1) = var(S(:));
        Nfin = tResults.Nfin{i};
        BiomMean_1(end + 1) = nanmean(Nfin(:));
        Rfin = tResults.Rfin{i};
        Resavailable_mean(end + 1) = nanmean(Rfin(:));
        Resavailable_var(end + 1) = var(Rfin(:));
    end
    indNonNaNBiom = ~isnan(BiomMean_1) & BiomMean_1 < 1000;
    
    %     subplot(1, 2, 1)
    %     loglog(ResSuppMean_1, ResSuppVar_1, 'o', 'Color', clrs(it, :));
    %     hold on
    %     f_Lbls('Mean resource', 'Variance')
    
    %%plot samples of resource distribution
    figure(fg1)
    subplot(2, 2, 4);
    loglog(ResSuppMean_1(indNonNaNBiom), ResSuppVar_1(indNonNaNBiom), '.', 'Color', clrs(it, :))
    hold on
    leg_alpha{end + 1} = sprintf('\\sigma^2(S)~ \\mu(S)^{%0.1f}', alphas(it));
    f_Lbls('Mean(resource suppl)', 'Var(res supplied)')
    if it == 4
        subplot(2, 2, it-1);
    else
        subplot(2, 2, it);
    end
    tResultsNotNan = tResults(indNonNaNBiom, :);
    PlotInd = unique(round(linspace(1, size(tResultsNotNan, 1), 5)));
    clrs2 = linspecer(length(PlotInd));
    for i = 1:length(PlotInd);
        S = tResultsNotNan.S{PlotInd(i)};
        S1 = S(:, :, 1);
        S2 = S(:, :, 2);

        Sfin = tResultsNotNan.Rfin{PlotInd(i)};
        Sfin1 = Sfin(:, :, 1);
        Sfin2 = Sfin(:, :, 2);
        
        figure(fg1);
        loglog(S1(:), S2(:), '.', 'Color', clrs2(i, :));
        hold on
        
        figure(fg3);
        subplot(2, length(SimulatIDs), it);    
        Counts = histcounts(log10(S1(:)./S2(:)), ResRatioRange_log);
        Counts = Counts./diff(ResRatioRange_log);
        Counts = Counts./(sum(Counts));
        plot(0.5 * (ResRatioRange_log(1:end-1) +ResRatioRange_log(2:end))', Counts, 'Color', clrs2(i, :), 'LineStyle', '-', 'Linewidth', 2)
        hold on

        subplot(2, length(SimulatIDs), it + length(SimulatIDs));    
        Counts_fin = histcounts(log10(Sfin1(:)./Sfin2(:)), ResRatioRange_log);
        Counts_fin = Counts_fin./diff(ResRatioRange_log);
        Counts_fin = Counts_fin./sum(Counts_fin);
        plot(0.5 * (ResRatioRange_log(1:end-1) +ResRatioRange_log(2:end))', Counts_fin, 'Color', clrs2(i, :), 'LineStyle', '-', 'Linewidth', 2)
        hold on
        if it == 1
            leg_muS{end + 1} = sprintf('\\mu(S)=%0.1e', mean(S1(:))); 
        end
    end
    figure(fg1);
    f_Lbls('Resource_1, supplied', 'Resource_2, supplied')
    
    figure(fg3);
    subplot(2, length(SimulatIDs), it);
    title(sprintf('\\sigma^2(S)~ \\mu(S)^{%0.1f}', alphas(it)))
    if it == 1
        f_Lbls('log_{10}(R1:R2), supplied', 'Count')
    else
        f_Lbls('log_{10}(R1:R2), supplied', '')
    end
    xlim([ResRatioRange_log(1), ResRatioRange_log(end)]);
    ylim([0, 1]);
    subplot(2, length(SimulatIDs), it + length(SimulatIDs));
    if it == 1
        f_Lbls('log_{10}(R1:R2), available', 'Count')
    else
        f_Lbls('log_{10}(R1:R2), available', '')
    end
    xlim([ResRatioRange_log(1), ResRatioRange_log(end)]);
    ylim([0, 0.8]);
    
    
    %%plot samples of resource distribution
    figure(fg11)
    subplot(2, 2, it);
    PlotInd = unique(round(linspace(1, size(tResults, 1), 5)));
    clrs2 = linspecer(length(PlotInd));
    for i = 1:length(PlotInd);
        Rfin = tResults.Rfin{PlotInd(i)};
        S1 = S(:, :, 1);
        S2 = S(:, :, 2);
        loglog(S1(:), S2(:), '.', 'Color', clrs2(i, :));
        hold on
    end
    f_Lbls('Resource_1, available', 'Resource_2, available')
    
    
    figure(fg2);
    ax1 = subplot(3, 1, 1);
    tResults.Simpson_ESN_glob = tFuncInd.Simpson_ESN_glob;
    tResults.RessuppliedGrp =  10.^(round(log10(ResSuppMean_1')*10)/10);
    tResults.BiomMeanGrp =  10.^(round(log10(BiomMean_1')*10)/10);

    tDivonRes = grpstats(tResults(:, {'Simpson_ESN_glob', 'RessuppliedGrp'}), {'RessuppliedGrp'}, {'mean', 'std'});
    hndls = shadedErrorBar(tDivonRes.RessuppliedGrp, tDivonRes.mean_Simpson_ESN_glob, tDivonRes.std_Simpson_ESN_glob, ...
            {'Color', clrs(it, :), 'LineWidth', 3}, 1);
        f_Set_scales({'logx'})
        
        LineHandles(end + 1) =  hndls.mainLine;
    %semilogx(ResSuppMean_1 , tFuncInd.Simpson_ESN_glob, 'o', 'MarkerFaceColor', clrs(it, :), 'MarkerEdgeColor', clrs(it, :));
    hold on
    %nonlinfun = @(beta, x) beta(1) * (exp(-x/beta(2))).*x + beta(3) * x./(x + 30);
%     nonlinfun = @(beta, x) beta(1) * x./(beta(2) * x.^2 + beta(3) * x + 1);
%     ind = ~isnan(tFuncInd.Simpson_ESN_glob);
%     nlm = fitnlm(ResSuppMean_1(ind)' , tFuncInd.Simpson_ESN_glob(ind), nonlinfun, [1,  1, 1]);
%     res = (1:5000)';
%     semilogx(res, predict(nlm, res), 'LineWidth', 2, 'Color', clrs(it, :));
    leg{it} = SimulatIDs{it};
    %leg{2*it} = '';
    f_Lbls('Mean resource, supplied', 'Biodiversity')
    hold on
    xlim([10, 5000])
    
    ax2 = subplot(3, 1, 2);
    plot(Resavailable_mean, tFuncInd.Simpson_ESN_glob, 'o', 'MarkerFaceColor', clrs(it, :), 'MarkerEdgeColor', clrs(it, :));
    hold on
    f_Lbls('Mean resource, available', 'Biodiversity')
    xlim([0, 20])
    
    ax3 = subplot(3, 1, 3);
    semilogx(BiomMean_1, tFuncInd.Simpson_ESN_glob, 'o', 'MarkerFaceColor', clrs(it, :), 'MarkerEdgeColor', clrs(it, :));
    xlim([1, 1000])
    hold on
    f_Lbls('Mean biomass', 'Biodiversity')
    nonlinfun = @(beta, x) beta(1) * x./(beta(2) * x.^2 + beta(3) * x + 1);
    nlm = fitnlm(BiomMean_1 , tFuncInd.Simpson_ESN_glob, nonlinfun, [1,  1, 1]);
    res = (1:1:1000)';
    plot(res, predict(nlm, res), 'LineWidth', 2, 'Color', clrs(it, :));
    leg{2*it - 1} = SimulatIDs{it};
    leg{2*it} = '';
    f_Lbls('Mean biomass', 'Biodiversity')
    hold on
    %    xlim([0, 180])
    
  
    
    
end

subplot(3, 1, 1);
legend(LineHandles, leg_alpha, 'Location', 'northeastoutside');
drawnow
ax2.Position(3) = ax1.Position(3);
ax3.Position(3) = ax1.Position(3);

figure(fg1);
    subplot(2, 2, 4);
 legend( leg_alpha, 'Location', 'Best');
 
 
figure(fg3)
subplot(2, 4, 4)
legend( leg_muS, 'Location', 'Best');

 
 %%plot N:P ratio distribution
 
 
 