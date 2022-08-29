%%
% Modeling spatial resource competiton
% Alexey Ryabov 2017
% This script makes plots figures for the modelling of effects of mean and variance in resource distribution 
% on biodiversity


%tResults
%tParams
%tFuncInd

iFgNr = 50;
%% Make a plot with isoclines 
iParam = 1;
ModParams = Params(iParam);
ModParams.RStar_All = ModParams.RStar;
%         [dx, Lx, Ly, tmin, tmax, ...
%             res, D, S, R, ResDyn, d, r, m, R_eq, c, N, ModParams] = ...
%             Model_Params (ModParams, iParam);
iFgNr = ResHeterogen_AddPlots('Fig2_ResPlane', iFgNr, struct('consumprate', ModParams.consumprate, 'RStar_All', Params(1).RStar(:, :)));


%% plot species distributsion
figure(10)
SpID = 1;
RStar = tParams.RStar{SpID};
subplot(2, 2, 1)
f_Resource_Plane_Plot(RStar, ModParams.consumprate, [1000, 1000], struct('Colors', true, 'LogScale', true, 'Number', 20));
hold on
S = tParams.S{SpID};
S1 = S(:, :, 1);
S2 = S(:, :, 2);
plot(S1(:), S2(:), '.k')
f_Lbls('Resource 1', 'Resource 2')
subplot(2, 2, 2)
f_DistribPlot( tResults.Nfin{SpID}, 1:tParams.Lx(SpID), 1:tParams.Ly(SpID), RStar(:, 1), false, 1)
title('Species dominance')
f_Lbls('X', 'Y')

subplot(2, 2, 3)
pcolor( 1:tParams.Lx(SpID), 1:tParams.Ly(SpID), tFuncInd.Simpson_ESN_loc_xy{SpID})
shading flat
colorbar
f_Lbls('X', 'Y')
        

%% Make Figure with Diversity vs Biomass 
Div_lcl = tFuncInd.Simpson_ESN_loc_avg;
Div_lcl = tFuncInd.Simpson_ESN_loc_xy{1};Div_lcl = Div_lcl(:);
Div_rgn = tFuncInd.Simpson_ESN_glob;
%Biomass = tFuncInd.BiomassTotal;
Biomass = tFuncInd.Biomass_xy{1}; Biomass= Biomass(:);

[iFgNr, fg] = ResHeterogen_AddPlots('Fign_DivBiom', iFgNr, struct('consumprate', ModParams.consumprate, 'RStar_All', Params(1).RStar(:, :)));

tiledlayout(2,3,'TileSpacing','compact')

%total biomass vs global diversity
nexttile(4)
scatter(Biomass, Div_rgn, 'o', 'filled'); f_Set_scales({'logx'}); xlim([1e3, 1e7])
f_Lbls('Total biomass', 'Regional diversity')

%total biomass vs global richness
nexttile(5)
scatter(Biomass, Div_lcl, 'o', 'filled'); f_Set_scales({'logx'}); xlim([1e3, 1e7])
f_Lbls('Total biomass', 'Local diversity')



%plot biodiversity as a function of maximal resource
fg10 = f_MakeFigure(10, [50, 50, 600, 350]);

ResSuppMean_1 = [];
ResSuppVar_1 = [];
for i = 1:size(tResults, 1)
    S = tResults.S{i};
    S = S(:, :, 1);
    ResSuppMean_1(end + 1) = mean(S(:));
    ResSuppVar_1(end + 1) = var(S(:));
end
subplot(1, 2, 1)
loglog(ResSuppMean_1, ResSuppVar_1, 'o');
f_Lbls('Mean resource', 'Variance')

subplot(1, 2, 2)
loglog(ResSuppMean_1 , tFuncInd.Simpson_ESN_glob, 'o');
f_Lbls('Mean resource', 'Biodiversity')


%%plot samples of resource distribution
fg11 = f_MakeFigure(11, [50, 50, 350, 350]);
PlotInd = unique(round(linspace(1, size(tResults, 1), 5)));
clrs = linspecer(length(PlotInd));
for i = 1:length(PlotInd);
    S = tResults.S{PlotInd(i)};
    S1 = S(:, :, 1);
    S2 = S(:, :, 2);
    loglog(S1(:), S2(:), '.', 'Color', clrs(i, :));
    hold on
end
f_Lbls('Resource_1', 'Resource_2')


%%
%% plot biomass dynamics for each experiment
rows = 4;
cols = 4;
%if rows < 7 && cols < 7
    Fg12 = figure(12);
    set(Fg12, 'Position', [20, 20, 1200, 1200]);
    sploti = 1;
    ind = unique(round(linspace(1, length(Params), 4)));
    for pi = 1:length(Params)
        if sploti > 16 
            break
        end
        subplot(rows, cols, sploti);
        Biom = tResults.GlobalBiomass_dyn{pi};
        S = tResults.S{pi};
        S = S(:, :, 1);
        BiomMax = max(Biom(:));
        plot(tResults.TSpan(pi, :), Biom);
        f_Lbls('Time', 'biomass');
        text(0, BiomMax, sprintf('\\mu(S)=%0.1f', mean(S(:))))
        sploti = sploti + 1;
        %ylim([1e-2, 1e5])
    end





