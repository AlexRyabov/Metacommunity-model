function iFgNrNew = DiffMaps_Adds_Calculations(PlotID, iFgNr, Data)

if iFgNr > 0
%create a figure if it does not persist otherwiese clean it
fg = findobj( 'Type', 'Figure', 'Name', PlotID );
if isempty(fg)
    fg = figure(iFgNr);
    fg.Name = PlotID;
else
    figure(fg);
end
iFgNrNew = iFgNr + 1;
clf
else 
    iFgNrNew = NaN;
end

switch PlotID


    case 'SimilarityLearningCurve'
        %% SimilarityLearningCurve plot learning curve for similarity matrix 
        k_max = Data.k_max;
        idxParam = Data.idxParam;
        SpUnStacked = Data.SpUnStacked;
        idxParamUnique = unique(idxParam);
        DistanceMetrics = {'Spearman', 'Pearson', 'NormzdEuc', 'CosSim'};
        DistanceMetricsNames = {'Spearman', 'Pearson', 'Euclidean, norm', 'Cosine similarity'};
        m_max = length(idxParamUnique);
        Ms = unique(ceil(logspace(log10(10), log10(m_max), 15)));
        [RelError] = f_DM_SimilartyTest(SpUnStacked, DistanceMetrics, k_max, Ms, idxParam);
        fg.Position = [10 10  400 400];
        %%
        clf
        clrs = f_Clrs_fresh(length(DistanceMetrics));
        for iDM = 1:length(DistanceMetrics)
            %loglog(Ms, RelError(:, iDM), '-o', 'Color', clrs(iDM, :), 'MarkerFaceColor', clrs(iDM, :));
            semilogx(Ms, RelError(:, iDM), '-o', 'Color', clrs(iDM, :), 'MarkerFaceColor', clrs(iDM, :));
            hold on
        end
        ylim([0, 0.9])
        %
        legend(DistanceMetricsNames, 'Location', 'Best')
        f_Lbls('Number of grids, m', 'cor(S^{(m)}, S^{(n)})')
        %
        
         lm = fitlm(log10(Ms(4:end)), log10(1-RelError(4:end, 1)));
%         M_0 = 10^(-lm.Coefficients.Estimate(1)/lm.Coefficients.Estimate(2));
%         xmin = 10;
%         xmax =  10^ceil(log10(M_0));
%         xlim([10,2500]);
%         ylim([0, 0.25])
%         xax = logspace(1, log10(M_0));
        hold on
        %p = plot(xax, 10.^lm.predict(log10(xax')),'--', 'Color', clrs(1, :));
        %xax = logspace(1, log10(800));
        %p = plot(xax, 3*xax.^(-0.5),'--', 'Color', 'k', 'LineWidth', 2);
        %p.Annotation.LegendInformation.IconDisplayStyle = 'off'
        %hold off
        %%
    otherwise
        error('Plot not defined in DiffMaps_Adds_Others')
end



