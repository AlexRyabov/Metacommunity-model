function iFgNrNew = DiffMaps_Adds_Plots(PlotID, iFgNr, Data)

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




switch PlotID
    case 'FigSI_DiffMapThresholds'
        %% Plot diffusion map for varios thresholds
        SpUnStacked = Data.SpUnStacked;
        RStar_All = Data.RStar_All;
        Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
        Cl = Mix;
        K_maxs = [1, 3,  5, 7, 10,  30, 50, 100,  200];
        fg.Position = [10 10  710 699];
        clf;
        tiledlayout('flow');
        for i = 1:length(K_maxs)
            k_max = K_maxs(i);
            [ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'Spearman', false, k_max);
            eEV =aEV(indPos, :);
            nexttile
            scatter(eEV(:, 1), -eEV(:, 2), 30, Cl, 'filled');
            xlabel('i-trait 1');
            ylabel('i-trait 2');
            f_FramePlotMy
            title([NS(k_max) ' trusted links'])
        end
        AddLetters2Plots('HShift', -0.06, 'VShift', -0.02);

    case 'FigSI_Sample&Map'
        %
        % %plot itrats as in Fig 3
        eEV = Data.aEV;
        RStar_All = Data.RStar_All;
        tResults = Data.tResults;
        
        fg.Position = [10 10  560 482];
        subplot(2, 2, 1);
        iSmpl = 4;
        sploti = 1;
        ip = 17
        TSpan = tResults.TSpan(ip, :);
        ind = TSpan >1100;
        Biom = tResults.GlobalBiomass_dyn{ip};
        plot(TSpan(ind), Biom(:, ind));
        f_Lbls('Time', 'Biomass');
        f_FramePlotMy
        
        subplot(2, 2, 2);
        Res = tResults.GlobalRes_dyn{ip};
        plot(TSpan(ind), Res(:, ind));
        f_Lbls('Time', 'Resouces');
        ylim([0, 7000]);
        legend({'R_1', 'R_2', 'R_3'}, 'Location', 'best')
        f_FramePlotMy
        
        
        % %plot R* values as in Fig 3
        %
        Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
        Cl = Mix;
        %
        
        ax1 = subplot(2, 2, 3);
        scatter3(RStar_All(:, 1), RStar_All(:, 2), RStar_All(:, 3), 20, Cl , 'filled');
        f_Lbls3D('R^*_1', 'R^*_2', 'R^*_3')
        view(ax1,[133.992 30]);
        f_FramePlotMy
        %         %
        subplot(2, 2, 4);
        %                scatter3(eEV(:, 1), -eEV(:, 2), eEV(:, 3), 30, Cl, 'filled');
        scatter(eEV(:, 1), -eEV(:, 2), 30, Cl, 'filled');
        xlabel('i-trait 1');
        ylabel('i-trait 2');
        zlabel('i-trait 3');
        f_FramePlotMy
        AddLetters2Plots
    
        %         subplot(2, 2, 3);
        %         %scatter(eEV(:, 1), eEV(:, 2), 20, Cl, 'filled');
        %         scatter3(eEV(:, 1), eEV(:, 2), eEV(:, 3), 20, Cl, 'filled');
        %         f_Lbls('\lambda_1^{-1}v_1', '\lambda_2^{-1}v_2')
        %
        %
        
    case 'SuppF3_DistanceCorr'
        %% SuppF3 plot correlations between distances in PCA space and in diffusion map
        % space. Supplimentary figure 3
        aS_PCA = Data.aS_PCA; %similarity matrix used for PCA
        D2_PCA = Data.D2_PCA;  %Distance in PCA space scaling 2 (V2 = sqrt(lambda)*V), each species is a row in matrix of eigenvectors
        aS = Data.aS;  %%Similarity matrix used for diffusion map, the same as as_PCA, but only 10 largest elements in each column >0
        DiffDist = Data.DiffDist; % diffusion distance;
        RStar_All = Data.RStar_All;  %R* values of species
        %set figure size
        fg.Position = [150 50  1152 315];
        clf
        
        %S_PCA vs distance in R*
        DRSt = squareform(pdist(RStar_All));
        
        MrkSz = 10;
        Ys = {aS_PCA(:), D2_PCA(:), DiffDist(:)};
        YLabels = {'Spearman correlation', 'PCA distances', 'Diffusion distances'};
        
        for iSP = 1:3
            subplot(1, 3, iSP)
            Y = Ys{iSP};
            %modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
            %lm = fitnlm(DRSt(:), Y(:), modelfun, [1, 1, -1]);
            %modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2 + b(4) * x.^3;
            %lm = fitnlm(DRSt(:), Y(:), modelfun, [1, 1, -1, 1]);
            modelfun = @(b,x) b(1)* x./(b(2) + x) + b(3);
            lm = fitnlm(DRSt(:), Y(:), modelfun, [1, 1, 1]);
            x = linspace(min(DRSt(:)), max(DRSt(:)), 100);
            plot(x, lm.predict(x'), '--', 'Color', [1,1,1]*0.25,'LineWidth', 2)
            hold on
            s2 = scatter(DRSt(:), Y(:), MrkSz, 'filled');
            s2.MarkerFaceAlpha = 0.01;
            % title(['R^2=' num2str(lm.Rsquared.Adjusted, 2)])
            xlabel('Distance, R* ')
            ylabel(YLabels{iSP});
            xlim([-0.5, 12])
        end
        
    case 'SuppF4_EigVectorSpectr'
        %% SuppF4 plot absolute values of eigenvectors for PCA and DiffusionMap
        % Supplimentary figure 4
        aEV = Data.aEV; %eigenvectors in Diff map space
        aEV2_PCA = Data.aEV2_PCA;  %eigenvectors in PCA space
        %set figure size
        fg.Position = [150 50  932 373];
        clf
        
        
        MrkSz = 10;
        Ys = {sqrt(sum(aEV.^2, 1)), sqrt(sum(aEV2_PCA.^2, 1))};
        YLabels = {'\lambda_i^{-1}', '\lambda_i^{0.5}'};
        legs = {'Diffusion map', 'PCA'};
        for iSP = 1:2
            subplot(1, 2, iSP)
            Y = Ys{iSP};
            IDs = 1:min(100, length(Y(:)));
            s2 = stem(IDs, Y(IDs), 'o', 'filled');
            xlabel('i')
            ylabel(YLabels{iSP});
            title(legs{iSP});
        end
        AddLetters2Plots('HShift', -0.08);
    case 'SuppF5_FractalDimension'
        %% SuppF5 Find fractal dimension
        % Supplimentary figure 5
        aEV = Data.aEV; %eigenvectors in Diff map space
        aEV2_PCA = Data.aEV2_PCA;  %eigenvectors in PCA space
        RStar_All = Data.RStar_All;  %species traits
        %set figure size
        fg.Position = [150 50  400 400];
        %%
        clf
        ev_ind = 1:3;
        %Dimension = f_FractalDimension(RStar_All, 'boxcounting2');
        DimRSt = f_FractalDimension(RStar_All, 'correlation', 'Epsilon', logspace(log10(0.1), log10(5), 50));
        %DimDM = f_FractalDimension(aEV(:, ev_ind), 'correlation', 'Epsilon', logspace(log10(0.15), log10(7), 50));
        DimDM = f_FractalDimension(aEV(:, ev_ind), 'correlation', 'Epsilon', logspace(log10(0.1), log10(8), 50));
        %DimPCA = f_FractalDimension(aEV2_PCA(:, ev_ind), 'correlation', 'Epsilon', logspace(log10(0.1), log10(0.7), 50));
        DimPCA = f_FractalDimension(aEV2_PCA(:, ev_ind), 'correlation', 'Epsilon', logspace(log10(0.1), log10(0.8), 50));
        set(gca, 'XScale','log', 'YScale','log')
        legend({'R*, data', ['R*, fit, d=' NS(DimRSt.Estimate, 2) '\pm' NS(DimRSt.SE, 1)],...
            'Diff map, data', ['Diff map, fit, d=' NS(DimDM.Estimate, 2) '\pm' NS(DimDM.SE, 1)],...
            'PCA, data', ['PCA, fit, d=' NS(DimPCA.Estimate, 2) '\pm' NS(DimPCA.SE, 1)]}, ...
            'Location', 'Best');
        f_Lbls('\epsilon', 'Correlation sum, C(\epsilon)')
        %%
        
    case 'Fig1'  %%Old version of Fig 1 with model data
        %% Fig 1
        eEV = Data.aEV;
        ev = Data.ev;
        aS = Data.aS;
        RStar_All = Data.RStar_All;
        
        Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
        Cl = Mix;
        
        %plot figure 1
        %(A) R*, (B) hist of outcomes
        %(C) eigen vator 1+2   (D) similar matrix
        
        SpUnStacked = Data.SpUnStacked;
        fg.Position = [10 10  660 500];
        ax1 = subplot(2, 2, 1);
        scatter3(RStar_All(:, 1), RStar_All(:, 2), RStar_All(:, 3), 20, Cl , 'filled');
        f_Lbls3D('R^*_1', 'R^*_2', 'R^*_3')
        view(ax1,[133.992 30]);
        
        %%hist of outcomes
        subplot(4, 2, 2);
        ind = SpUnStacked > 0.1;
        Richness = sum(ind, 1);
        [~, ind] = sort(Richness, 'desc');
        subplot(4, 2, 2);
        bar(SpUnStacked(:, ind(1)));
        f_Set_scales({'x', 'logy'});
        f_Lbls('', 'Biomass')
        f_AdjTickLbls('Y', 5)
        ax41 = subplot(4, 2, 4);
        if size(SpUnStacked, 2) > 8995
            bar(SpUnStacked(:, 8995));
        else
            bar(SpUnStacked(:, 50));
        end
        f_Set_scales({'x', 'logy'});
        f_AdjTickLbls('Y', 5)
        f_Lbls('Species ID', 'Biomass')
        
        %
        ax4 = subplot(2, 2, 4);
        [~, indEV2] = sort(eEV(:, 1));  %sort along the first eigenvector
        pcolor(aS(indEV2, indEV2));
        cm = jet;
        cm = cm(end:-1:1, :);
        cm(1, :) = [1, 1, 1];
        colormap(cm);
        colorbar
        shading flat
        f_Lbls('Species ID', 'Species ID')
        ax4.Position(3) = ax41.Position(3);
        
        subplot(2, 2, 3);
        %scatter(eEV(:, 1), eEV(:, 2), 20, Cl, 'filled');
        scatter3(eEV(:, 1), eEV(:, 2), eEV(:, 3), 20, Cl, 'filled');
        f_Lbls('\lambda_1^{-1}v_1', '\lambda_2^{-1}v_2')
        %%
    case {'Fig3', 'Fig3_3D'}
        %% Fig 3: RSt space, eigenvalue space and functional diversity
        fg.Position = [50 50  1152 315];
        eEV = Data.aEV;
        RStar_All = Data.RStar_All;
        
        Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
        Cl = Mix;
        
        clf
        %plot figure 3
        %(A) R*, (B) eigen vector 1+2+3
        %(C) functional diversity
        
        
        ax1 = subplot(1, 3, 1);
        scatter3(RStar_All(:, 1), RStar_All(:, 2), RStar_All(:, 3), 30, Cl , 'filled');
        f_Lbls3D('R^*_1', 'R^*_2', 'R^*_3')
        view(ax1,[133.992 30]);
        f_FramePlotMy;
        
        
        subplot(1, 3, 2);
        switch PlotID
            case 'Fig3_3D'
                scatter3(eEV(:, 1), -eEV(:, 2), eEV(:, 3), 30, Cl, 'filled');
                zlabel('i-trait 3');
            otherwise
                scatter(eEV(:, 1), -eEV(:, 2), 30, Cl, 'filled');
        end
        xlabel('i-trait 1');
        ylabel('i-trait 2');
        f_FramePlotMy
        
        subplot(1, 3, 3);
        tResultsVal = Data.tResultsVal;
        
        Clrs = lines(size(tResultsVal, 1));
        Clrs3 = lines(3);
        vLocalObrv = [];
        vLocalPred = [];
        for ip = 1:size(tResultsVal, 1)
            %s2 = scatter(tResultsVal.FuncDivRaoRstSample{ip}, tResultsVal.FuncDivRaoSample{ip}, 10, Clrs(ip, :), 'filled');
            s2 = scatter(tResultsVal.FuncDivRaoRstSample{ip}, tResultsVal.FuncDivRaoSample{ip}, 10, Clrs3(1, :), 'filled');
            %s2 = scatter(tResultsVal.FuncDivRaoRstSample{ip}, tResultsVal.FuncDivRaoSample{ip}, 5, [0, 130, 200]/256, 'filled');
            s2.MarkerFaceAlpha = 0.1;
            vLocalObrv = [vLocalObrv; tResultsVal.FuncDivRaoRstSample{ip}];
            vLocalPred = [vLocalPred; tResultsVal.FuncDivRaoSample{ip}];
            hold on
        end
        plot(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', Clrs3(3, :), 'LineWidth', 1);
        
        %plot(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 0.5*[1,1,1]);
        %modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2+ b(4) * x.^3;
        %lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1, -1, 1])
        modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
        lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1, -1]);
        %modelfun = @(b,x) b(1) + b(2) * x ;
        %lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1]);
        x = linspace(min(tResultsVal.FuncDivRaoRst), max(tResultsVal.FuncDivRaoRst), 100);
        plot(x, lm.predict(x'), 'Color', [210, 245, 60]/256, 'LineWidth', 3)
        title(['R^2=' num2str(lm.Rsquared.Adjusted, 2)])
        xlabel('Fdiv Rao, ground truth ')
        ylabel('Fdiv Rao, diffusion map')
        f_FramePlotMy
        
        modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
        lm_loc = fitnlm(vLocalObrv, vLocalPred, modelfun, [1, 1, -1]);
        AddLetters2Plots(fg,  {'A', 'B', 'C'})

        %%
    case 'Fig5'
        %% Fig 5 effect of the training set
        fg.Position = [10 10 1200 1200];
        MetaDatas = Data.MetaDatas;
        %make tiles 3x3
        tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact')
        Clrs3 = lines(3);
        
        for iD = 1:length(MetaDatas)
            %% for each Data cell
            ModData = MetaDatas{iD};
            %Get values
            aEV = ModData.aEV;
            RStar_All = ModData.RStar_All;
            tResultsVal = ModData.tResultsVal;
            DiffDist = ModData.DiffDist;
            Y = DiffDist(:);
            
            
            %plot eigevectors
            ax1 = nexttile(iD);
            Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
            Cl = Mix;
            scatter(aEV(:, 1), aEV(:, 2), 30, Cl, 'filled');
            xlabel('Trait_1');
            if iD == 1
                ylabel('Trait_2');
            end
            
            
            %plot Distance
            nexttile(iD + 3);
            DRSt = squareform(pdist(RStar_All)); %distance in R*
            MrkSz = 10;
            modelfun = @(b,x) b(1)* x./(b(2) + x) + b(3);
            lm = fitnlm(DRSt(:), Y(:), modelfun, [1, 1, 1]);
            x = linspace(min(DRSt(:)), max(DRSt(:)), 100);
            plot(x, lm.predict(x'), '--', 'Color', [1,1,1]*0.25,'LineWidth', 2)
            hold on
            s2 = scatter(DRSt(:), Y(:), MrkSz, Clrs3(1, :), 'filled');
            s2.MarkerFaceAlpha = 0.01;
            xlabel('Distance, R* ')
            if iD == 1
                ylabel('Diffusion distances');
            end
            xlim([-0.5, 12])
            title(['R^2' NS(lm.Rsquared.Adjusted, 2)])
            
            %plot functional diversity
            ax1 = nexttile(iD + 6);
            for ip = 1:size(tResultsVal, 1)
                s2 = scatter(tResultsVal.FuncDivRaoRstSample{ip}, tResultsVal.FuncDivRaoSample{ip}, 10, Clrs3(1, :), 'filled');
                s2.MarkerFaceAlpha = 0.1;
                hold on
            end
            plot(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', Clrs3(3, :), 'LineWidth', 1);
            
            modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
            lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1, -1]);
            x = linspace(min(tResultsVal.FuncDivRaoRst), max(tResultsVal.FuncDivRaoRst), 100);
            plot(x, lm.predict(x'), 'Color', [210, 245, 60]/256, 'LineWidth', 3)
            title(['R^2=' num2str(lm.Rsquared.Adjusted, 2)])
            xlabel('Fdiv Rao, R* ')
            if iD == 1
                ylabel('Fdiv Rao, Diffusion map')
            end
            
            
        end
        
    case 'EigenVectorsAndR*'
        %% plot Eig vectors and colored  R*1, R*2
        fg.Position = [10 10 900 600];
        aEV = Data.aEV;
        RStar_All = Data.RStar_All;
        %% eig vectors R*1, 2 R2 as color
        subplot(2, 2, 1);
        colormap(jet);
        scatter(aEV(:, 1), aEV(:, 2), 20, RStar_All(:, 1), 'filled');
        cb = colorbar;
        f_Lbls('ev_1', 'ev_2');
        cb.Label.String = 'R^*_1';
        subplot(2, 2, 2);
        colormap(jet);
        scatter(aEV(:, 1), aEV(:, 2), 20, RStar_All(:, 2), 'filled');
        f_Lbls('ev_1', 'ev_2')
        cb = colorbar;
        cb.Label.String = 'R^*_2';
        
        subplot(2, 2, 3)
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, RStar_All(:, 1), 'filled');
        f_Lbls3D('ev_1', 'ev_2', 'ev_3')
        cb = colorbar;
        cb.Label.String = 'R^*_1';
        
        subplot(2, 2, 4)
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, RStar_All(:, 2), 'filled');
        f_Lbls3D('ev_1', 'ev_2', 'ev_3')
        cb = colorbar;
        cb.Label.String = 'R^*_2';
    case 'EigenVectorsAndRStAsColor' 
        %% plot 2 first eigenvectors and R* as color (same as in Fig. 1)
        fg.Position = [10 10 400 400];
        eEV = Data.aEV;
        RStar_All = Data.RStar_All;
        
        Mix =  (RStar_All(:, :)/max(RStar_All(:))).^1.3;
        Cl = Mix;
        scatter(eEV(:, 1), eEV(:, 2), 20, Cl, 'filled');
        f_Lbls('\lambda_1^{-1}v_1', '\lambda_2^{-1}v_2')
        
    case '3FirstEigenVectors'
        fg.Position = [10 10 400 400];
        aEV = Data.aEV;
        colormap(jet)
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3),20, aEV(:, 3), 'o', 'filled')
        f_Lbls3D('\lambda_1^{-1}v_1', '\lambda_2^{-1}v_2', '\lambda_3^{-1}v_3')
        c = colorbar;
        c.Label.String = '\lambda_3^{-1}v_3';
        c.Location = 'northoutside';
    case 'EigenVectorsDistanceRSt'
        %plot Eig vectors colored with distance from centrer
        fg.Position = [10 10 900 600];
        aEV = Data.aEV;
        RStar_All = Data.RStar_All;
        %% eig vectors R*1, 2 R2 as color
        subplot(2, 2, 1);
        colormap(jet);
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, pdist2(RStar_All(:, :), mean(RStar_All, 1)), 'filled');
        cb = colorbar;
        f_Lbls3D('ev_1', 'ev_2', 'ev_3');
        cb.Label.String = 'Distance from center in R^*';
        
        subplot(2, 2, 3);
        colormap(jet);
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, pdist2(RStar_All(:, :), mean(RStar_All, 1)), 'filled');
        view(90, 0)
        cb = colorbar;
        f_Lbls3D('ev_1', 'ev_2', 'ev_3');
        cb.Label.String = 'Distance from center in R^*';
        
        subplot(2, 2, 2);
        colormap(jet);
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, pdist2(RStar_All(:, :), mean(RStar_All, 1)), 'filled');
        view(0, 0)
        cb = colorbar;
        f_Lbls3D('ev_1', 'ev_2', 'ev_3');
        cb.Label.String = 'Distance from center in R^*';
        
        subplot(2, 2, 4);
        colormap(jet);
        scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, pdist2(RStar_All(:, :), mean(RStar_All, 1)), 'filled');
        view(0, 90)
        cb = colorbar;
        f_Lbls3D('ev_1', 'ev_2', 'ev_3');
        cb.Label.String = 'Distance from center in R^*';
    case 'DoninancePerGrid_R*,EigVec'
        %% show the most dominant species in space of R* and eigvectors
        fg.Position = [10 10 897.0000  377];
        aEV = Data.aEV;
        PlAbnd = Data.PlAnd;
        meanAbnd = mean(PlAbnd, 2);
        RStar_All = Data.RStar_All;
        subplot(1, 2, 1);
        scatter(aEV(:, 1), aEV(:, 2), 10, [1,1,1]*0.75, 'filled')
        hold on
        ind =  meanAbnd(:)>0;
        scatter(aEV(ind, 1), aEV(ind, 2), meanAbnd(ind), 'r', 'filled')
        f_Lbls3D('ev_1', 'ev_2')
        
        subplot(1, 2, 2);
        scatter(RStar_All(:, 1), RStar_All(:, 2), 10, [1,1,1]*0.75, 'filled')
        hold on
        ind =  meanAbnd(:)>0;
        scatter(RStar_All(ind, 1), RStar_All(ind, 2), meanAbnd(ind), 'r', 'filled')
        f_Lbls3D('R^*_1', 'R^*_2')
    case 'MakePCAPlot'
        %plot the length of eigenvectors
        fg.Position = [10 10 900 500];
        aEV = Data.aEV;
        subplot(1, 2, 1);
        plot(sqrt(sum((aEV).^2, 1)));
        %make PCA
        
        %         a = rand(10, 1);
        %         b = a + 0.1 * rand(10, 1);
        %         plot(a, b, '.')
        [coeff,score,latent,tsquared,explained,mu] = pca(aEV)
        subplot(1, 2, 2);
        plot(explained);
        %plot the lenght of vector in PCA
    case 'ClusterFuncDiv'
        %make clusters in R* map them into clusters in eigenvalues
        %find functional diversity in R* and eigenvalues.
        ClNum = [1:20];
        fg.Position = [10 10 900 900];
        ax1 = subplot(2, 2, 1);
        ax2 = subplot(2, 2, 2);
        ax3 = subplot(2, 2, 3);
        
        aEV = Data.aEV;
        RStar_All = Data.RStar_All;
        DiffDist = Data.DiffDist;
        
        for iCl = 1:length(ClNum)
            clrs = jet(iCl);
            %make clustering
            idx = kmeans(aEV(:,1:4), iCl);
            subplot(2, 2, 1);
            gscatter(RStar_All(:, 1), RStar_All(:, 2), idx, clrs);
            f_Lbls('R^*_1', 'R^*_2')
            subplot(2, 2, 2);
            
            gscatter(aEV(:, 1), aEV(:, 2), idx, clrs);
            f_Lbls('\lambda_1\Psi_1', '\lambda_2\Psi_2')
            %make a loop over clusters and find functional diveristy
            PlAnd = ones(size(RStar_All, 1), 1);
            subplot(2, 2, 3)
            %cla
            for ik = 1:iCl
                indPl =  idx== ik;
                if sum(indPl) > 3  %cluster includes more than 3 species
                    FuncDivRao = f_DiversityMetrics(PlAnd(indPl), aEV(indPl, :), 'RaoDist',  DiffDist(indPl, indPl));
                    FuncDivHull90 = f_DiversityMetrics(PlAnd(indPl), aEV(indPl, 1:2), 'Hull90');
                    FuncDivRaoRst =    f_DiversityMetrics(PlAnd(indPl),RStar_All(indPl, :), 'Rao');
                    FuncDivHull90Rst = f_DiversityMetrics(PlAnd(indPl), RStar_All(indPl, 1:2), 'Hull90');
                    
                    %plot func diversity
                    subplot(2, 2, 3);
                    plot(FuncDivRaoRst, FuncDivRao, '.', 'MarkerSize', 20, 'Color', clrs(ik, :));
                    f_Lbls('Fdiv Rao, R* ', 'Fdiv Rao, DiffMap')
                    hold on
                    subplot(2, 2, 4);
                    plot(FuncDivHull90Rst, FuncDivHull90, '.', 'MarkerSize', 20, 'Color', clrs(ik, :));
                    f_Lbls('95%-hull side, R^*', '95%-hull side, diff map')
                    hold on
                    drawnow
                end
            end
            % pause
        end
    otherwise
        error('Plot not defined in DiffMaps_Adds_Others')
end

%plot Eig vectors and colored with distance from centrer
% and ev1, ev2, R*1, R*2
end

