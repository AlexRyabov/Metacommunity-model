%% Modeling spatial resource competiton
%% Alexey Ryabov 2020
%% processing results of ith diffusion maps



%% plot species distributios
% Nfin = tResults.Nfin{1};
% %Nfin = N;
% figure(1)
% for i = 1:size(Nfin, 3)
%     pcolor(Nfin(:, :, i));
%     shading flat
%     colorbar
%     drawnow
%     caxis([0, max([1, caxis])])
%     title(NS(i))
%     pause
% end

%% Run through tResults and mark all szmulations where only one species survived
%% We cannot use this simulation for comparison, in particular in this case we
%% cannot normilize diveding by variance in a sample

ID2stay = false(length(Params), 1);
for ID = 1:length(Params)
    GlobalBiomass_dyn = tResults.GlobalBiomass_dyn{ID};
    %Get number of species in tResutls for given ID
    Nsp_ID = sum(GlobalBiomass_dyn(:, end) > 1e-1); % count only species with total biomass > 0.1;
    %if one species in tResults mark this record
    ID2stay(ID) = Nsp_ID > 2;
end
Params = Params(ID2stay);
tResults = tResults(ID2stay, :);
tFuncInd = tFuncInd(ID2stay, :);

%%
iFgNr = 50;  %for addign new figures

fg01 =    f_MakeFigure(1, [10, 10, 600, 600]);
tDifMapDiv = table();

%'Params', 'tParams', 'tResults'
%create a list of replicas, and different experiments
%and make a loop over this list
SpUnStacked = [];
idxParam = [];
for ID = 1:length(Params)
    ID
    tResults.ID(ID) = ID;  %%Add an ID for each simulation
    %Create a table with species biomasses for the given replica
    N_SpaceTimeDyn = tResults.N_SpaceTimeDyn{ID};
    GlobalBiomass_dyn = tResults.GlobalBiomass_dyn{ID};
    
    %ix = GlobalBiomass_dyn(:, end) > 0.1;
    %sum(ix)
    %Create a tables with species paramteres Add species ID into this table
    RSt = tParams.RStar{ID};
    RSt1 = RSt(:, 1);
    RSt2 = RSt(:, 2);
    SpecID = tParams.specID(ID, :)';
    if iscell(SpecID)
        SpecID = SpecID{1:end}';
    end
    MaxGrowthRate = tParams.MaxGrowthRate{ID};
    tSpec = table(SpecID, MaxGrowthRate, RSt1 , RSt2);
    
    
    %remove species with average biomass less than 10^-5 from both tables
    %tSpec = tSpec(ix, :);
    %N_SpaceTimeDyn = N_SpaceTimeDyn(:, :, ix, :);
    %unstack the table species x (time+location) (for each moment of time unwrap locations
    SpUnStackedParam = []; %results for given paramter set
    for iSp = 1:size(N_SpaceTimeDyn, 3)
        switch SimulatID
            %include only last state (without perturb)
            case {'DeltaSRandShift', 'DeltaSRandShiftR3', 'DeltaSRandShiftR3Reg','DeltaSRandShiftR3L', 'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3Gr', 'DeltaSRandShiftR3D03', 'DeltaSLarge_R3', 'DeltaSWideR3D03'} %
                SpUnStackedBlock =  N_SpaceTimeDyn(:, :, iSp, end);
            otherwise
                %include all states except the first
                SpUnStackedBlock =  N_SpaceTimeDyn(:, :, iSp, 2:end);
        end
        SpUnStackedParam = [SpUnStackedParam, SpUnStackedBlock(:)];
    end
    
    %take 100 random results from this replica
    %     if size(SpUnStackedParam, 1) > 300
    %         indSmpl = randperm(size(SpUnStackedParam, 1), 300);
    %     else
    %         indSmpl = 1:size(SpUnStackedParam, 1);
    %     end
    %and transpose the table)
    %SpUnStackedParam = SpUnStackedParam(indSmpl, :)';
    SpUnStackedParam = SpUnStackedParam(:, :)';
    SpUnStackedParamAllSp = zeros(length(Params(ID).SpecID_All), size(SpUnStackedParam, 2));
    %SpUnStackedParamAllSp(Params(ID).specID(ix), :) = SpUnStackedParam;
    SpUnStackedParamAllSp = SpUnStackedParam;
    SpUnStacked = [SpUnStacked, SpUnStackedParamAllSp];
    Params(ID).ID = ID;
    idxParam = [idxParam, ones(1, size(SpUnStackedParam, 2))*ID]; %list of ParamIndexes corresponding to columns in SpUnStacked
end

SpUnStacked(SpUnStacked <0.01) = 0;

%make validation sets tResultsVal, ParamsVal, SpUnStackedVal
if MainModelParams.Validate < 1  %what fraction of data will used for diffusion map
    m = ceil(length(idxParam)*MainModelParams.Validate);
    idxParamVal = idxParam(end-m+1:end);
    m = ceil(size(SpUnStacked, 2)*MainModelParams.Validate);
    SpUnStackedVal = SpUnStacked(:, end-m+1:end);
    m = ceil(length(Params)*MainModelParams.Validate);
    ParamsVal = Params(end-m+1:end);
    tResultsVal = tResults(end-m+1:end, :);
else
    SpUnStackedVal = SpUnStacked;
    idxParamVal = idxParam;
    ParamsVal = Params;
    tResultsVal = tResults;
end

%put into tResults, Params and SpUnStacked only the training part
if MainModelParams.Train < 1 %what fraction of data will used for diffusion map
    m = floor(length(idxParam)*MainModelParams.Train);
    idxParam = idxParam(1:m);
    m = floor(size(SpUnStacked, 2)*MainModelParams.Train);
    SpUnStacked = SpUnStacked(:, 1:m);
    m = floor(length(Params)*MainModelParams.Train);
    Params = Params(1:m);
    tResults = tResults(1:m, :);
end


%%
figure
pcolor((SpUnStacked))
shading flat
cm = gray;
cm = cm(end:-1:1, :);
colormap(cm);
ca = caxis;
caxis([0, 0.1*ca(2)])
%%
%% make diffusion map
k_max = 10;
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'StndzdEuc', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'StndzClmnsdEuc', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'VectProd', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'NormzdEuc', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'StndzClmnsdEuc', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'NormzdGaus', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'NormzdSpearman', true, k_max);
[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'Spearman', true, k_max);
%[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked, 'NormzdPearson', true, k_max);

%Results assigned for the structure, which will be returned by
Data.aEV = aEV;
Data.indPos = indPos;
Data.RStar_All = Params(1).RStar_All;
Data.DiffDist = DiffDist;

%plot Diffusion map for varios thresholds 
iFgNr = DiffMaps_Adds_Plots('FigSI_DiffMapThresholds', iFgNr,  struct( 'SpUnStacked', SpUnStacked, 'RStar_All', Params(1).RStar_All(indPos, :)));


%Get representation in PCA space
%[aEV_PCA, aEV2_PCA, score_PCA, ev_PCA, aS_PCA, D2_PCA, explained_PCA, mu_PCA] = f_PCA(SpUnStacked'); indPos_PCA = 1:size(ev);
%aEV = aEV2_PCA; DiffDist = D2_PCA;  %to use PCA distances

%%Plot diffusion distance vs species traits
DiffMaps_Adds_SimTraits(DiffDist(indPos, indPos), Params(1).RStar_All(indPos, :), figure(8))

%%Plot similarity vs species traits
%DiffMaps_Adds_SimTraits(aS, Params(1).RStar_All(indPos, :), figure(8))

% plot learning curve for similarity matrix
%iFgNr = DiffMaps_Adds_Calculations('SimilarityLearningCurve', iFgNr,  struct( 'k_max', k_max, 'idxParam', idxParam, 'SpUnStacked', SpUnStacked));

% plot correlations between distances in PCA space and in diffusion map
% space. Supplimentary figure 3
%iFgNr = DiffMaps_Adds_Plots('SuppF3_DistanceCorr', iFgNr, struct('aS',aS(indPos, :), 'aS_PCA', aS_PCA, 'D2_PCA', D2_PCA, 'DiffDist', DiffDist, 'RStar_All', Params(1).RStar_All(indPos, :)));

% plot absolute values of eigenvectors for PCA and DiffusionMap
% Supplimentary figure 4
%iFgNr = DiffMaps_Adds_Plots('SuppF4_EigVectorSpectr', iFgNr, struct('aEV',aEV(indPos, :), 'aEV2_PCA', aEV2_PCA));

%% Find fractal dimension
%iFgNr = DiffMaps_Adds_Plots('SuppF5_FractalDimension', iFgNr, struct('aEV',aEV(indPos, :), 'aEV2_PCA', aEV2_PCA, 'RStar_All', Params(1).RStar_All(indPos, :)));

%% Make MDS of similatiry matrix and compare it with species traits
% DiffMaps_Adds_SimMDS(aS, Params(1).RStar_All(indPos, :), figure(9))
% Dd=1-aS;
% aEV = MDS(Dd,5);
% aEV = aEV';
% aEV_fin = zeros(size(SpUnStacked, 1), size(aEV, 2));
% aEV_fin(indPos, :) = aEV;
% aEV = aEV_fin;


%iFgNr = DiffMaps_Adds_Plots('MakePCAPlot', iFgNr, struct('aEV',aEV(indPos, :)));
%[coeff,score,latent,tsquared,explained,mu] = pca(aEV);
%aEV = coeff .* repmat(explained', size(aEV, 1), 1);

%plot Eig vectors and colored  R*1, R*2
iFgNr = DiffMaps_Adds_Plots('EigenVectorsAndR*', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :)));
%plot Eig vectors colored with distance from centrer
iFgNr = DiffMaps_Adds_Plots('EigenVectorsDistanceRSt', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :)));


%make clusters in R* map them into clusters in eigenvalues
%find functional diversity in R* and eigenvalues.
%  iFgNr = DiffMaps_Adds_Plots('ClusterFuncDiv', iFgNr, ...
%      struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), 'DiffDist', DiffDist));

iFgNr = DiffMaps_Adds_Plots('Fig1', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), 'aS', aS, 'SpUnStacked', SpUnStacked, 'ev', ev));

%plot 3 first eigenvectors in 3D
iFgNr = DiffMaps_Adds_Plots('3FirstEigenVectors', iFgNr, struct('aEV',aEV(indPos, :)));

%plot 2 first eigenvectors and R* as color (same as in Fig. 1)
iFgNr = DiffMaps_Adds_Plots('EigenVectorsAndRStAsColor', iFgNr, struct('aEV',aEV(indPos, :), 'RStar_All', Params(1).RStar_All(indPos, :)));



iFgNr = DiffMaps_Adds_Plots('Fig1', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), 'aS', aS, 'SpUnStacked', SpUnStacked, 'ev', ev));

%plot biomasses as a funciton of time 
%take them from 
% %%
%%
iFgNr = DiffMaps_Adds_Plots('FigSI_Sample&Map', iFgNr, ...
    struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), 'ev', ev, 'tResults', tResults));
%    f_exportfig_my(gcf , '.\figs\DM\SI_Fig3Oscill.png', 1.1);
        %
%% plot eigenvectors vs r* and mu
figure(1);
clf
iSupPl = 1;
for iEv = 1:4
    subplot(4, 3, iSupPl);
    hold on
    plot(aEV(indPos, iEv), Params(1).RStar_All(indPos, 1), '.');
    iSupPl = iSupPl + 1;
    f_Lbls(sprintf('ev%.f', iEv), 'R*_1');
    
    subplot(4, 3, iSupPl);
    hold on
    plot(aEV(indPos, iEv), Params(1).MaxGrowthRate_All(indPos), '.');
    f_Lbls(sprintf('ev%.f', iEv), '\mu_{max}');
    iSupPl = iSupPl + 1;
    
    subplot(4, 3, iSupPl);
    hold on
    plot(aEV(indPos, 1), aEV(indPos, iEv+ 1), '.')
    f_Lbls(sprintf('ev%.f', 1), sprintf('ev%.f', iEv+ 1));
    iSupPl = iSupPl + 1;
    
end
drawnow





%% find the resulting functional diversity for each sample
fg06 = figure(6);
fg06.Position = [10 10 400 400];
clf

%tResults ->> tResultsVal, becuase it is a validation array
%add ID there as a column
%Why if we run it for 25 training set we get only 27 in the validation???
for ip = 1:length(ParamsVal)
    %
    switch SimulatID
        case {'DeltaSRandShift', 'DeltaSRandShiftR3', 'DeltaSRandShiftR3Reg', 'DeltaSRandShiftR3L', 'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3D03', 'DeltaSRandShiftR3Gr', 'DeltaSLarge_R3', 'DeltaSWideR3D03', 'DeltaSRandShiftR3L1SpecDyn'}  %1D trait space
            TraitsEV = aEV(:, 1:20);
            TraitsRst = ParamsVal(ip).RStar_All;
        case {'DeltaSConstGlOpportResh', 'DeltaSRandShiftGlOpportResh'}
            TraitsEV = aEV(:, 1:4);
            TraitsRst = [ParamsVal(ip).RStar_All, 5*(ParamsVal(ip).MaxGrowthRate)];
        otherwise
            error('Traits was not defined');
    end
    
    indParsSel = find(idxParamVal == ParamsVal(ip).ID);
    %indSmpl = randperm(length(indParsSel), 100);
    indSmpl = indParsSel;% (indSmpl);
    tFDiv = table();
    warning('off','all');
    for i = 1:length(indSmpl)
        PlAnd = SpUnStackedVal(:, indSmpl(i));
        indPl = PlAnd>0;
        %tFDiv.FuncDivRao(i) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'Rao',  DiffDist(indPl, indPl));
        tFDiv.FuncDivRao(i) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'RaoDist',  DiffDist(indPl, indPl));
        %tFDiv.FuncDivRao(i) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'FDDist',  DiffDist(indPl, indPl));
        %tFDiv.FuncDivHull90(i) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'Hull90');
        tFDiv.FuncDivRaoRst(i) =    f_DiversityMetrics(PlAnd(indPl),TraitsRst(indPl, :), 'Rao');
        %tFDiv.FuncDivRaoRst(i) =    f_DiversityMetrics(PlAnd(indPl),TraitsRst(indPl, :), 'FD');
        %tFDiv.FuncDivHull90Rst(i) = f_DiversityMetrics(PlAnd(indPl), TraitsRst(indPl, :), 'Hull90');
        tFDiv.SimpsonESN(i) = f_DiversityMetrics(PlAnd(indPl), [], 'SimpsonESN');
        tFDiv.Richness(i) = f_DiversityMetrics(PlAnd(indPl), [], 'Richness');
    end
    warning('on','all');
    lm = fitlm(tFDiv.FuncDivRaoRst, tFDiv.FuncDivRao);
    tResultsVal.FuncDivRaoFitEst1(ip) = lm.Coefficients.Estimate(1);
    tResultsVal.FuncDivRaoFitEst2(ip) = lm.Coefficients.Estimate(2);
    tResultsVal.FuncDivRaoFitR2(ip) = lm.Rsquared.Adjusted;
    
    %find func diversity across the gridd
    PlAnd = mean(SpUnStackedVal(:, indSmpl), 2);
    %SpUnStackedVal(:, indSmpl(i));
    indPl = PlAnd>0;
    
    %tResultsVal.FuncDivRao(ip) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'Rao',  DiffDist(indPl, indPl));
    tResultsVal.FuncDivRao(ip) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'RaoDist',  DiffDist(indPl, indPl));
    %tResultsVal.FuncDivRao(ip) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'FDDist',  DiffDist(indPl, indPl));
    %tResultsVal.FuncDivHull90(ip) = f_DiversityMetrics(PlAnd(indPl), TraitsEV(indPl, :), 'Hull90');
    tResultsVal.FuncDivRaoRst(ip) =    f_DiversityMetrics(PlAnd(indPl),TraitsRst(indPl, :), 'Rao');
    %tResultsVal.FuncDivRaoRst(ip) =    f_DiversityMetrics(PlAnd(indPl),TraitsRst(indPl, :), 'FD');
    %tResultsVal.FuncDivHull90Rst(ip) = f_DiversityMetrics(PlAnd(indPl), TraitsRst(indPl, :), 'Hull90');
    tResultsVal.SimpsonESN(ip) = f_DiversityMetrics(PlAnd(indPl), [], 'SimpsonESN');
    tResultsVal.Richness(ip) = f_DiversityMetrics(PlAnd(indPl), [], 'Richness');
    
    
    tResultsVal.FuncDivRaoSample{ip} = tFDiv.FuncDivRao;
    tResultsVal.FuncDivRaoRstSample{ip} = tFDiv.FuncDivRaoRst;
    tResultsVal.SampleCols{ip} = indSmpl';
    tResultsVal.Species{ip}= tSpec;
    tResultsVal.EigenVect{ip}= TraitsEV(:, :);
    tResultsVal.Similarity{ip}= aS;
    tResultsVal.EigenVal{ip}= ev;
    %     figure(6)
    
    %
    %     hold on
    %     pl = plot(tFDiv.FuncDivRaoRst, tFDiv.FuncDivRao, '.', 'MarkerSize', 20);
    %     hold on
    %     xlabel('Fdiv Rao, R* ')
    %     ylabel('Fdiv Rao, DiffMap')
    
    %     figure(8)
    %     subplot(1, 2, 2)
    %     p = ParamsVal(ip);
    %
    %
    %
    %     scatter(1:3, p.ResHigh, 50, 'r', 'filled');
    %     hold on
    %     scatter(1:3, p.ResLow, 50, 'b', 'filled');
    %     hold off
    %       xlabel('Resource ID');
    %     ylabel('R_{min}-R_{max}');
    
    %show the most dominant species in space of R* and eigvectors
    %     iFgNr = DiffMaps_Adds_Plots('DoninancePerGrid_R*,EigVec', iFgNr, ...
    %     struct('aEV',aEV(:, :),...
    %         'RStar_All', Params(1).RStar_All(:, :), ...
    %         'PlAnd', SpUnStackedVal(:, indSmpl)) );
    
    
    %pause
    %       drawnow
    %      figure(6)
    %      delete(pl)
    %      plot(tFDiv.FuncDivRaoRst, tFDiv.FuncDivRao, '.', 'Color', [1,1,1]*0.9);
end

%
iFgNr = DiffMaps_Adds_Plots('Fig3', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), ...
    'tResultsVal', tResultsVal));
%f_exportfig_my(gcf , '.\figs\DM\Fig3.png', 1.1);
iFgNr = DiffMaps_Adds_Plots('Fig3_3D', iFgNr, struct('aEV',aEV(indPos, :),'RStar_All', Params(1).RStar_All(indPos, :), ...
    'tResultsVal', tResultsVal));
%f_exportfig_my(gcf , '.\figs\DM\SI_Fig3in3D.png', 1.1);

%%
figure(5)
clf
subplot(2, 2, 1)
yyaxis left
plot(tResultsVal.SimpsonESN);
hold on
plot(tResultsVal.Richness);
ylabel('Diversity')
yyaxis right
plot(tResultsVal.FuncDivRao);
plot(tResultsVal.FuncDivRaoRst);
ylabel('Funcitonal diversity')
legend({'Eff sp number', 'Richness', 'FDiv, DM', 'FDiv, R*'}, 'Location', 'northwest')
xlabel('Setup ID')


figure(6)
hold on
plot(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 0.75*[1,1,1]);
hold on
%modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2+ b(4) * x.^3;
%lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1, -1, 1]);
modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
lm = fitnlm(tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRao, modelfun, [1, 1, -1]);
%[Bl, FitInfo] = lasso([ones(length(tResultsVal.FuncDivRaoRst), 1), tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRaoRst.^2, tResultsVal.FuncDivRaoRst.^3],tResultsVal.FuncDivRao, 'CV',10)
%lassoPlot(Bl,FitInfo,'PlotType','CV');
%legend('show') % Show legend ->> lambda = 0.0055
%[Bl, FitInfo] = lasso([ones(length(tResultsVal.FuncDivRaoRst), 1), tResultsVal.FuncDivRaoRst, tResultsVal.FuncDivRaoRst.^2, tResultsVal.FuncDivRaoRst.^3],tResultsVal.FuncDivRao, 'CV',10)

x = linspace(min(tResultsVal.FuncDivRaoRst), max(tResultsVal.FuncDivRaoRst), 100);
plot(x, lm.predict(x'), 'k', 'LineWidth', 3)
title(['R^2=' NS(lm.Rsquared.Adjusted)])
xlabel('Fdiv Rao, R* ')
ylabel('Fdiv Rao, diff map')

figure(7)
%plot(tResultsVal.FuncDivHull90Rst, tResultsVal.FuncDivHull90, '.', 'MarkerSize', 20);
%hold on
%modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
%lm = fitnlm(tResultsVal.FuncDivHull90Rst, tResultsVal.FuncDivHull90, modelfun, [1, 1, -1]);
%x = linspace(min(tResultsVal.FuncDivHull90Rst), max(tResultsVal.FuncDivHull90Rst), 100);
%plot(x, lm.predict(x'), 'LineWidth', 3)
%title(['R^2=' NS(lm.Rsquared.Adjusted)])
%xlabel('95%-hull side, R* ')
%ylabel('95%-hull side, diff map')


%return species paramters, eigenvectors, eigenvalues,
%functional diversities for each sample and replica
%Table with results for validate dataset
Data.tResultsVal = tResultsVal;



