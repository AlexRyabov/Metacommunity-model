MainModelParams.ReshuffleStep = 50;
MainModelParams.TraitRangesLowHigh = [1, 10];

%Meta analysis of DeltaSRandShiftR3L1Spec to produce figure 5
Calculate = 0;
MetaTrainSets = [100, 200, 600];
MetaDatas = {};
for iTS = 1:length(MetaTrainSets)
    MainModelParams.Train     = MetaTrainSets(iTS)/688;  %what fraction of data will used for diffusion map
    MainModelParams.Validate  = 100/688;  %what fraction of data will used for diffusion map
    [tResults, tParams, tFuncInd, tRes_BEF, Data] = ...
        NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3L1Spec', Calculate, MainModelParams);
    MetaDatas{iTS} = Data;
end
%plot figure 5
fg = figure;
iFgNr = DiffMaps_Adds_Plots('Fig5', fg.Number, struct('MetaDatas',{MetaDatas}));
AddLetters2Plots