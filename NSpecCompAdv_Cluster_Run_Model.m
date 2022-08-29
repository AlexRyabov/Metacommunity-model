%% Modeling spatial resource competiton
%% Alexey Ryabov 2013-2019
% Dorothe Hodapp 2014,
% Barbara Bauer 2019
% - introduced conversion factor (c_i) in growth limitation term
% - changed parametrization so that c_i values correspond to molar N and P
% proportions
%%
%% Parameters
% the script is divided into blocks which
% 1. initialize the model
% 2. run the model
% 3. process the results
% The user should setup only write code for the initialization
% and processing blocks.
% |ScriptID|  is the prefix of the script files describing the model setup and process of the results
% these scripts should be called ScriptID_Setup.m and ScriptID_ProcessResults.m
% e.g. for ScriptID = 'EcolReply', the names should be EcolReply_Setup.m and EcolReply_ProcessResults.m

% in the initialization file you need to define a cell array
% |Params| which includes records
% with parameters for a set of simulations. The results will be saved also in
% cell arrays. So each elment of the resulting cell array will match a
% corresponding element in the |Params| array

% |SimulatID| is the simulation ID in the scripts, so in each setup and
% processing block you can have a switch statement to setup slightly different
% variants of the model
%
% |Calculate| defines wheather the model is calculated ( = 1) or should be loaded from a saved file and plot results (=0)
% results will be saved in a subfolder ./Data/ScriptID

% |MainModelParams| not a necessary paremeter. It is a structure of the main model paramters, which will
% be used in the setup and process script. The structure and the way to use
% it are completly defined by the user. The pupose of this aurgument to
% send some dznamic data to the script which setups the model. The user
% should of course define the way to use this date in the setup script.

%output paramters 
%tResults table with modelling results
%tParams  array of records with model and species paramters which were assigned in each
%simulation
%tFuncInd  table with functional diversity and diversity indexes 
%tRes_BEF some parameters
%Data    a place holder for user defined data structure


%% Example
% NSpecCompAdv_Cluster_Run_Model('EnvBiom', 'Stnr_dyn', 0)  - scripts EnvBiom_EcolReply.m and EnvBiom_EcolReply.m
% will be used to setup the model and process the results, where the block 'Stnr_dyn' will be used.
% Simulation will not run, but saved data will be loaded, so the simulation should be done before.

% To run on a cluster (uni Oldenburg):
% sched = findResource('scheduler', 'Configuration', 'HERO');
% ReplytoEcol_l2 = batch(sched,'NSpecCompAdv_Cluster_Run_Model(''EcolReply'',''ReplytoEcol_l'', 1)','matlabpool', 8, 'FileDependencies',{'NSpecCompAdv_Cluster_Run_Model', 'EcolReply_Setup', 'run_MetaCom.m','SpatLocalChange.m','SpatBoundary.m','Laplace2D_2x2.m','Model_Params.m', 'MC_GetHalfSaturation.m', 'randresshuffle.m', 'resourceshuffle.m'})
%               job = batch(sched,'NSpecCompAdv_Cluster_Run_Model(''EcolReply'',''Diff001'',1)','matlabpool', 8, 'FileDependencies',{'NSpecCompAdv_Cluster_Run_Model', 'SetupAlex', 'run_MetaCom.m','SpatLocalChange.m','SpatBoundary.m','Laplace2D_2x2.m','Model_Params.m', 'MC_GetHalfSaturation.m'});

function [tResults, tParams, tFuncInd, tRes_BEF, Data] = ...
    NSpecCompAdv_Cluster_Run_Model(ScriptID, SimulatID, Calculate, MainModelParams)


close all


%% default paramters for all simulations
NumReplica = 1; % set the number of replicas of one simulation
GenerateMonocultures = 0;   %1 for generate monocultures
ReshuffleStep = 0;
Params = [];

ResDistr = 'g';     %Default resource distribution g - gradient, l - random locations, s - random supplies

ParamDefault.GetFinalSolultion = false; %set true to get the final solution
%for which the right hand side of the model equals zero
%works terribly slow. change the
%method?
ParamDefault.Solver = 'ode45'; %define the ode solver   {'Euler', 'ode45', 'cvode'}
%ParamDefault.Solver = 'cvode'; %define the ode solver   {'Euler', 'ode45', 'cvode'}

ParamDefault.BoundaryCond = 0;      % default boundary condition. 0 -- zero flux, 1 -- periodical boundary
ParamDefault.comptype = 'C';        % P for positive correlation between K and c, N for negative
ParamDefault.RStar = [];            % set default for R*; if empty then the TraitDist will be used to generate R*
ParamDefault.TraitDist = 'eT';      % default distribution of R* in the trait range. eT - equidistantly; rT - randomly;
ParamDefault.MaxGrowthRate = [];    % set default for maximum growth rate, if empty all growth rates = 1;
ParamDefault.m = [];                % set default for mortality, if empty mortality  = 0.25;

ParamDefault.TreatmentPlot = 0;       % 0 for a control plot, >0 for others;

%default grid setup
ParamDefault.DispRes = 0;                 %Diffusion of resources

ParamDefault.S = [];                % default resource supplies over the grid. If empty they will be generated in test_RCSetup
ParamDefault.ReSup = 'eS';          % default distribution for supply points. eS - equidistantly; rS - randomly; rLogNorm - lognormal, see the rest in Model_Params.m
ParamDefault.R2toR1SupplRatio = 1;  % set default for resource supply ratio
ParamDefault.SSupplActivity = @(t) 1;  %Dilution rate activity
ParamDefault.SSupplResourceShuffle = @(t, S, Params) resourceshuffle(S); %function for random shaffling of the reosources

p.Env.ReshuffleTimes = 0;           %a vector or constant of resource reshuffle times.

FertilizationTime = 0;                           %  Fertilization time, default = 0, no fertilization
% if positive you should  also define
% ParamDefault.FertTime,
% ParamDefault.Fertilized, ParamDefault.SAddition
ParamDefault.FertTime = FertilizationTime;       % Fertilization time for each experiment
ParamDefault.Fertilized = false;                 % a flag, which becomes true after the fertilization
ParamDefault.SAddition = [];                     % Resource addition to the supply

ParamDefault.RAddition = [];                      % max of random addition to the local resources

TimeSteps2Save = 40;                             % the number of time steps which will be saved

ParamDefault.PlotSolutions = false;              %show solutions in run time. Works only when the system is not preturbed

%what shell we save ~~~~~~~~~~~~~~~~~~~~~~~~~
SaveSpeciesSpaceTimeDynamics = 1;    %save species dynamics in time in each point of the grid
SaveSpeciesTimeDynamics = 1;         %save species dynamics in time averaged across the grid

%%initialize output results
Data = [];  %user defined data structure.



%Call the script which sets up the model-------------
eval([ScriptID '_RCSetup']);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END OF PARAMETER Section ~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Modellling spatial dynamics~~~~~~~~~~~~~~~~

%profile on -history
if Calculate
    Replicas = [];
    Results = [];
    for iParam = 1:length(Params)
        %%use parfor loop to run on many cores
        %    parfor iParam = 1:10%length(Params)
        iParam
        %initialize the model
        ModParams = Params(iParam);
        [dx, Lx, Ly, tmin, tmax, ...
            res, D, S, R, ResDyn, d, r, m, R_eq, c, N, ModParams] = ...
            Model_Params (ModParams, iParam);
        Params(iParam).S = ModParams.S;
        Params(iParam).consumprate = c;
        % ModParams = NewParam;
        %setup TSpan
        if ~exist('TSpan', 'var')
            TSpan = linspace(tmin, tmax, TimeSteps2Save);
            TSpan = [linspace(TSpan(1), TSpan(2), 20), TSpan(3:end)];
        end
        %the start of fertilization
        
        %% run the model
        tic
        FertStart_ind = find(TSpan < FertilizationTime, 1, 'last');
        [TSpan, N_Dyn, R_Dyn] = run_MetaCom...
            (N, c, R_eq, r, m, d, R, S, D, TSpan, dx, Lx, Ly, ModParams);
        toc
        %~~~~the format of the results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % N_Dyn(x_bins, y_bins, nSpecies_number, nTimeSteps) -- species
        % dynamics in each cell
        % R_Dyn(x_bins, y_bins, nResources, nTimeSteps) -- resource
        % dyncamis in each cell
        % example
        % N_Dyn(x, y, Sp_number, :) the dynamics of species Sp_number
        % at the point x, y
        % N_Dyn(:, :, Sp_number, end) the final distribution of species
        % Sp_number accross the grid
        
        %%the structure, which will be saved for further analysis
        Results(iParam).Params =  ModParams;   %%paramters
        Results(iParam).TSpan  =  TSpan    ;        %%time span
        if FertStart_ind > 1
            Results(iParam).b4Fert =  N_Dyn(:, :, :, FertStart_ind); %biomass distribution just before fertilization
        end
        Results(iParam).Nfin   =  N_Dyn(:, :, :, end);  %final biomass  distribution
        Results(iParam).Rfin   =  R_Dyn(:, :, :, end);  %final resource distribution
        Results(iParam).R_eq   =  R_eq;                 %R* values of the species
        Results(iParam).c      =  c;
        Results(iParam).S      =  S     ;              %supply points
        Results(iParam).Ninit   =  N_Dyn(:, :, :, 1);    %initial biomass distribution
        
        if SaveSpeciesSpaceTimeDynamics  %%if SaveSpeciesDynamics == 1 then save the complete spatial dynamics
            Results(iParam).N_SpaceTimeDyn   =  N_Dyn(:, :, :, :);
            Results(iParam).R_SpaceTimeDyn   =  R_Dyn(:, :, :, :);
        end
        
        if SaveSpeciesTimeDynamics
            %calculate the dynamics of the total biomass and total resource
            Results(iParam).GlobalBiomass_dyn = squeeze(sum(sum(N_Dyn, 1), 2));
            Results(iParam).GlobalRes_dyn     = squeeze(sum(sum(R_Dyn, 1), 2));
        end
        %%
        
        
        if mod(iParam, 100) == 0
            [pathstr, name, ext] = fileparts(DataFileName);
            mkdir(pathstr);
            %Convert data into tables
            tResults = struct2table(Results, 'AsArray', true);
            tParams = struct2table(Params, 'AsArray', true);
            save(DataFileName, 'Params', 'tParams', 'tResults');  %save results and paramters
            
            %% some sample plots which can be done
            figure(1)
            subplot(2, 3, 1)
            plot(TSpan(2:end),  Results(iParam).GlobalBiomass_dyn(:, 2:end)', '-', 'Marker', '.'); % species biomass over time
             yl = ylim;
            f_Set_scales('y');
            ylim([1e-5, yl(2)]);
            xlabel('Time'); ylabel('biomass');
            subplot(2, 3, 2)
            %scatter(ModParams.MaxGrowthRate, Results(iParam).GlobalBiomass_dyn(:, end))
            %f_Set_scales({'x', 'logy'})
            %xlabel('\mu_{max}'); ylabel('biomass');
            plot(TSpan(2:end),  squeeze(R_Dyn(5, 5, :, 2:end)) , '-', 'Marker', '.'); % species biomass over time
            xlabel('Time'); ylabel('Resource');
            f_Set_scales('logy');
            subplot(2, 3, 3)
            scatter( ModParams.RStar(:, 1), ModParams.MaxGrowthRate)
            xlabel('R_1*'); ylabel('\mu_{max}');
            subplot(2, 3, 4)
            cla
            s = plot( ModParams.RStar(:, 1), ModParams.RStar(:, 2), '.')
            hold on
            s = scatter( ModParams.RStar(:, 1), ModParams.RStar(:, 2),  Results(iParam).GlobalBiomass_dyn(:, end)', 'filled')
            hold off
            xlabel('R_1*'); ylabel('R_2*');
            
            subplot(2, 3, 5)
            plot(TSpan(2:end),  squeeze(N_Dyn(5, 5, :, 2:end)) , '-', 'Marker', '.'); % species biomass over time
            yl = ylim;
            f_Set_scales('y');
            ylim([1e-5, yl(2)]);
            subplot(2, 3, 6)
            semilogy(sort(Results(iParam).GlobalBiomass_dyn(:, end), 'descend'), '.')
            drawnow
        end
    end
    %% profresults = profile('info')
    
    [pathstr, name, ext] = fileparts(DataFileName);
    mkdir(pathstr);
    %Convert data into tables
    tResults = struct2table(Results, 'AsArray', true);
    tParams = struct2table(Params, 'AsArray', true);
    save(DataFileName, 'Params', 'tParams', 'tResults');  %save results and paramters
    %save(DataFileName, 'Params', 'Results');
    %return;
else
    %ParamsOld = Params;
    load(DataFileName, '-mat');    %load data
    %     %Params = ParamsOld;
         Params = Params(1:size(tResults, 1));
         tParams = tParams(1:size(tResults, 1), :);
end



%% examples 
%%how to select all control plots
%ind_cont_plot = tParams.TreatmentPlot == 0;
%%take results for the control plots
%tResults(ind_cont_plot, :)
%%how to select all monocultures for experimental plost
%ind_cont_plot = tParams.TreatmentPlot > 0 & tParams.Monoculture == 1;


%% plot biomass dynamics for each experiment
rows = 4;
cols = 4;
%if rows < 7 && cols < 7
Fg01 = figure(1);
set(Fg01, 'Position', [20, 20, 1200, 1200]);
sploti = 1;
for ip = unique(round(linspace(1, size(tResults, 1), 16)))
    if sploti > 16
        break
    end
    subplot(rows, cols, sploti);
    Biom = tResults.GlobalBiomass_dyn{ip};
    plot(tResults.TSpan(ip, :), Biom);
    f_Lbls('Time', 'biomass');
    text()
    sploti = sploti + 1;
    %ylim([1e-2, 1e5])
end


%% Find indexes of functioning (biomass, biodiversity etc for each run)
for ri = 1:length(Params)
    if ri == 1
        FuncInd(length(Params)) = GetBiomassIndexes(tResults.Nfin{ri}); %% initialize the struct array
        FuncInd(1) = FuncInd(length(Params));
    else
        FuncInd(ri) = GetBiomassIndexes(tResults.Nfin{ri}); %
    end
end
tFuncInd = struct2table(FuncInd, 'AsArray', true);


%find average BEF indexes for different configurations
if GenerateMonocultures
    Configurations = unique(tParams.ConfigurationID);
    for ci = 1:length(Configurations)
        ind = tParams.ConfigurationID == Configurations(ci);
        From_i = find(ind > 0, 1, 'first');
        BEF = GetSelCompl(tFuncInd(ind, :), tParams(ind, :));
        BEF.ResRange = tParams.ResHigh(From_i) - tParams.ResLow(From_i);
        RStar = cell2mat(tParams.RStar(From_i));
        BEF.RstRange = abs(RStar(1) - RStar(2));
        BEF.Diff    = tParams.Diff(From_i);
        if ci == 1
            Res_BEF(length(Configurations)) = BEF;
            Res_BEF(1)  = BEF;
        else
            Res_BEF(ci) = BEF;
        end
        
    end
    
    tRes_BEF = struct2table(Res_BEF);
else
    tRes_BEF = [];
end

%~~~~~~~~~~~~~~~~~~~~~~~~END Of SPATIAL DYNAMICS~~~~~~~~~~~~~~~~~~~~~~~~~~
switch ScriptID
    case 'EcolReply'
        EcolReply_ProcessResults;
    case 'NutNet'
        
    otherwise
        eval([ScriptID '_RCPlot']);
end

%save(DataFileName, 'Params', 'tParams', 'tResults');  %save results and
%paramters again after calculating all indexes. In general not a good
%idea.


end


function FuncInd = GetBiomassIndexes(Nxyi)
Richness_border = 1e-6 * max(Nxyi(:));                                        %take the maximal local abuandance
SpNumber = size(Nxyi, 3);                                                     %total number of species in the model

FuncInd.Cells = size(Nxyi, 1) * size(Nxyi, 2);
FuncInd.BiomassTotal_i = squeeze(sum(sum(Nxyi, 1), 2));                     %total biomass of every species
FuncInd.BiomassMean_i  = FuncInd.BiomassTotal_i/FuncInd.Cells;              %mean biomass of every species
FuncInd.BiomassTotal   = sum(FuncInd.BiomassTotal_i);                       %total biomass of all species

FuncInd.Simpson_glob = sum((FuncInd.BiomassTotal_i/FuncInd.BiomassTotal).^2);  %Simpson index
FuncInd.Simpson_ESN_glob = 1/FuncInd.Simpson_glob;                             %Simpson ESN

%% Local diversity Simpson = sum_i(p_xyi.^2)
SpecPresAnbs_xyi = Nxyi > Richness_border;                                     %presence/anabsence data
FuncInd.Richness_xy = sum(SpecPresAnbs_xyi, 3);                                %local richness
FuncInd.Richness_glob = sum(sum(sum(SpecPresAnbs_xyi, 2), 1) > 0);               %global richness

Nxy_total = sum(Nxyi(:, :, :), 3);                                             %total biomass in every cell
FuncInd.Biomass_xy = Nxy_total;

EmptyCells = Nxy_total < Richness_border;                                       %Exclude empty cells
Nxy_total(EmptyCells) = NaN;

Pxyi = Nxyi(:, :, :)./repmat(Nxy_total, 1, 1, SpNumber);
FuncInd.Simpson_loc_xy = sum(Pxyi.^2, 3);                                      %local Simpson index
FuncInd.Simpson_ESN_loc_xy = 1./FuncInd.Simpson_loc_xy;
FuncInd.Simpson_loc_avg =     nanmean(FuncInd.Simpson_loc_xy(:));              %local Simpson index
FuncInd.Simpson_ESN_loc_avg = nanmean(FuncInd.Simpson_ESN_loc_xy(:));

ExtinctSpec_xyi = Pxyi < 1e-15;                                             %remove extinct species
Pxyi(ExtinctSpec_xyi) = NaN;
FuncInd.Shannon_loc_xy = nansum(-Pxyi .* log(Pxyi), 3);                     %local Shannon index
%FuncInd.Shannon_loc_xy(FuncInd.Shannon_loc_xy == 0) = NaN;                 %remove empty cells
FuncInd.Shannon_ESN_loc_xy = exp(FuncInd.Shannon_loc_xy);                   %local Shannon ESN
FuncInd.Shannon_loc_avg = nanmean(FuncInd.Shannon_ESN_loc_xy(:));           %average local Shannon index
FuncInd.Shannon_ESN_loc_avg = nanmean(FuncInd.Shannon_ESN_loc_xy(:));          %average local Shannon ESN

FuncInd.Pielou_even_loc_xy = FuncInd.Shannon_loc_xy./log(FuncInd.Richness_xy);  %local eveness

end


function BEF = GetSelCompl(tFuncInd, tParams)
Cells = tFuncInd.Cells(1);
ind_mono = tParams.Monoculture == 1;
ind_mix = ~ind_mono;
%%observed yeld for yeach species
Yoi = squeeze(cell2mat(( tFuncInd.TotalBiomassSp(ind_mix))))/Cells;  %%replica x species
%species ID in the mixture
Sp_mix = cell2mat(tParams.specID(ind_mix));                    %%replica x species
%monoculture yelds
Mono_i = cell2mat(tFuncInd.TotalBiomassSp(ind_mono))/Cells;
N_in_mixture = size(Sp_mix, 2);

Compl_eff = zeros(size(Sp_mix, 1), 1);
Sel_eff   = Compl_eff;
Delta_Y   = Compl_eff;

for ri = 1:size(Sp_mix, 1)
    %monoculture yelds of the species in mixture
    Sp_mono_mix = Mono_i(Sp_mix(ri, :))';
    %relative observed yeld
    RYoi = Yoi(ri, :)./Sp_mono_mix;
    %RYei(end+1)= Yei(i)/sum(Yei); %!! Change instead of M/sum(M) we use 1/n
    RYei = 1/N_in_mixture; %!! Change instead of M/sum(M) we use 1/n
    deltaRY = RYoi - RYei;
    mean_deltaRY = mean(deltaRY);
    Compl_eff(ri) = N_in_mixture * mean_deltaRY * mean(Sp_mono_mix);
    covmat =  mean( (deltaRY - mean(deltaRY)).*(Sp_mono_mix - mean(Sp_mono_mix)));
    Sel_eff(ri) = N_in_mixture * covmat;
    Delta_Y(ri) = sum(Yoi(ri, :)) - mean(Sp_mono_mix);
end
BEF.Compl_eff = mean(Compl_eff);
BEF.Sel_eff   = mean(Sel_eff);
BEF.Delta_Y   = mean(Delta_Y);
BEF.SpNumbSimpson_glob = mean(tFuncInd.SpNumbSimpson_glob(ind_mix));
BEF.SpNumbSimpson_loc  = mean(tFuncInd.SpNumbSimpson_loc(ind_mix));
end


function Rst_proxy = GetRstProxies(Nfin,Nfin_ctrl,Nb4Fert,Nb4Fert_ctrl,sp)
for i= 1:sp
    Rst_proxy = size(sp,1);
    Rst_proxy(i) = log(sum(sum(Nfin(:,:,i)))/sum(sum(Nb4Fert(:,:,i))))-log(sum(sum(Nfin_ctrl(:,:,i)))/sum(sum(Nb4Fert_ctrl(:,:,i))));
end
end
