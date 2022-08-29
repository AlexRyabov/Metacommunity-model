%% Modeling spatial resource competiton
%% Alexey Ryabov 2013-2017
%% This script makes setup of modelling of effects of mean and variance in resource distribution
%% on biodiversity

%% define your root folder for data
if Calculate == 1  %% calculate on the cluster
    RootDataFolder = 'Data/';
else
    RootDataFolder = 'Data/';
end

%% define parameters of your simulations
%switch SimulatID
%Random maximal growth rate, random R* values with a trade-off in
%resource requirements. The half-saturation constants will be
%calculated from the R*, mu_max and m as H = R* m/(mu_max - m); 10
%replicas with different R* and growth rates in the same range
%Sigma=Mu^2     --  variance of both resources will be proportional to
%mean squared
%Sigma=Mu^2_dyn50 -- the same, but resource supplies will be redistributed over the
%grid with the period given by ReshuffleStep
%    case {'Sigma=Mu^2', 'Sigma=Mu^2_dyn50', 'Sigma=Mu^2_dyn05'}
%Time
if contains(SimulatID,'Grid') 
    tmax = 5000;
    TimeSteps2Save = 1000;
    TSpan = [0:TimeSteps2Save:tmax];
elseif contains(SimulatID,'Single')
    tmax = 5000;
    TimeSteps2Save = 1000;
    TSpan = [0:TimeSteps2Save:tmax];
else
    tmax = 10000;
    TimeSteps2Save = 1000;
    TSpan = [0:TimeSteps2Save:tmax];
end


%Model
GenerateMonocultures = 0;
NumReplica = 1;

%Species
sp = 20; %number of species
GrowthRateRange = [0.3,2];
GlOpprtTradeOff = 0; %% trade-off between Rst and growth rate (1) or random trait combination (0)
Diff = [0.1];        %%Set dispersal for species

%Resources
%the ranges of resource supply points
switch SimulatID
    case {'Sigma=Mu^1.5', 'GridSigma=Mu^1.5', 'Sigma=Mu^1.5_dyn50', 'Sigma=Mu^1.5_dyn05'}
        Exponent = 1.5;
    case {'Sigma=Mu^2', 'GridSigma=Mu^2', 'Sigma=Mu^2_dyn50', 'Sigma=Mu^2_dyn05'}
        Exponent = 2;
    case {'Sigma=Mu^2.5', 'GridSigma=Mu^2.5', 'Sigma=Mu^2.5_dyn50', 'Sigma=Mu^2.5_dyn05'}
        Exponent = 2.5;
    case {'Sigma=Mu^3', 'GridSigma=Mu^3', 'Sigma=Mu^3_dyn50', 'Sigma=Mu^3_dyn05'}
        Exponent = 3;
    case 'Single'
        Exponent = 2;
    otherwise
        error('wrong ID');
end
if contains(SimulatID,'Grid') 
    ResDistr = 's';   %%s for random supplies  ('g' -- gradient; 'l' -- random locations)
    SMean = logspace(log10(5), 3.5, 30)';  %mean value
    ReSup = 'rLogNorm';   
else  %Single simulation
    Diff = [0.05];  
    SMean = 500;  %mean value
    ReSup = 'rLogNorm';   
    ResDistr = 's';   %%s for random supplies  ('g' -- gradient; 'l' -- random locations)
    
end
SVar =  0.8 * sqrt(10.^(Exponent * log10(SMean)));  %variance
ResourceRangesLowHigh = ...
    [SMean - SVar/2, SMean + SVar/2];
TraitRangesLowHigh = [1,10;];   %the minimal and maximal R1* and R2* values
FertilizationTime = 0;     %time when new nutrients will be added


%Names
DataFileName = [SimulatID '/' num2str(sp) '_'  ResDistr num2str(tmax) '.dat'];
ParamName_field = 'ResR';            %% the name of the varying variable in the script
ParamName_plot   = 'Resource range'; %% the name of this variable in figures

switch SimulatID
    case    {'Sigma=Mu^1.5_dyn50', 'Sigma=Mu^2_dyn50', 'Sigma=Mu^2.5_dyn50', 'Sigma=Mu^3_dyn50'}
        ReshuffleStep = 50;  %reshuffle resource supply points every ReshuffleStep steps
    case {'Sigma=Mu^1.5_dyn05', 'Sigma=Mu^2_dyn05', 'Sigma=Mu^2.5_dyn05', 'Sigma=Mu^3_dyn05'}
        ReshuffleStep = 05;  %reshuffle resource supply points every ReshuffleStep steps
end

if contains(SimulatID,'Grid') ||  contains(SimulatID,'Single')  % we symulate a grid of supply clouds
    %% assing paramters of the simulations
    ConfigurationID = 1;                             %unique ID for each configuration
    for j = 1:size(ResourceRangesLowHigh, 1)
        for k = 1:size(ResourceRangesLowHigh, 1)
            for i = 1:size(TraitRangesLowHigh, 1)
                p = ParamDefault;                        %initialize with default values
                if ReshuffleStep > 0
                    p.Env.ReshuffleTimes = ReshuffleStep:ReshuffleStep:tmax;
                else
                    p.Env.ReshuffleTimes = 0;
                end
                
                p.ReSup = ReSup;                    %lognormal dsitribution of resource supplies
                p.TreatmentPlot = 0;                     %0 for a control plot
                p.Lx = 20;                               %x dimension
                p.Ly = 25;                               %y dimension
                p.BoundaryCond = 0;                      %periodical boundary conditions
                p.ResLow(1) = ResourceRangesLowHigh(j, 1);  %the minimal values of the resource supplies
                p.ResLow(2) = ResourceRangesLowHigh(k, 1);  %the minimal values of the resource supplies
                p.ResHigh(1) = ResourceRangesLowHigh(j, 2); %the maximal values of the resource supplies
                p.ResHigh(2) = ResourceRangesLowHigh(k, 2); %the maximal values of the resource supplies
                
                p.TraitLow = TraitRangesLowHigh(i, 1);   %the minimal R* value
                p.TraitHigh = TraitRangesLowHigh(i, 2);  %the maximal R* value
                p.ResDistr = ResDistr;                   %the resource distribution
                p.Diff     = Diff;                       %the species dispersal rate
                p.tmax = tmax;                           %the simulation time
                p.tfert = FertilizationTime;
                p.ConfigurationID = ConfigurationID;
                ConfigurationID = ConfigurationID + 1;

                %%generate replicas
                for k=1:NumReplica
                    
                    %Randomly distributes R* values with a trade-off in resource requirements
                    RStar(:, 1) = p.TraitLow + rand(sp, 1) * (p.TraitHigh - p.TraitLow);        %R* values
                    RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R* values
                    
                    MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) * (GrowthRateRange(2) - GrowthRateRange(1));
                    
                    p.Monoculture = 0;
                    p.Replica     = k;
                    %
                    
                    p.specID        = 1:sp;        %here we take all species, but you can take you can define e.g. p.specID = [1, 5, 6]
                    %to select which species will be used
                    p.MaxGrowthRate = MaxGrowthRate(p.specID);
                    p.RStar         = RStar(p.specID, :);
                    p.RStar_All     = RStar; %entire ensemble
                    p.sp            = length(p.specID);
                    if isempty(Params)
                        Params = p;
                    else
                        Params(end + 1) = p;
                    end
                end
            end
        end
    end
else
    %% assing paramters of the simulations
    ConfigurationID = 1;                             %unique ID for each configuration
    for j = 1:size(ResourceRangesLowHigh, 1)
        for i = 1:size(TraitRangesLowHigh, 1)
            p = ParamDefault;                        %initialize with default values
            p.ReSup = ReSup ;                    %lognormal dsitribution of resource supplies
            p.TreatmentPlot = 0;                     %0 for a control plot
            p.Lx = 20;                               %x dimension
            p.Ly = 25;                               %y dimension
            p.BoundaryCond = 0;                      %periodical boundary conditions
            p.ResLow(1:2) = ResourceRangesLowHigh(j, 1);  %the minimal values of the resource supplies
            p.ResHigh(1:2) = ResourceRangesLowHigh(j, 2); %the maximal values of the resource supplies
            
            p.TraitLow = TraitRangesLowHigh(i, 1);   %the minimal R* value
            p.TraitHigh = TraitRangesLowHigh(i, 2);  %the maximal R* value
            p.ResDistr = ResDistr;                   %the resource distribution
            p.Diff     = Diff;                       %the species dispersal rate
            p.tmax = tmax;                           %the simulation time
            p.tfert = FertilizationTime;
            p.ConfigurationID = ConfigurationID;
            ConfigurationID = ConfigurationID + 1;
            
            %%generate replicas
            for k=1:NumReplica
                
                %Randomly distributes R* values with a trade-off in resource requirements
                RStar(:, 1) = p.TraitLow + rand(sp, 1) * (p.TraitHigh - p.TraitLow);        %R* values
                RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R* values
                
                MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) * (GrowthRateRange(2) - GrowthRateRange(1));
                
                p.Monoculture = 0;
                p.Replica     = k;
                %
                
                p.specID        = 1:sp;        %here we take all species, but you can take you can define e.g. p.specID = [1, 5, 6]
                %to select which species will be used
                p.MaxGrowthRate = MaxGrowthRate(p.specID);
                p.RStar         = RStar(p.specID, :);
                p.RStar_All     = RStar; %entire ensemble
                p.sp            = length(p.specID);
                if isempty(Params)
                    Params = p;
                else
                    Params(end + 1) = p;
                end
            end
        end
    end
end

DataFileName = [RootDataFolder DataFileName ];