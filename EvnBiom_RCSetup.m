%% Modeling spatial resource competiton
%% Alexey Ryabov 2013-2017
%% This script makes setup for different calcualtions to evaluate the relationships
%% between the final biomass and eveness of species distribution

%% define your root folder for data
if Calculate == 1  %% calculate on the cluster
    RootDataFolder = 'Data/';
else
    RootDataFolder = 'Data/';
end

%% define paramters of your simulations
switch SimulatID
    %Random maximal growth rate, random R* values with a trade-off in
    %resource requirements. The half-saturation constants will be
    %calculated from the R*, mu_max and m as H = R* m/(mu_max - m); 10
    %replicas with different R* and growth rates in the same range
    %Stnr_dyn   -- stationary conditions
    %MixSrs_dyn -- the resource supplies will be redistributed over the
    %grid with the period given by ReshuffleStep 
    case {'Stnr_dyn', 'MixSrs_dyn', 'MixSrs_dyn05'}
        %Time
        tmax = 5000;
        TimeSteps2Save = 100;
        
        %Model
        GenerateMonocultures = 0; 
        NumReplica = 3;
        
        %Species
        sp = 20; %number of species
        GrowthRateRange = [0.3,2];
        GlOpprtTradeOff = 0; %% trade-off between Rst and growth rate (1) or random trait combination (0)
        Diff = [0.1];        %%Set dispersal for species

        %Resources
        %the ranges of resource supply points
        ResDistr = 's';   %%s for random supplies  ('g' -- gradient; 'l' -- random locations)
        switch SimulatID
            case 'MixSrs_dyn05'
                dS = [5:10:35]'; dS = dS(end:-1:1);  %the range of resource variability
            otherwise
                dS = [5:5:30]'; dS = dS(end:-1:1);  %the range of resource variability
        end
        ResourceRangesLowHigh = ...
            [20.5 - dS/2, 20.5 + dS/2];
        TraitRangesLowHigh = ...   %the minimal and maximal R1* and R2* values
            [1,10; 3,8];
        FertilizationTime = 0;     %time when new nutrients will be added
        
        
        %Names
        DataFileName = [SimulatID '/' num2str(sp) '_'  ResDistr num2str(tmax) '.dat'];
        ParamName_field = 'ResR';            %% the name of the varying variable in the script
        ParamName_plot   = 'Resource range'; %% the name of this variable in figures

        switch SimulatID
            case    'MixSrs_dyn'
                ReshuffleStep = 50;  %reshuffle resource supply points every ReshuffleStep steps
            case 'MixSrs_dyn05'
                ReshuffleStep = 05;  %reshuffle resource supply points every ReshuffleStep steps
        end
        %% assing paramters of the simulations
        ConfigurationID = 1;                             %unique ID for each configuration
        for j = 1:size(ResourceRangesLowHigh, 1)
            for i = 1:length(TraitRangesLowHigh)
                p = ParamDefault;                        %initialize with default values
                p.TreatmentPlot = 0;                     %0 for a control plot
                p.Lx = 20;                               %x dimension
                p.Ly = 25;                               %y dimension
                p.BoundaryCond = 0;                      %periodical boundary conditions
                p.ResLow = ResourceRangesLowHigh(j, 1);  %the minimal values of the resource supplies
                p.ResHigh = ResourceRangesLowHigh(j, 2); %the maximal values of the resource supplies
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
                    
                    %%bad competitor -- good competitor (no trade-off)
                    %    RStar = [linspace(p.TraitLow, p.TraitHigh, sp)', ...        %R* values
                    %         linspace(p.TraitLow, p.TraitHigh, sp)'];

                    %%Trade_off in R*   equidistantly distributed R* values
                    %%on a line
                    %    RStar = [linspace(p.TraitLow, p.TraitHigh, sp)', ...        %R* values
                    %         linspace(p.TraitHigh, p.TraitLow, sp)'];
                    
                    %Randomly distributes R* values with a trade-off in resource requirements
                    RStar(:, 1) = p.TraitLow + rand(sp, 1) * (p.TraitHigh - p.TraitLow);        %R* values
                    RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R* values
                    
                    GlOpprtTradeOff = 0;
                    if GlOpprtTradeOff
                        %gleaner-ooportunist trade-off
                        MaxGrowthRate = GrowthRateRange(1) + (GrowthRateRange(2)...
                            - GrowthRateRange(1)) * (sum(RStar,2)-2*p.TraitLow)/(2*(p.TraitHigh-p.TraitLow));
                    else
                        % random growth rate is a random value 
                        MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) * (GrowthRateRange(2) - GrowthRateRange(1));
                    end
                    
                    %                         if GenerateMonocultures
                    %                             %%generate monoculture experiments
                    %                             for si = 1:sp
                    %                                 p.specID = si;   %contains IDs of all species
                    %                                 p.Monoculture = 1;
                    %                                 p.Replica = 1;   %Number of replica
                    %                                 p.sp   = 1;      %number of species
                    %                                 p.MaxGrowthRate = MaxGrowthRate(si);
                    %                                 p.RStar = RStar(si, :);
                    %                                 if isempty(Params)
                    %                                     Params = p;
                    %                                 else
                    %                                     Params(end + 1) = p;
                    %                                 end
                    %                             end
                    %                         end
                    %
                    p.Monoculture = 0;
                    p.Replica     = k;
                            %
                
                    p.specID        = 1:sp;        %here we take all species, but you can take you can define e.g. p.specID = [1, 5, 6]
                                                   %to select which species will be used
                    p.MaxGrowthRate = MaxGrowthRate(p.specID);
                    p.RStar         = RStar(p.specID, :);
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