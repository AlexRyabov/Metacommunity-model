%% Modeling spatial resource competiton
% Alexey Ryabov 2013-2022
% This script makes setup for the diffusion map analysis
% The goal is to generate the dynamics of different species groups for
% 1. a wide and narrow trait range
% 2. two and three clusters of species
% 3. two clusters of species including B and C competitors
% 4. do it with small, moderate and large ReshuffleStep

%before calling this script you should define
%MainModelParams.ReshuffleStep (e.g. 50)
%MainModelParams.NumberOfSpecies  (e.g. 20)


%% define your root folder for data
if Calculate == 1  %% calculate on the cluster
    RootDataFolder = 'Data/';
else
    RootDataFolder = 'Data/';
end

DataFileName = [RootDataFolder 'DiffMap/' ScriptID '_' SimulatID '_' ...
    sprintf('rs=%.f.mat', ...
    MainModelParams.ReshuffleStep)];

p = ParamDefault;                        %initialize parmaeteres with default values

RAddition =[];
% here

%% define paramters of your simulations
%Random maximal growth rate, random R* values with a trade-off in
%resource requirements. The half-saturation constants will be
%calculated from the R*, mu_max and m as H = R* m/(mu_max - m); 10
%replicas with different R* and growth rates in the same range
%Stnr_dyn   -- stationary conditions
%MixSrs_dyn -- the resource supplies will be redistributed over the
%grid with the period given by ReshuffleStep
tmax = 6000;
TimeSteps2Save = 100;

tmax1 = tmax-3000; %skip first 1000 days
TDiff = (randn(1, tmax1/TimeSteps2Save) + 1);
TSpan = sort(cumsum(TDiff));
TSpan = TSpan/max(TSpan)*tmax1;
TSpan = [TSpan + 3000];

ReshuffleStep = MainModelParams.ReshuffleStep;
%Model
GenerateMonocultures = 0;
NumReplica = 5;

%Species
sp = 1000; %maximal number of species
GrowthRateRange = [1,1];
GlOpprtTradeOff = 0; %% trade-off between Rst and growth rate (1) or random trait combination (0)

%Resources
%the ranges of resource supply points
ResDistr = 's';   %%s for random supplies  ('g' -- gradient; 'l' -- random locations)
ResourceRangesLowHigh = ...
    [1, 39];
switch SimulatID
    case 'DeltaRsIncr'
        %increase the trait range of competitor
        TraitRangeDelta = linspace(0.5, 4.5, 20)';
        TraitRangesLowHigh  = [5-TraitRangeDelta, 5 + TraitRangeDelta ];
        SpPerSimulation = 100; %the number of species used in each simulation
    case 'NIncr'    
        SpPerSimulation = round(logspace(1, log10(400), 20)); %species number changes from 10 to 400
        TraitRangesLowHigh = [ 0.5000    9.5000];
    case 'DeltaSIncr'  %chage resource range
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangeDelta = linspace(0.5, 19, 40)';
        ResourceRangesLowHigh  = [20-ResourceRangeDelta, 20 + ResourceRangeDelta ];
        NumReplica = 2;  
    case 'DeltaSIncrStabCentr'  %chage resource range, no reshafling
        ReshuffleStep = 0;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangeDelta = linspace(0.5, 19, 40)';
        ResourceRangesLowHigh  = [20-ResourceRangeDelta, 20 + ResourceRangeDelta ];
        NumReplica = 3;
        tmax = 6000;
        TimeSteps2Save = (tmax-3000)/10;
        TSpan = 3000:TimeSteps2Save:tmax;
        
    case {'DeltaSRandShift', 'DeltaSRandShiftR3', 'DeltaSRandShiftR3Reg', ...
            'DeltaSRandShiftR3L', 'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3D03', ...
            'DeltaSRandShiftR3Gr', 'DeltaSRandShiftR3L1SpecPeriod', 'DeltaSRandShiftR3L1SpecDyn', ...
            'DeltaSRandShiftR1L1SpecGlOp'}  
        %varios random resource range and independ for each res
        %DeltaSRandShiftR3 the same but for 3 resources
        %DeltaSRandShiftR3L1Spec
        %DeltaSRandShiftR3L1SpecPeriod periodical supply rate of resources
        ReshuffleStep = 0;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangeDelta = linspace(0.5, 19, 40)';
        ResourceRangesLowHigh  = [1, 39];
        switch SimulatID
            case {'DeltaSRandShiftR3L', 'DeltaSRandShiftR3L1Spec'}
                NumReplica = 800;
                %NumReplica = 50;
                tmax = 10000;
                TSpan = [1, 6000, tmax]; %used for the paper!!
                %TSpan = [1:30:5000];  %high resolution data 
            case  {'DeltaSRandShiftR3L1SpecPeriod', 'DeltaSRandShiftR1L1SpecGlOp'}
                NumReplica = 500; 
                tmax = ceil(365*5.5);
                TSpan = linspace(tmax-365, tmax, 13);
                TSpan = 0.5*(TSpan(1:end-1) + TSpan(2:end)); %12 months
                TSpan = [1 TSpan];
                GrowthRateRange = [0.5,2];                
                GlOpprtTradeOff = 1;
                %p.SSupplActivity = @(t)((cos(2*pi*t/365)+1)/2).^2*0.8 + 0.2;
                %p.SSupplActivity = @(t)((cos(2*pi*t/365)+1)/2).^2*1 *0.9 + 0.1;
                p.SSupplActivity = @(t)((cos(2*pi*t/365)+1)/2);
                
            case 'DeltaSRandShiftR3L1SpecDyn'
                %p.SSupplResourceShuffle = @(t, S, Params) f_SupplySetup(t, Params);
                p.SSupplResourceShuffle = @(t, S, Params) f_SupplySetupRndWindow(t, Params);
                NumReplica = 100; 
                ReshuffleStep = 50;
                tmax = 2000;
               TSpan = [0, 1000:30:tmax];  %
              %  TSpan = [0, 1000:5:tmax]; %just to plot sample dynamics 
%                 tmax = 5000;
%                 TSpan = [3990:30:tmax];
            otherwise
                NumReplica = 60;
                tmax = 10000;
                TSpan = [1, 6000, tmax];
        end
        
    case {'DeltaSWideR3D03'}  
        %one large grid with different variations of 3 resources 
        %small diffusion
        ReshuffleStep = 0;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangesLowHigh  = [1, 39];
%         NumReplica = 10;
%         tmax = 1000;
%         TSpan = [1, linspace(600, tmax, 10)];
         NumReplica = 100;
         tmax = 10000;
         TSpan = [1, linspace(6000, tmax, 10)];
    case {'DeltaSLarge_R3'}  
        %one large grid with different variations of 3 resources 
        %small diffusion
        ReshuffleStep = 0;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangesLowHigh  = [1, 39];
        NumReplica = 1;
        
        tmax = 10000;
        TSpan = [1, linspace(6000, tmax, 10)];
    case {'DeltaSLarge_R3_Tune'}  
        %one large grid with different variations of 3 resources 
        %small diffusion
        ReshuffleStep = 0;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangesLowHigh  = [1, 39];
        NumReplica = 1;
        tmax = 8000;
        TSpan = 1:100:5000;
        p.PlotSolutions = true;
    case {'DeltaSRandShiftGlOpportResh', 'DeltaSConstGlOpportResh'}  
        %reshafling, 
        %random  res range 
        %gleaner opport
        
        ReshuffleStep = 100;
        sp = 200; %maximal number of species
        TraitRangesLowHigh = [ 0.5000    9.5000];
        ResourceRangeDelta = linspace(0.5, 19, 40)';
        ResourceRangesLowHigh  = [1, 39];
        NumReplica = 50;
        
        %gleraner opportuninst
        GrowthRateRange = [0.5,2];
        GlOpprtTradeOff = 1;
        RAddition = [50, 50];
        %save 10 times at random steps from t_relax to tmax
        tmax = 20000;  
        t_relax = 7000;
        tmax1 = tmax-t_relax; %skip first 6000 days
        TDiff = (randn(1, 30) + 1);
        TSpan = sort(cumsum(TDiff));
        TSpan = TSpan/max(TSpan)*tmax1;
        TSpan = [0, TSpan + t_relax];
    otherwise
        TraitRangesLowHigh = [ 0.5000    9.5000];
        
end
FertilizationTime = 0;     %time when new nutrients will be added

%% assing paramters of the simulations
p.TreatmentPlot = 0;                     %0 for a control plot
switch SimulatID
    case 'DeltaSLarge_R3'
        p.Lx = 100;                               %x dimension
        p.Ly = 101;                               %y dimension
    case {'DeltaSRandShiftR3D03', 'DeltaSWideR3D03'}
        p.Lx = 20;                               %x dimension
        p.Ly = 21;                               %y dimension
    case {'DeltaSLarge_R3_Tune'}
        p.Lx = 20;                               %x dimension
        p.Ly = 21;                               %y dimension
    otherwise
        p.Lx = 10;                               %x dimension
        p.Ly = 12;                               %y dimension
end
switch SimulatID
    case {'DeltaSRandShiftR3D03', 'DeltaSWideR3D03'}
        p.DispRes = 0.3;
end
Diff = [0.1];                            %%Set dispersal for species
p.Diff     = Diff;                       %the species dispersal rate
p.BoundaryCond = 0;                      %periodical boundary conditions
p.ResLow = ResourceRangesLowHigh(1, 1);  %the minimal values of the resource supplies
p.ResHigh = ResourceRangesLowHigh(1, 2); %the maximal values of the resource supplies
p.TraitLow = TraitRangesLowHigh(end, 1);   %the minimal R* value
p.TraitHigh = TraitRangesLowHigh(end, 2);  %the maximal R* value
p.ResDistr = ResDistr;                   %the resource distribution
p.tmax = tmax;                           %the simulation time
p.tfert = FertilizationTime;
if ReshuffleStep > 0
    p.Env.ReshuffleTimes = ReshuffleStep:ReshuffleStep:tmax;
else
    p.Env.ReshuffleTimes = 0;
end
p.RAddition = RAddition;

switch SimulatID
    case {'UniformRst', 'DeltaRsIncr', 'NIncr', 'DeltaSIncr', ...
            'DeltaSIncrStabCentr', 'DeltaSRandShift', 'DeltaSRandShiftGlOpportResh', 'DeltaSConstGlOpportResh'}
        RStar(:, 1) = p.TraitLow + rand(sp, 1) * (p.TraitHigh - p.TraitLow);        %R1* values
        RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R2* values
    case {'DeltaSRandShiftR3', 'DeltaSRandShiftR3L', 'DeltaSRandShiftR3D03', 'DeltaSRandShiftR3Gr', 'DeltaSLarge_R3', 'DeltaSWideR3D03', 'DeltaSLarge_R3_Tune'}
        %RStar(:, :) = rand(sp, 3);        %Random values
        %RStar = RStar./repmat(sum(RStar, 2), 1, 3);  %normilize to have (R1* + R2* + R3*). 
        %Note that we get a higher proportion of intermediate species compared with extreme species
        RStar = f_RandConstSum(sp, 3); %new code, with uniform distribution of R*
        RStar = RStar * (p.TraitHigh - p.TraitLow) + p.TraitLow;
    case {'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3L1SpecPeriod', 'DeltaSRandShiftR3L1SpecDyn'}
        RStar = f_RandConstSum(10*sp, 3); %new code, with uniform distribution of R*
        alpha = 0.5;
        RStar = RStar(RStar(:, 1)< alpha & RStar(:, 2)< alpha & RStar(:, 3)< alpha, :);
        RStar = RStar(1:sp, :)/alpha; %if you do not get enought species generate more with f_RandConstSum above
        RStar = RStar * (p.TraitHigh - p.TraitLow) + p.TraitLow;
        
%         RStar = rand(10*sp, 2)*(p.TraitHigh - p.TraitLow) + p.TraitLow;
%         RStar(:, 3) = 19.5 - RStar(:, 1) - RStar(:, 2);
%         alpha = 9.5;
%         RStar = RStar( RStar(:, 3)< alpha, :);
        %figure; scatter3(RStar(:, 1), RStar(:, 2), RStar(:, 3))
    case {'DeltaSRandShiftR1L1SpecGlOp'}
        RStar(:, 1) = p.TraitLow + rand(sp, 1) * (p.TraitHigh - p.TraitLow);        %R1* values
        RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R2* values
    case {'DeltaSRandShiftR3Reg'}
        steps = 20;
        Rst1 = linspace(0, 1, steps);
        Rst2 = linspace(0, 1, steps);
        [RST1, RST2] = meshgrid(Rst1, Rst2);
        RST3 = 1-RST1-RST2;
        RStar = [RST1(:), RST2(:), RST3(:)];
        RStar = RStar(RST3(:)>=0, :);
        RStar = RStar * (p.TraitHigh - p.TraitLow) + p.TraitLow;
        sp = size(RStar, 1);
    case 'TwoClust'  %remove the middle group
        RSt1 = p.TraitLow  +  rand(sp/2, 1) * 1/3. * (p.TraitHigh - p.TraitLow);
        RSt2 = p.TraitLow  +  2/3. * (p.TraitHigh - p.TraitLow) + rand(sp/2, 1) * (p.TraitHigh - (p.TraitLow  +  2/3. * (p.TraitHigh - p.TraitLow)));
        RStar(:, 1) = [RSt1; RSt2];
        RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R* values
    case 'ThreeClust'  %remove the middle group
        qDelta = (p.TraitHigh - p.TraitLow)/5;
        Qs = linspace(p.TraitLow, p.TraitHigh, 6);
        RSt1 = Qs(1) +  rand(round(sp/3), 1) * (Qs(2) - Qs(1));
        RSt2 = Qs(3) +  rand(round(sp/3), 1) * (Qs(4) - Qs(3));
        RSt3 = Qs(5) +  rand(sp - 2 * round(sp/3), 1) * (Qs(6) - Qs(5));
        RStar(:, 1) = [RSt1; RSt2; RSt3];
        RStar(:, 2) = p.TraitHigh - RStar(:, 1) + p.TraitLow;        %R* values
    otherwise
        error('This ScriptID is not defined in the setup file')
end
RStar = sortrows(RStar);
if GlOpprtTradeOff
    %gleaner-ooportunist trade-off
    switch SimulatID
        case 'DeltaSRandShiftR1L1SpecGlOp'
            MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) *  (GrowthRateRange(2) - GrowthRateRange(1));
            delta = repmat(MaxGrowthRate, 1, size(RStar, 2))/GrowthRateRange(2);
            RStar(:, :) = RStar(:, :) .* (1 * (delta - 1) + 1);
%             ind = RStar(:, 1) >3 | RStar(:, 2)>3;
%             RStar = RStar(ind, :);
%             sp = size(RStar, 1);
        otherwise
            MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) *  (GrowthRateRange(2) - GrowthRateRange(1));
            delta = repmat(MaxGrowthRate, 1, size(RStar, 2))/GrowthRateRange(2);
            RStar(:, :) = RStar(:, :) .* (0.2 * (delta - 1) + 1);
    end
    %scatter(RStar(:, 1), RStar(:, 2), 3, MaxGrowthRate);
    %GrowthRateRange(1) + (GrowthRateRange(2)...
    %    - GrowthRateRange(1)) * (sum(RStar,2)-2*p.TraitLow)/(2*(p.TraitHigh-p.TraitLow));
else
    % random growth rate is a random value
    MaxGrowthRate = GrowthRateRange(1) + rand(sp, 1) * (GrowthRateRange(2) - GrowthRateRange(1));
end

%save paramters of the entire ensemple
p.RStar_All = RStar;
p.MaxGrowthRate_All= MaxGrowthRate;
p.SpecID_All= 1:sp;

ConfigurationID = 1;                             %unique ID for each configuration
switch SimulatID
    case 'DeltaRsIncr'
        for i = 1:size(TraitRangesLowHigh, 1)
            p.ConfigurationID = ConfigurationID;
            ConfigurationID = ConfigurationID + 1;
            for k=1:NumReplica
                p.Monoculture = 0;
                p.Replica     = k;
                %
                ind_RStar_sel = TraitRangesLowHigh(i, 1)  < RStar(:, 1)  & RStar(:, 1) < TraitRangesLowHigh(i, 2);
                p.TraitLow = TraitRangesLowHigh(i, 1);   %the minimal R* value
                p.TraitHigh = TraitRangesLowHigh(i, 2);  %the maximal R* value
                
                RStar_sel = RStar(ind_RStar_sel, :);
                specID_sel = p.SpecID_All(ind_RStar_sel);
                specID_sel = sort(specID_sel(randperm(length(specID_sel), SpPerSimulation)));
                
                p.specID        = specID_sel;        %here we take all species, but you can take you can define e.g. p.specID = [1, 5, 6]
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
    case {'DeltaSIncr', 'DeltaSIncrStabCentr', 'DeltaSRandShift', ...
            'DeltaSRandShiftR3', 'DeltaSRandShiftR3Reg', 'DeltaSRandShiftR3L', ...
            'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3L1SpecPeriod', 'DeltaSRandShiftR3L1SpecDyn',...
            'DeltaSRandShiftR3D03', 'DeltaSRandShiftR3Gr', ...
            'DeltaSRandShiftGlOpportResh', ...
            'DeltaSConstGlOpportResh', 'DeltaSLarge_R3', 'DeltaSWideR3D03', 'DeltaSLarge_R3_Tune', ...
            'DeltaSRandShiftR1L1SpecGlOp'}
        for i = 1:size(ResourceRangesLowHigh, 1)
            p.ConfigurationID = ConfigurationID;
            ConfigurationID = ConfigurationID + 1;
            for k=1:NumReplica
                p.Monoculture = 0;
                p.Replica     = k;
                %
                RandRange =@(minVal, maxVal) sort(minVal + (maxVal -minVal) * rand(1, 2));
                %RandRange =@(minVal, maxVal) sort(minVal + (maxVal -minVal) * [0, 1]);
                switch SimulatID
                    case {'DeltaSRandShift', 'DeltaSRandShiftGlOpportResh'} %shifted resource supply points
                        p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies
                        p.ResLow(1) = rr(1);
                        p.ResHigh(1) = rr(2);  %
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies
                        p.ResLow(2) = rr(1);
                        p.ResHigh(2) = rr(2);  %
                    case {'DeltaSRandShiftR3', 'DeltaSRandShiftR3Reg', 'DeltaSRandShiftR3L', ...
                            'DeltaSRandShiftR3L1Spec', 'DeltaSRandShiftR3L1SpecPeriod', ...
                            'DeltaSRandShiftR3D03'} %3 resources
                        p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies R1
                        p.ResLow(1) = rr(1);
                        p.ResHigh(1) = rr(2);  %
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies for R2
                        p.ResLow(2) = rr(1);
                        p.ResHigh(2) = rr(2);  %
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies  for R3
                        p.ResLow(3) = rr(1);
                        p.ResHigh(3) = rr(2);  %
                        if strcmp(SimulatID, 'DeltaSRandShiftR3L1SpecPeriod')
                            p.SSupplActivity = @(t)((cos(2*pi*t/50)+1)/2).^2;
                        end
%                     case 
%                         p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
%                         p.ResLow(1:3) = ResourceRangesLowHigh(1);
%                         p.ResHigh(1:3) = ResourceRangesLowHigh(2);  %
                    case 'DeltaSRandShiftR3L1SpecDyn'
                    p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
                        p.ResLow(1) = ResourceRangesLowHigh(1);
                        p.ResHigh(1) = ResourceRangesLowHigh(2);  %
                        p.ResLow(2) = ResourceRangesLowHigh(1);
                        p.ResHigh(2) = ResourceRangesLowHigh(2);  %
                        p.ResLow(3) = ResourceRangesLowHigh(1);
                        p.ResHigh(3) = ResourceRangesLowHigh(2);  %
                                        
                    case 'DeltaSRandShiftR1L1SpecGlOp'
                        p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies R1
                        p.ResLow(1) = rr(1);
                        p.ResHigh(1) = rr(2);  %
                        rr = RandRange(ResourceRangesLowHigh(1), ResourceRangesLowHigh(2));   %the range of resource supplies for R2
                        p.ResLow(2) = rr(1);
                        p.ResHigh(2) = rr(2);  %
                    case {'DeltaSLarge_R3', 'DeltaSWideR3D03', 'DeltaSLarge_R3_Tune'}
                        p.ReSup = 'rSh'; %ResLow and ResMin defined for both resources 
                        p.ResLow(1:3) = ResourceRangesLowHigh(1);
                        p.ResHigh(1:3) = ResourceRangesLowHigh(2);
                    otherwise
                        p.ResLow = ResourceRangesLowHigh(i, 1);   %the range of resource supplies
                        p.ResHigh = ResourceRangesLowHigh(i, 2);  %
                end
                p.TraitLow = TraitRangesLowHigh(1, 1);   %the minimal R* value
                p.TraitHigh = TraitRangesLowHigh(1, 2);  %the maximal R* value
                
                
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
    case 'NIncr' %species number increases, but DeltaR* does not change
        for i = 1:length(SpPerSimulation)
            p.ConfigurationID = ConfigurationID;
            ConfigurationID = ConfigurationID + 1;
            for k=1:NumReplica
                p.Monoculture = 0;
                p.Replica     = k;
                %
                ind_RStar_sel = TraitRangesLowHigh(1, 1)  < RStar(:, 1)  & RStar(:, 1) < TraitRangesLowHigh(1, 2);
                p.TraitLow = TraitRangesLowHigh(1, 1);   %the minimal R* value
                p.TraitHigh = TraitRangesLowHigh(1, 2);  %the maximal R* value
                
                RStar_sel = RStar(ind_RStar_sel, :);
                specID_sel = p.SpecID_All(ind_RStar_sel);
                specID_sel = sort(specID_sel(randperm(length(specID_sel), SpPerSimulation(i))));
                
                p.specID        = specID_sel;        %here we take all species, but you can take you can define e.g. p.specID = [1, 5, 6]
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





