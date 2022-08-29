% Calculate Community Characteristics
% date: 23.10.2015
% author: Dahu
% dependencies: NSpecCompAdv_Cluster_Run_Model

%% calculate community values

% create params vector from main script
Params_vector = zeros(size(Params));
for  iParam = 1:length(Params)
    Params_vector(iParam) = Params(iParam).(ParamName_field);
end

%% start of calculation

%% calculate relative biomasses
% save in strucure object - Results

% determine max biomass of all species at a certain timestep
for d2t = 1:length(Results(iParam).TSpan)
    for i2sp = 1:sp
        SpMaxBioCol = max(Results(iParam).N_Dyn(:, :, i2sp, d2t)); % max biomass as row vector over the columns of each species
        SpMaxBio(i2sp) = max(SpMaxBioCol); % max biomass of the previous vector saved for each species
    end
    ComMaxBio(d2t) = max(SpMaxBio); % maximal biomass of a undefined species at a certain timestep  
end

% calculation of the relative biomasses
for iParam = 1:length(Params)
    for d3t = 1:length(Results(iParam).TSpan)
        for i3sp = 1:sp
            Results(iParam).RelBiomass(:,:,i3sp,d3t) = Results(iParam).N_Dyn(:,:,i3sp,d3t)/ComMaxBio(d3t);            
        end
    end 
end

%% different characteristics

SpeciesBiom_Param = zeros(length(Params), sp);
for  iParam = 1:length(Params)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate the global biomass at the beginning the simulation of each species
    StartingBiomass = zeros(Results(iParam).Params.sp, 1);
    for s = 1:Results(iParam).Params.sp
        Sp_BegDistrib = Results(iParam).N_Dyn(:, :, s, 1);
        StartingBiomass(s)  = sum(Sp_BegDistrib(:));
    end
    Results(iParam).Sp_Total_Biomass_Start = sum(Results(iParam).N_Dyn(:, :, :, end), 3);
    Results(iParam).StartingBiomass =  StartingBiomass;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate the final global biomass of each species at the end of the
    % simulation
    FinalBiomass = zeros(Results(iParam).Params.sp, 1);
    for s = 1:Results(iParam).Params.sp
        Sp_FinalDistrib = Results(iParam).N_Dyn(:, :, s, end);
        FinalBiomass(s)  = sum(Sp_FinalDistrib(:));
    end
    Results(iParam).Sp_Total_Biomass = sum(Results(iParam).N_Dyn(:, :, :, end), 3);
    Results(iParam).FinalBiomass =  FinalBiomass;
       
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate the dynamics of the global biomass
    % used in plot: 1, 3 (Sp_Distrib)
    GlobalBiomass_dyn = zeros(Results(iParam).Params.sp,  length(Results(iParam).TSpan));
    for t = 1:length(Results(iParam).TSpan)
        for s = 1:Results(iParam).Params.sp
            Sp_Distrib = Results(iParam).N_Dyn(:, :, s, t);
            GlobalBiomass_dyn(s, t)  = sum(Sp_Distrib(:));
        end
    end
    Results(iParam).GlobalBiomass_dyn =  GlobalBiomass_dyn;
    Results(iParam).Sp_Distrib = Sp_Distrib;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate biomass of a species for a certain dispersal rate at the
    % end of the simulation
    % used in plot: 2
    for s = 1:Results(iParam).Params.sp
        SpBiomassPerCell = Results(iParam).N_Dyn(:, :, s, end);
        SpBiomassPerGrid(s) = sum(SpBiomassPerCell(:));        
    end
    Results(iParam).SpBiomassPerGrid = SpBiomassPerGrid;
    SpeciesBiom_Param(iParam, :) = SpBiomassPerGrid'; % sp biomass per grid at the end of the simulation
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate relative biomass of a species for a certain dispersal rate at the
    % end of the simulation
    for s = 1:Results(iParam).Params.sp
        SpRelBiomassPerCell = Results(iParam).RelBiomass(:, :, s, end);
        SpRelBiomassPerGrid(s) = sum(SpRelBiomassPerCell(:));        
    end
    Results(iParam).SpRelBiomassPerGrid = SpRelBiomassPerGrid;
    SpeciesRelBiom_Param(iParam, :) = SpRelBiomassPerGrid'; % sp biomass per grid at the end of the simulation
    
     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate biomass of a species for a certain dispersal rate at the
    % beginning of the simulation
    for s = 1:Results(iParam).Params.sp
        SpBiomPerCell = Results(iParam).N_Dyn(:, :, s, 1);
        SpBiomPerGridBeg(s) = sum(SpBiomPerCell(:));        
    end
    Results(iParam).SpBiomPerGridBeg = SpBiomPerGridBeg;
    SpeciesBiom_ParamBeg(iParam, :) = SpBiomPerGridBeg'; % sp biomass per grid at the beginning of the simulation
   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate of the final global biomass of all species at the end of
    % time range
    % used in plot: 5
    % calculation of the final biomass of each species 
    for s = 1:Results(iParam).Params.sp
        Sp_FinalDistrib = Results(iParam).N_Dyn(:, :, s, end);
        FinalBiomass(s)  = sum(Sp_FinalDistrib(:));
    end
    % calculation of the final global biomass of all species at the end
    % of the simulation time
    GlobalBiomassAllSp = sum(FinalBiomass);
    Results(iParam).GlobalBiomassAllSp = GlobalBiomassAllSp;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % calculate of the standardized biomass
    % of species at the end of the simulation
    MaxBioMassPerSp = max(FinalBiomass);
    FinBioPerSpStand = bsxfun(@rdivide, FinalBiomass, MaxBioMassPerSp); 
    Results(iParam).FinBioPerSpStand = FinBioPerSpStand; 
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
 
end % end of loop iParam

%% calculation of community characteristics out of the iParam loop
%  calculation of biodiversity indices
% now only results

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculate the standardized biomass of species per dispersal rate at the
% end of the simulation
    
% pre-calculation
MaxBioMassPerSpGrid = max(SpeciesBiom_Param);
MaxBioMassPerSpGridTrans = MaxBioMassPerSpGrid';
SpBiom_ParamTrans = SpeciesBiom_Param';
% calculation of the standardized biomass
FinBioPerSpPerDispStand = bsxfun(@rdivide, SpBiom_ParamTrans, MaxBioMassPerSpGridTrans);
Results(iParam).FinBioPerSpPerDispStand =  FinBioPerSpPerDispStand;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% calculate the shannon index - H' (H'=-sum(n(i)/n)*ln(n(i)/n))

% calculate the global biomass of each dispersal rate
for n = 1:length(Diff)
    GloBioDispRa(n) = sum(SpeciesBiom_Param(n,:));
    GloBioDispRaTrans = GloBioDispRa';
end
Results(iParam).GloBioDispRaTrans = GloBioDispRaTrans;
% start with the calculation
ShannonPi = bsxfun(@rdivide, SpeciesBiom_Param, GloBioDispRaTrans); % relative abundance of a species
ShannonLog = (ShannonPi.*(log(ShannonPi)))';
ShannonIndex = -sum(ShannonLog);
Results(iParam).ShannonIndex = ShannonIndex;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% calculate the shannon index - H' (H'=-sum(n(i)/n)*ln(n(i)/n))

% relative global biomass of each dispersal rate
for n = 1:length(Diff)
    GloRelBioDispRa(n) = sum(SpeciesRelBiom_Param(n,:));
    GloRelBioDispRaTrans = GloRelBioDispRa';
end
Results(iParam).GloRelBioDispRaTrans = GloRelBioDispRaTrans;
% start with the calculation
ShannonPiRel = bsxfun(@rdivide, SpeciesRelBiom_Param, GloRelBioDispRaTrans); % relative abundance of a species
ShannonLogRel = (ShannonPiRel.*(log(ShannonPiRel)))';
ShannonIndexRel = -sum(ShannonLogRel);
Results(iParam).ShannonIndexRel = ShannonIndexRel;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculate the simpson index - D (D=n(i)*(n(i)-1)/n*(n-1))
% 
% for n = 1:length(Diff)
%     Simpson_Numerator = SpeciesBiom_Param*(bsxfun(@minus,GloBioDispRaTrans, 1))
% end
% Simpson_Fraction = GloBioDispRaTrans*(bsxfun(@minus, GloBioDispRaTrans, 1))
% 
% 
% 
% Results(iParam).SimpsonIndex = SimpsonIndex;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculate the eveness (Pielou) - J (J=H'/H'max)

MaxShannonIndex = max(ShannonIndex);
Eveness = ShannonIndex./MaxShannonIndex;
Results(iParam).Eveness = Eveness;

%----------------------------------------------
    % calculation of species rank abundance
    for s = 1:Results(iParam).Params.sp;
        % Berechnung 
    end