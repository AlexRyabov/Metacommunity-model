
MainModelParams.ReshuffleStep = 50;%
MainModelParams.TraitRangesLowHigh = [1, 10];
MainModelParams.Train     = 1;  %what fraction of data will used for diffusion map
MainModelParams.Validate  = 1;  %what fraction of data will used for diffusion map

Calculate =0;


%% Main paper simulation
%800 replica with random variations in the cloud of supply points. Stable, take the last
%values 3 resources. Specialist are the best competitors for one resource
%best variant
%   MainModelParams.Train     = 688/688;  %what fraction of data will used for diffusion map
%   MainModelParams.Validate  = 100/688;  %what fraction of data will used for diffusion map
%   NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3L1Spec', Calculate, MainModelParams);%

 %1000 replica with the same params as DeltaSRandShiftR3L1Spec, but with
 %reshafling of resources
  MainModelParams.Train     = 80/100;  %what fraction of data will used for diffusion map
  MainModelParams.Validate  = 20/100;  %what fraction of data will used for diffusion map
  NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3L1SpecDyn', Calculate, MainModelParams);%
% 

 %100 replica with the same params as DeltaSRandShiftR3L1Spec, but with
 %periodic oscillations 
%  MainModelParams.Train     = 80/100;  %what fraction of data will used for diffusion map
%  MainModelParams.Validate  = 20/100;  %what fraction of data will used for diffusion map
%  NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3L1SpecPeriod', Calculate, MainModelParams);%
% 

 
 %1000 replica with for gleaner opportunist in periodic env
%  rng(1);
%  MainModelParams.Train     = 80/100;  %what fraction of data will used for diffusion map
%  MainModelParams.Validate  = 20/100;  %what fraction of data will used for diffusion map
%  NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR1L1SpecGlOp', Calculate, MainModelParams);%



%60 replica with random variations in the cloud of supply points but regular R*. Stable, take the last
%values 3 resources. 
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3Reg', Calculate, MainModelParams);%


%100 replica with wide resource ranges. Stable, take the last
%values 3 resources. Diffusion of resources = 0.3
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSWideR3D03', Calculate, MainModelParams);%

%60 replica with random variations in the cloud of supply points. Stable, take the last
%values 3 resources. Diffusion of resources = 0.3
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3D03', Calculate, MainModelParams);%





%Random variations in the cloud of supply points. Stable, take the last
%values 3 resources, Gradient.
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftR3Gr', Calculate, MainModelParams);%

%A larger variability in resources, grid 30x31, small dispersal. Stable, take the last
%values 3 resources D=0.01
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSLargerSmallD_R3', Calculate, MainModelParams);%

%A larger variability in resources, grid 100x101, small dispersal. Stable, take the last
%values 3 resources D=0.1, Grid 100x101
%very unlear pattersn, too much noise, probably need to add diffusion of
%resources
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSLarge_R3', Calculate, MainModelParams);%

%%tune resource diffusion rate to get a smooth distribution of resources
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSLarge_R3_Tune', Calculate, MainModelParams);%

%Random variations in the cloud of supply points.  Permanent random addition of resources.
%Gleaner Opportunist trade off
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSRandShiftGlOpportResh', Calculate, MainModelParams);%

%large resource and trait range. Permanent random addition of resources.
%Gleaner Opportunist trade off
%NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaSConstGlOpportResh', Calculate, MainModelParams);%



% %for sp = [5,7, 10, 15, 20, 25, 30, 40 50, 75, 100, 200, 250, 300, 350, 400
% for sp = [50, 75, 100, 200, 250, 300, 350, 400]
%     MainModelParams.NumberOfSpecies = sp;
%     %[tResults, tParams, tFuncInd, tRes_BEF] = ...
%         NSpecCompAdv_Cluster_Run_Model('DiffMaps', 'DeltaRsIncr', Calculate, MainModelParams);
% end



