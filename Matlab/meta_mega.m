%%meta mega analysis

Calculate = 0;

%%initialize all result variables as cell arrays
vec_tResults = {};
vec_tParams  = {}; 
vec_tFuncInd = {}; 
vec_tRes_BEF = {};
vec_tRst = {};

%%one variant
[vec_tResults{end + 1}, vec_tParams{end + 1}, vec_tFuncInd{end + 1}, vec_tRes_BEF{end + 1}, vec_tRst{end + 1}] =...
    NSpecCompAdv_Cluster_Run_Model('NutNet', 'Rst_EffSize_R1', Calculate);

[vec_tResults{end + 1}, vec_tParams{end + 1}, vec_tFuncInd{end + 1}, vec_tRes_BEF{end + 1}, vec_tRst{end + 1}] =...
    NSpecCompAdv_Cluster_Run_Model('NutNet', 'Rst_EffSize_R1_glopp', Calculate);


[vec_tResults{end + 1}, vec_tParams{end + 1}, vec_tFuncInd{end + 1}, vec_tRes_BEF{end + 1}, vec_tRst{end + 1}] =...
    NSpecCompAdv_Cluster_Run_Model('NutNet', 'Rst_EffSize_R2', Calculate);

[vec_tResults{end + 1}, vec_tParams{end + 1}, vec_tFuncInd{end + 1}, vec_tRes_BEF{end + 1}, vec_tRst{end + 1}] =...
    NSpecCompAdv_Cluster_Run_Model('NutNet', 'Rst_EffSize_R2_glopp', Calculate);

tRst_R1 = vec_tRst{1};
tRst_R2 = vec_tRst{2};

tParams_R1 = vec_tParams{1};
tParams_R2 = vec_tParams{2};

tResults_R1 = vec_tResults{1}; 
tResults_R2 = vec_tResults{2};


%% plot figures
Fg1 = figure(1);
set(Fg1, 'Position', [100, 100, 900, 600]);
% plot species abundance dynamics for fertilized plots 
plot(1:size(tResults_R1.N_SpaceTimeDyn{63},4),squeeze(sum(sum(tResults_R1.N_SpaceTimeDyn{63}(:,:,1,:)))))
hold on
for i = 2:tParams_R1.sp(63)
    plot(1:size(tResults_R1.N_SpaceTimeDyn{63},4),squeeze(sum(sum(tResults_R1.N_SpaceTimeDyn{63}(:,:,i,:)))))
end

Fg2 = figure(2);
set(Fg2, 'Position', [600, 100, 900, 600]);
% plot species abundance dynamics for control plots
plot(1:size(tResults_R1.N_SpaceTimeDyn{11},4),squeeze(sum(sum(tResults_R1.N_SpaceTimeDyn{11}(:,:,1,:)))))
hold on
for i = 2:tParams_R1.sp(11)
    plot(1:size(tResults_R1.N_SpaceTimeDyn{11},4),squeeze(sum(sum(tResults_R1.N_SpaceTimeDyn{11}(:,:,i,:)))))
end


Fg3 = figure(3);
set(Fg3, 'Position', [800, 200, 900, 600])
% find indices for certain table entries (here: only species mixtures, with nutrient addition, for smaller trait range)
ind = (tParams_R1.Monoculture == 0 & tParams_R1.TreatmentPlot == 1 & tParams_R1.TraitHigh == 8);
trt = ind(size(tParams_R1,1)/2+1:size(tParams_R1,1),:); % fertilized plots

% plot calculated and specified R* values (both resources) for all species
%for i = 1:length(trt)
RStar_R1 = tParams_R1.RStar{trt(1)}(:,1);
RStar_R2 = tParams_R2.RStar{trt(1)}(:,2);
ES_R1 = transpose(tRst_R1.Rst{trt(1)});
ES_R2 = transpose(tRst_R2.Rst{trt(1)});

col = [1:tParams_R1.sp(trt(1)),2,1]; % define colors by index value
c = linspace(1,10,length(RStar_R1));
scatter(RStar_R1,RStar_R2,[],c,'filled')
hold on
scatter(ES_R1,ES_R2,[],c,'filled')
hold off
%end










% 
% %%another variant
% vSimulatID = {'Rst_EffSize_R1', 'Rst_EffSize_R2'};
% for iSim = 1:length(vSimulatID)
%     [tResults(end + 1), tParams(end + 1), tFuncInd(end + 1), tRes_BEF(end + 1), tRst(end + 1)] =... 
%         NSpecCompAdv_Cluster_Run_Model('NutNet', vSimulatID{iSim}, Calculate);
% end


%%processing
% figure(1)
% for iSim = 1:numel(vec_tResults)
%     tResults = vec_tResults{iSim};
%     tParams = vec_tParams{iSim};
%     tRst = vec_tRst{iSim};
%     plot(Rstars...)
%         hold on;
%     
%     
% end





