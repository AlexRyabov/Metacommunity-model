Calc = 0;

SimulatIDs = {'Sigma=Mu^1.5', 'Sigma=Mu^2',    'Sigma=Mu^2.5',      'Sigma=Mu^3', ...
    'Sigma=Mu^1.5_dyn50', 'Sigma=Mu^2_dyn50', 'Sigma=Mu^2.5_dyn50', 'Sigma=Mu^3_dyn50', ...
    'Sigma=Mu^1.5_dyn05', 'Sigma=Mu^2_dyn05', 'Sigma=Mu^2.5_dyn05', 'Sigma=Mu^3_dyn05'};
%SimulatIDs = {  'Sigma=Mu^3_dyn50', ...
%    'Sigma=Mu^2.5_dyn05', 'Sigma=Mu^3_dyn05'};

for i = 1:length(SimulatIDs)
    try
        i
        NSpecCompAdv_Cluster_Run_Model('ResHeterogen', SimulatIDs{i}, Calc)
    catch ME
       display(ME) 
    end
end


%NSpecCompAdv_Cluster_Run_Model('ResHeterogen', 'Sigma=Mu^2', Calc)
%NSpecCompAdv_Cluster_Run_Model('ResHeterogen', 'Sigma=Mu^2.5', Calc)
