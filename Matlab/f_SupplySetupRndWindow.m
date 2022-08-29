function S = f_SupplySetupRndWindow(t, ModParams)
%define resource supply distribution within a randomly choosen window 


res = length(ModParams.ResHigh);
S = NaN(ModParams.Ly, ModParams.Lx, res);
iRanges = 3;
for ir = 1:res
    Edgges = linspace(ModParams.ResLow(ir), ModParams.ResHigh(ir), iRanges + 1);
    iRndRange = randi(3);
    Low = Edgges( iRndRange);
    High = Edgges( iRndRange + 1);
    S(:, :, ir)  = Low + (High -Low) * rand(ModParams.Ly, ModParams.Lx);
end