function S = f_SupplySetup(t, ModParams)
res = length(ModParams.ResHigh);
S = NaN(ModParams.Ly, ModParams.Lx, res);
for ir = 1:res
    S(:, :, ir)  = ModParams.ResLow(ir) + (ModParams.ResHigh(ir) -ModParams.ResLow(ir)) * rand(ModParams.Ly, ModParams.Lx);
end