function MetaCom_plotResults(Params, Tsol, Xsol, fg_Handlers)
fg_res = fg_Handlers{1};clf;
fg_sp = fg_Handlers{2};clf;
fg_tr_biom = fg_Handlers{3};clf;

Xsol = reshape(Xsol, length(Tsol), Params.Ly, Params.Lx, Params.nSp + Params.nRes);
Xsol = permute(Xsol, [2, 3, 4, 1]); %set time as the last dimension in the array;
N_Dyn = Xsol(:, :, 1:Params.nSp, :);
R_Dyn = Xsol(:, :, Params.nSp + (1:Params.nRes), :);
figure(fg_res)
%% plot resources
Res = [];
R = R_Dyn(:, :, :, end);
MaxR = max(R(:));
for ir = 1:Params.nRes
    subplot(2, 2, ir)
    pcolor(R_Dyn(:, :, ir, end));
    Resi = R_Dyn(:, :, ir, end);
    Res = [Res, Resi(:)];
    colorbar
    shading flat
    title(['Resource' NS(ir)])
    caxis([0, MaxR])
end
 subplot(2, 2, 4)
 scatter3(Res(:, 1), Res(:, 2), Res(:, 3));
 
 
%% plot species 
figure(fg_sp)
N = N_Dyn(:, :, :, end);
MaxN = max(N(:));
for iSp = 1:min(25, Params.nSp)
    subplot(5, 5, iSp)
    pcolor(N_Dyn(:, :, iSp, end));
    caxis([0, 0.1*MaxN])
    shading flat
    title(['Sp' NS(iSp)])
end

%% spatial similarity
ix = 5;
jy = 5;
figure(fg_tr_biom)
iSP = 1;
colormap(jet(200))
for j = jy+2:-2:jy-2
    for i = ix-2:2:ix+2
        subplot(3, 3, iSP)
        iSP = iSP + 1;
        ind = squeeze( N_Dyn(j, i, :, end)) > 1e-5;
        scatter3(Params.RStar(ind , 1), Params.RStar(ind , 2), log(squeeze(N_Dyn(j, i, ind , end))), 60, log(squeeze(N_Dyn(j, i, ind , end))), '.');
        view(2)
    end
end



