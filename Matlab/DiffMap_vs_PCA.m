
rng(5)
%rng(3)
m = 20000;
x = sort(rand(1, m));

sigma = 0.01;
F = @(x, z) exp(-(x-z).^2/(2 * sigma^2));


n = 50; %number of species
Traits = sort(rand(n, 1));
%Traits = linspace(0, 1, 10);

%fill a matrix of observations
mObs = NaN(n, m);
for iz = 1:n
    mObs(iz, :) = F(x, Traits(iz));
end
%mObs(mObs(:) <0.00001) = 0;
fg1 = figure(1);
fg1.Position(3:4) = [ 672 916];
clf
tiledlayout('flow')
nexttile
clrs = parula(n);
for iz = 1:n
    plot(x, mObs(iz, :), 'Color', clrs(iz, :));
    hold on
end

title('Distribution')
f_Lbls('x', 'Biomass')
DM_Distance = 'Spearman';   %use euqlidean distances 
%DM_Distance = 'Pearson';   %use euqlidean distances 

%find similarity
%aS = (corr(mObs','Type', 'Spearman', 'Rows', 'pairwise')+1)/2;
aS = (corr(mObs','Type', DM_Distance, 'Rows', 'pairwise')+1)/2;
for i = 1:n
    aS(i, i) = 0;
end
ax = nexttile
sp_ax = 1:n;
pcolor_central(sp_ax, sp_ax, aS)
colormap(1-gray)
shading flat
title('Similarity, S')
f_Lbls('Species ID', 'Species ID')
pos = ax.Position;
colorbar
drawnow
ax.Position = pos;
%make Diff map 
k_max = 10;
[ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(mObs, DM_Distance, false, k_max);
  [U_pca, Usc2_pca, score_pca, lambdas_pca, S_pca, D2_pca, explained_pca, mu_pca] = f_PCA(mObs');
%  ev_pca = lambdas_pca;
  aEV_pca = Usc2_pca;
%  DiffDist_pca = D2_pca; 
%  aS_pca = S_pca;
figure(1)
nexttile
scatter(aEV(:, 1), aEV(:, 2), 20, clrs, 'filled')
f_Lbls('Trait_1', 'Trait_2')
title('Diffusion map')
nexttile
scatter(aEV_pca(:, 1), aEV_pca(:, 2), 20, clrs, 'filled')
f_Lbls('Trait_1', 'Trait_2')
title('PCA')

%% find distances
transp = 0.1;
D0 = squareform(pdist(Traits, 'euclidean'));
nexttile
s1 = scatter(D0(:), DiffDist(:), 20, 'filled');
s1.MarkerFaceAlpha = transp;
f_Lbls('Orig distance', 'Reconst. distance');
nexttile
s1 = scatter(D0(:), D2_pca(:), 20, 'filled')
s1.MarkerFaceAlpha = transp;
f_Lbls('Orig distance', 'Reconst. distance');
