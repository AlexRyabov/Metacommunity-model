%% Define R* as in Metacom
close all
rng(1)
sp = 2000;
Dims = 3;
RStar = rand(sp, Dims);        %Random values
RStar = RStar./repmat(sum(RStar, 2), 1, Dims);  %normilize to have (R1* + R2* + R3*).

%Note that we get a higher proportion of intermediate species
%compared with extreme species
TraitHigh = 1;
TraitLow = 10;
RStar = RStar * (TraitHigh - TraitLow) + TraitLow;

ind = RStar(:, 1)<6 | RStar(:, 1)>7;
RStar = RStar(ind, :);

RStar = sortrows(RStar, 2);

figure(1)
clf
%subplot(2, 2, 1)
scatter3(RStar(:, 1), RStar(:,2), RStar(:, 3),  10, RStar(:, 3), '.')
axis([0, 15, 0, 15, 0, 15])
f_Lbls3D('R1', 'R2', 'R3')
colormap(jet)
colorbar
view(gca,[145.54 48.24]);

%% Diff map it
k_max = 10;
%[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'VectProd', true, k_max);
[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'NormzdEuc', true, k_max);
%[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'StndzdEuc', true, k_max);
%[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'NormzdGaus', true, k_max);
%[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'NormzdPearson', true, k_max);
%[ev, aEV, aSClnd, indPos] = f_DM_DMit(RStar, 'NormzdSpearman', true, k_max);

%%
%plot eigval
figure(6)
scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3),20, RStar(:, 1), '.')
colormap(jet)
hold on
ind = RStar(:, 1) < 5;
scatter3(aEV(ind, 1), aEV(ind, 2), aEV(ind, 3),20, RStar(ind, 1), 'o', 'filled')
hold off
c = colorbar
c.Label.String = 'R*_1'
f_Lbls3D('ev_1', 'ev_2', 'ev_3')

figure(7)
scatter3(aEV(:, 3), aEV(:, 4), aEV(:, 5),20, RStar(:, 1), '.')
colormap(jet)
hold on
ind = RStar(:, 1) < 5;
scatter3(aEV(ind, 3), aEV(ind, 4), aEV(ind, 5),20, RStar(ind, 1), 'o', 'filled')
hold off
c = colorbar
c.Label.String = 'R*_1'
f_Lbls3D('ev_3', 'ev_4', 'ev_5')
%%
figure(8)
D1 = 3;
D2 = 6;
D3 = 8;
scatter3(aEV(:, D1), aEV(:, D2), aEV(:, D3),20, RStar(:, 1), '.')
colormap(jet)
hold on
ind = RStar(:, 1) > 6;
scatter3(aEV(ind, D1), aEV(ind, D2), aEV(ind, D3),20, RStar(ind, 1), 'o', 'filled')
hold off
c = colorbar
c.Label.String = 'R*_1'
f_Lbls3D(['ev_{' NS(D1) '}'], ['ev_{' NS(D2) '}'], ['ev_{' NS(D3) '}'])

%%


figure()
clf
%subplot(2, 2, 1)
[~, ev1_order] = sort(aEV(:, 1));
scatter3(RStar(:, 1), RStar(:,2), RStar(:, 3), 10, ev1_order, '.')
axis([0, 15, 0, 15, 0, 15])
f_Lbls3D('R1', 'R2', 'R3')
colorbar
