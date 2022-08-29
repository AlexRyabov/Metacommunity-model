%%Make MDS of similatiry matrix and compare it with species traits
function  DiffMaps_Adds_SimMDS(SimMatrix, indPos, Dimensions, Traits, Fig)
if ~exist('Fig')
    Fig = figure();
end
%%
Dd=1-SimMatrix;
Y = MDS(Dd,Dimensions);

Mean = mean(Traits);
Dist = Traits - repmat(Mean, size(Traits, 1), 1);
Dist = sqrt(sum(Dist.^2, 2));

figure(Fig)
clf
subplot(2, 2, 1)
scatter3(Y(1,:), Y(2,:), Y(3,:), 20, Dist, 'o', 'filled');
colormap(jet)

subplot(2, 2, 2)
scatter3(Y(1,:), Y(2,:), Y(3,:), 20, Traits(:, 1), 'o', 'filled');
colormap(jet)
subplot(2, 2, 3)
scatter3(Y(1,:), Y(2,:), Y(3,:), 20, Traits(:,2), 'o', 'filled');
colormap(jet)
subplot(2, 2, 4)
scatter3(Y(1,:), Y(2,:), Y(3,:), 20, Traits(:, 3), 'o', 'filled');
colormap(jet)



end