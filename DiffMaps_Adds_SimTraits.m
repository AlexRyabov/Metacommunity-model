function DiffMaps_Adds_SimTraits(DiffDist, Traits, Fig)
%plot the correlation between similarities between species and
%TraitDistance
if ~exist('Fig')
   Fig = figure();
end

figure(Fig)
clf
%find distances 
D = squareform(pdist(Traits));


%scatter plot similarity vs distance
scatter(D(:), DiffDist(:), '.');
ind = DiffDist(:)>0;


%modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2 + b(4) * x.^3;
%lm = fitnlm(D(ind), DiffDist(ind), modelfun, [1, 1, -1, 1]);
%modelfun = @(b,x) b(1) + b(2) * x + b(3) * x.^2;
%lm = fitnlm(D(ind), DiffDist(ind), modelfun, [1, 1, -1]);
modelfun = @(b,x) b(1) + b(2) * x;
lm = fitnlm(D(ind), DiffDist(ind), modelfun, [1, 1]);
x = linspace(min(D(ind)), max(D(ind)), 100);
hold on
plot(x, lm.predict(x'), 'Color', [210, 245, 60]/256, 'LineWidth', 3)

title(['R^2_{adj}=' NS(lm.Rsquared.Adjusted)])
f_Lbls('Distance in R^* space', 'Diffusion distance')