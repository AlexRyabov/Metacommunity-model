%% Modeling spatial resource competiton
%% Alexey Ryabov 2013-2017
%% This script plots final results

%% plot eveness versus biomass
fg1 = figure(1); %f_MakeFigure(1, [10, 10, 1200, 400]);
clf;
clrs = linspecer(size(ResourceRangesLowHigh, 1));   %use different colors for different resource ranges
subplots = size(TraitRangesLowHigh, 1);
leg = {};
for iTr = 1  %% :size(TraitRangesLowHigh, 1); use only the largerst trait range
    for iRes = 1:size(ResourceRangesLowHigh)
        ind = tParams.ResLow == ResourceRangesLowHigh(iRes, 1) & tParams.TraitLow == TraitRangesLowHigh(iTr, 1);
        x = []; x2 = []; x3 = []; y= [];
        for i = 1:length(ind)
            if ind(i)
                x = [x; tFuncInd.Pielou_even_loc_xy{i}];
                x2 = [x2; tFuncInd.Richness_xy{i}];
                x3 = [x3; tFuncInd.Simpson_ESN_loc_xy{i}];
                y = [y; tFuncInd.Biomass_xy{i}];
            end
        end
        ax1 = subplot(1, 2, 1);
        semilogy(x(:), y(:), '.', 'Color', clrs(iRes, :), 'MarkerSize', 10);
        f_Lbls('Eveness, local', 'Biomass, local');
        hold on;
        leg{iRes, iTr} = ['\Delta R* = ' num2str( TraitRangesLowHigh(iTr, 2) - TraitRangesLowHigh(iTr, 1)) ...
            ', \Delta R = ' num2str( ResourceRangesLowHigh(iRes, 2) - ResourceRangesLowHigh(iRes, 1))];

%         ax2 = subplot(1, 3, 2);
%         semilogy(x2(:) + (iRes/6-0.5)*0.5, y(:), '.', 'Color', clrs(iRes, :), 'MarkerSize', 10);
%         hold on;
%         f_Lbls('Richness, local', 'Biomass, local');
% 
        ax3 = subplot(1, 2, 2);
        semilogy(x3(:), y(:), '.', 'Color', clrs(iRes, :), 'MarkerSize', 10);
        hold on;
        f_Lbls('Diversity, local', 'Biomass, local');
    end
end
%  subplot(1, 3, 3);
%% 
  legend(leg{:, 1}, 'Location', 'northeastoutside');
  drawnow;
   w = mean([ax1.Position(3), ax3.Position(3)]) * 1.1;
   ax1.Position(3) = w;
%   ax2.Position(3) = w;
   ax3.Position(3) = w;
   offset = 0.05;
   ax1.Position(1) = ax1.Position(1) - offset;
%   ax2.Position(1) = ax2.Position(1) - offset;
   ax3.Position(1) = ax3.Position(1) - offset;

