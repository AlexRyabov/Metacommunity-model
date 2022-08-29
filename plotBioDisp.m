% plot 2
% barplot of final biomass of each species
% author: Daniel Huber

%%
figure(2);

if length(Diff) < 24
    
    for iParam = 1:length(Params)
        subplot(4,6,iParam);
        hold on;
        bar(Results(iParam).SpBiomassPerGrid, 'EdgeColor', 'k');
        title([ParamName_plot ' = ' num2str(Params(iParam).(ParamName_field))]);
        xlabel('species id');
        ylabel('biomass');
        xlim([0, sp+1]);
        hold off;
    end    
else
    for iParam = 1:length(Params)
        subplot(5,5,iParam);
        hold on;
        bar(Results(iParam).SpBiomassPerGrid, 'EdgeColor', 'k');
        title([ParamName_plot ' = ' num2str(Params(iParam).(ParamName_field))]);
        xlabel('species id');
        ylabel('biomass');
        xlim([0, sp+1]);
        hold off;
    end    
end

% %special case - plot for only one barplot
%     bar(Results(23).SpBiomassPerGrid, ); % set special dispersal rate 
%     title([ParamName_plot ' = ' num2str(Params(23).(ParamName_field))], 'FontSize', 18);
%     xlabel('species id', 'FontSize', 18, 'FontWeight', 'bold');
%     ylabel('biomass', 'FontSize', 18, 'FontWeight', 'bold');
%     xlim([0, sp+1]);
%     hold off;