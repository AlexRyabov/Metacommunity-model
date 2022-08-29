% DaHu
% plot 1
% plot of the biomass dynamic of each species over simulation time
figure(1);
% set(Fig1, 'Position', [100, 100, 400, 800]); % set figure to a certain
% position

if length(Diff) > 24
    
    for iParam = 1:length(Params)
        subplot(5,5,iParam);
        clrs = jet(Results(iParam).Params.sp);
        
        hold on
        for s = 1:Results(iParam).Params.sp
            plot(Results(iParam).TSpan, Results(iParam).GlobalBiomass_dyn(s, :), 'Color', clrs(s, :))
        end
        title([ParamName_plot '=' num2str(Params(iParam).(ParamName_field))]);
        %    legend(sp,'Location','southeast');
        hold off
        axis tight
        
    end
    
else
    
    for iParam = 1:length(Params)
        subplot(4,6,iParam);
        clrs = jet(Results(iParam).Params.sp);
        
        hold on
        for s = 1:Results(iParam).Params.sp
            plot(Results(iParam).TSpan, Results(iParam).GlobalBiomass_dyn(s, :), 'Color', clrs(s, :))
        end
        title([ParamName_plot '=' num2str(Params(iParam).(ParamName_field))]);
        %    legend(sp,'Location','southeast');
        hold off
        axis tight
        
    end
    
end