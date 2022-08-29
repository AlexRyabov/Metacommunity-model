function f_Resource_Plane_Plot(RstAll, ConsRates, XYMax, Params)
    function PlotIsoc(Rst, XYMax, Clr)
        x1 = linspace(Rst(1), XYMax(1), 10);
        y1 = x1*0 + Rst(2);
        y2 = linspace(XYMax(2), Rst(2), 10);
        x2 = y2*0 + Rst(1);
        plot([x2, x1], [y2, y1], 'Color', Clr, 'LineWidth', 1, 'LineStyle', '-');
    end
    function   PlotCons(XY0, XYMax, Slope, Clr)
        DXY = XYMax-XY0;
        x1 = logspace(log10(XY0(1)), log10(XYMax(1)), 100);
        y1 = XY0(2) + Slope * (x1 - XY0(1));
        if y1(end) >  XYMax(2)
            y1 = logspace(log10(XY0(2)), log10(XYMax(2)), 100);
            x1 = (y1 - XY0(2))/Slope  + XY0(1);
        end
        plot(x1, y1, 'Color', Clr, 'LineWidth', 2)
    end
%sort Rstr values from small Rst 1 to large Rst
[RstAll,index] = sortrows(RstAll);
ConsRates = ConsRates(index, :);

if ~isnan(Params.Number)
 ind = round(linspace(1, size(RstAll, 1), Params.Number));
 RstAll = RstAll(ind, :);
 ConsRates = ConsRates(ind, :);
end

if Params.Colors
    Colors = linspecer(size(RstAll, 1), 'sequential');
    %Colors = f_Clrs_fresh(size(RstAll, 1)); % linspecer(size(RstAll, 1), 'sequential');
else
    Colors = zeros((size(RstAll, 1)), 3);
end
Slopes = ConsRates(:, 2)./ConsRates(:, 1);
PlotIsoc(RstAll(1, :), XYMax, Colors(1, :))
hold on
for i = 2:size(RstAll, 1)
    PlotIsoc(RstAll(i, :), XYMax, Colors(i, :))
    %find intersection points
    XY0 =  [RstAll(i, 1), RstAll(i-1, 2)];
    PlotCons(XY0, XYMax, Slopes(i-1), Colors(i-1, :));
    PlotCons(XY0, XYMax, Slopes(i), Colors(i, :));
end
if Params.LogScale
    f_Set_scales({'logx', 'logy'});
end
end
