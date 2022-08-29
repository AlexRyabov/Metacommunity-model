function [iFgNrNew, fg] = ResHeterogen_AddPlots(PlotID, iFgNr, Data)

%create a figure if it does not persist otherwiese clean it

fg = findobj( 'Type', 'Figure', 'Name', PlotID );
if isempty(fg)
    fg = figure(iFgNr);
    fg.Name = PlotID;
else
    figure(fg);
end
iFgNrNew = iFgNr + 1;
clf




switch PlotID
    case 'Fig2_ResPlane'
        %% Resource plane with consumtion vectors
        RStar_All = Data.RStar_All;  %R* values of species
        consumprate = Data.consumprate;  %R* values of species
        %set figure size
        fg.Position = [150 50  1048 433];
        subplot(1, 2, 1);
        f_Resource_Plane_Plot(RStar_All, consumprate, [100, 100], struct('Colors', true, 'LogScale', false, 'Number', 5));
        f_Lbls('Resource 1', 'Resource 2')
        subplot(1, 2, 2);
        f_Resource_Plane_Plot(RStar_All, consumprate, [1000, 1000], struct('Colors', true, 'LogScale', true, 'Number', 5));
        f_Lbls('Resource 1', 'Resource 2')
    case 'Fign_DivBiom'
        fg.Position = [150 50   1148  690];
end