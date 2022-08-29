function f_DistribPlot(AbundanceXYi, XScale, YScale, Trait, PlotDetail, Power, Name)
arguments
    AbundanceXYi double
    XScale (1,:) double
    YScale (1,:) double
    Trait  double = NaN
    PlotDetail  = false
    Power double = 2
    Name =  'Species'
end

if isnan(Trait)
    Trait = 1:size(AbundanceXYi, 3); %no sorting of species
end

%find total abundance in each cell 
switch ndims(AbundanceXYi)
    case {1,2}
        AbundanceXY = AbundanceXYi;
        FindNorm = false;
    case 3
        FindNorm = true;
        AbundanceXY = sum(AbundanceXYi, 3);
        %Sotr by trait
        [Trait, SortIndex] = sort(Trait);
        AbundanceXYi = AbundanceXYi(:, :, SortIndex);
    otherwise
        error('AbundanceXYi should be 2D or 3D matrix')
end
[XSize, YSize, SpNum] = size(AbundanceXYi);

%set NaNs where Abuandance < 1e-10
FilledCells_idx = AbundanceXY < 1e-10;
NormAbundXYi = NaN(size(AbundanceXYi));
for iSp = 1:SpNum
    tmpAb = AbundanceXYi(:, :, iSp);
    tmpAb(FilledCells_idx) = NaN;
    %Normilize
    if FindNorm
        tmpAb = tmpAb./AbundanceXY;
    end
    NormAbundXYi(:, :, iSp) = tmpAb;
end

%Take to power
NormAbundXYi = NormAbundXYi.^Power;

%Generare species colors (distinct or continuous)
Colors_i = linspecer(SpNum);
%Get color distribution for each species
SpColorDistr = NaN(XSize, YSize, SpNum, 3);
for iSp = 1:SpNum
    for iCl = 1:3
        SpColorDistr(:, :, iSp, iCl) = NormAbundXYi(:, :, iSp) * Colors_i(iSp, iCl);
    end
end

if ~PlotDetail
    %Get Average color
    SpAvrColorDistr = squeeze(sum(SpColorDistr(:, :, :, :), 3));
    imagesc(XScale, YScale, SpAvrColorDistr);
    %change axis direction
    ax = gca;
    ax.YAxis.Direction = 'normal';
else
    %make 4x4 subplots for each species with largest biomass
    %Get species biomasses
    Abundance_i = squeeze(sum(sum(AbundanceXYi, 1), 2));
    %find 16 species with maximal biomass
    [~, FilledCells_idx] = maxk(Abundance_i, 16);
    SnNum2Plot = length(FilledCells_idx);
    CCol = ceil(sqrt(SnNum2Plot));
    CRows = ceil(SnNum2Plot/CCol);
    cMax = 0;
    for iSp = 1:SnNum2Plot
        subplot(CCol, CRows, iSp)
        pcolor(XScale, YScale, NormAbundXYi(:, :, FilledCells_idx(iSp)));
        shading flat
        colorbar
        ca = caxis;
        caxis([0, ca(2)])
        title([Name '=', NS(Trait(FilledCells_idx(iSp)))]);
        
    end
end
end