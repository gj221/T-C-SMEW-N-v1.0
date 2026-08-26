function [setups,summary,figure_handle] = Plot_Field_Grid( ...
        shape_file,cell_size_m,feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot one calibrated grid and summary for each disconnected polygon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
code_directory = fullfile(v2_directory,'T&C_Code');
addpath(code_directory);

use_bruce_farm_default = nargin < 1 || isempty(shape_file);
if use_bruce_farm_default
    shape_file = fullfile(v2_directory,'Spatial inputs', ...
        'Bruce Farm Virginia','Whole farm', ...
        'AdamWilbourne_Control_Fields.shp');
end
if nargin < 2; cell_size_m = []; end
if nargin < 3
    feature_indices = [];
    if use_bruce_farm_default; feature_indices = [1 2]; end
end

setups = Build_Shapefile_Grid_Setups( ...
    shape_file,cell_size_m,feature_indices);
features = shaperead(shape_file);
field_name = strings(numel(setups),1);
source_indices = strings(numel(setups),1);
for k = 1:numel(setups)
    indices = setups(k).source_feature_indices;
    source_indices(k) = strjoin(string(indices),',');
    field_name(k) = feature_label(features,indices,k);
end

cell_size = reshape([setups.cell_size_m],[],1);
grid_columns = reshape([setups.main_grid_nx],[],1);
grid_rows = reshape([setups.main_grid_ny],[],1);
cell_count = reshape([setups.N_cells],[],1);
field_area = reshape([setups.field_area_m2],[],1);
simulated_area = reshape([setups.simulated_area_m2],[],1);
outside_area = simulated_area-field_area;
summary = table((1:numel(setups))',source_indices,field_name, ...
    cell_size,grid_columns,grid_rows,cell_count,field_area, ...
    field_area/4046.8564224,simulated_area,outside_area, ...
    outside_area./simulated_area, ...
    'VariableNames',{'Polygon','SourceFeatureIndices','Field', ...
    'CellSizeM','GridColumns','GridRows','CellCount','FieldAreaM2', ...
    'FieldAreaAcres','SimulatedAreaM2','OutsideAreaM2','OutsideFraction'});

column_count = ceil(sqrt(numel(setups)));
row_count = ceil(numel(setups)/column_count);
figure_handle = figure('Color','w','Position',[100 100 ...
    max(700,600*column_count) max(600,560*row_count)]);
layout = tiledlayout(figure_handle,row_count,column_count, ...
    'TileSpacing','compact','Padding','compact');
for k = 1:numel(setups)
    topo = setups(k);
    axes_handle = nexttile(layout);
    axes_handle.Tag = 'FieldSetupAxes';
    axes_handle.Color = 'w';
    axes_handle.XColor = 'k';
    axes_handle.YColor = 'k';
    hold(axes_handle,'on');
    for c = 1:topo.N_cells
        x = [topo.cell_x_lo(c) topo.cell_x_hi(c) ...
            topo.cell_x_hi(c) topo.cell_x_lo(c)];
        y = [topo.cell_y_lo(c) topo.cell_y_lo(c) ...
            topo.cell_y_hi(c) topo.cell_y_hi(c)];
        patch(axes_handle,x,y,topo.cell_overlap_fraction(c), ...
            'EdgeColor',[0.35 0.35 0.35],'LineWidth',0.7, ...
            'Tag','GridCell');
    end
    plot(axes_handle,topo.field_boundary_x_m,topo.field_boundary_y_m, ...
        '-','Color',[0.85 0.25 0.05],'LineWidth',2.5, ...
        'Tag',sprintf('FieldOutline%02d',k));
    axis(axes_handle,'equal','tight');
    box(axes_handle,'on');
    clim(axes_handle,[0 1]);
    xlabel(axes_handle,'Rotated grid x [m]','Color','k');
    ylabel(axes_handle,'Rotated grid y [m]','Color','k');
    title(axes_handle,sprintf('%s: %d by %d cells at %.2f m', ...
        field_name(k),topo.main_grid_nx,topo.main_grid_ny, ...
        topo.cell_size_m),'FontWeight','normal','Color','k');
    subtitle(axes_handle,sprintf( ...
        'Field %.0f m^2, grid %.0f m^2, outside %.0f m^2', ...
        topo.field_area_m2,topo.simulated_area_m2, ...
        topo.simulated_area_m2-topo.field_area_m2), ...
        'Color','k','FontSize',9);
end
colormap(figure_handle,parula);
colorbar_handle = colorbar(axes_handle,'Ticks',0:0.2:1);
colorbar_handle.Layout.Tile = 'east';
colorbar_handle.Label.String = 'Cell area inside polygon [-]';
colorbar_handle.Color = 'k';
colorbar_handle.Label.Color = 'k';
title(layout,sprintf('Independent full-grid setups for %d polygons', ...
    numel(setups)),'FontWeight','normal','Color','k');
xlabel(layout,'Each polygon is simulated and aggregated separately', ...
    'Color','k');

[~,shape_name] = fileparts(shape_file);
shape_name = regexprep(shape_name,'[^A-Za-z0-9]+','_');
feature_tag = '';
if ~isempty(feature_indices)
    parts = arrayfun(@(x) sprintf('%d',x),feature_indices, ...
        'UniformOutput',false);
    feature_tag = sprintf('_features%s',strjoin(parts,'_'));
end
if isempty(cell_size_m)
    resolution = 'geometry';
else
    resolution = sprintf('cell%gm',cell_size_m);
    resolution = strrep(resolution,'.','p');
end
output_name = sprintf('Field_Grid_%s%s_%s', ...
    shape_name,feature_tag,resolution);
figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
exportgraphics(figure_handle,fullfile(figure_directory,[output_name '.png']), ...
    'Resolution',300);
savefig(figure_handle,fullfile(figure_directory,[output_name '.fig']));
writetable(summary,fullfile(figure_directory,[output_name '_summary.csv']));

end

function label = feature_label(features,indices,setup_index)
if isfield(features,'Field')
    values = string({features(indices).Field});
    values = unique(values(strlength(values) > 0),'stable');
    if ~isempty(values)
        label = strjoin(values,', ');
        return
    end
end
label = sprintf('Polygon %d',setup_index);
end
