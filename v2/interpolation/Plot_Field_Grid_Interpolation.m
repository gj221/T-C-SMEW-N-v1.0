function [setups,summary,figure_handle,interpolations] = ...
        Plot_Field_Grid_Interpolation(shape_file,cell_size_m, ...
        feature_indices,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the calibrated grid of each disconnected polygon, with the kriged
%% attribute surface underneath it.
%%
%% Called with a shapefile alone this is the plain grid figure: retained
%% full cells shaded by the fraction of their area inside the polygon, the
%% exact boundary over the top, and a per-polygon summary CSV. Given a
%% sample table, an attribute and a stage, each polygon that holds three or
%% more usable samples also carries its interpolated surface, with the
%% measurements drawn on the same colour scale and partial boundary cells
%% outlined. Polygons without data at that stage stay in grid form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(fullfile(v2_directory,'T&C_Code'));

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

parser = inputParser;
parser.addParameter('Points',table());
parser.addParameter('Attribute','');
parser.addParameter('SampleStage','');
parser.addParameter('StageColumn','Sample Type');
parser.addParameter('XColumn','lng (X)');
parser.addParameter('YColumn','lat (Y)');
parser.addParameter('VariogramModel','exponential');
parser.addParameter('VariogramParameters',struct());
parser.addParameter('GridResolution',200);
parser.addParameter('OutputName','');
parser.parse(varargin{:});
settings = parser.Results;
interpolating = ~isempty(settings.Points) && ...
    ~isempty(settings.Attribute) && ~isempty(settings.SampleStage);

setups = Build_Shapefile_Grid_Setups( ...
    shape_file,cell_size_m,feature_indices);
features = shaperead(shape_file);
field_name = strings(numel(setups),1);
source_indices = strings(numel(setups),1);
for k = 1:numel(setups)
    indices = setups(k).source_feature_indices;
    source_indices(k) = strjoin(string(indices),',');
    field_name(k) = Field_Setup_Label(features,indices,k);
end
summary = Field_Grid_Summary_Table(setups,source_indices,field_name);

interpolations = cell(numel(setups),1);
value_limits = [Inf -Inf];
if interpolating
    for k = 1:numel(setups)
        try
            interpolations{k} = Interpolate_Field_Attribute( ...
                settings.Points,setups(k),settings.Attribute, ...
                settings.SampleStage,'StageColumn',settings.StageColumn, ...
                'XColumn',settings.XColumn,'YColumn',settings.YColumn, ...
                'VariogramModel',settings.VariogramModel, ...
                'VariogramParameters',settings.VariogramParameters, ...
                'GridResolution',settings.GridResolution);
            value_limits = [min(value_limits(1), ...
                interpolations{k}.value_limits(1)) ...
                max(value_limits(2),interpolations{k}.value_limits(2))];
        catch interpolation_error
            if ~startsWith(interpolation_error.identifier, ...
                    'Interpolate_Field_Attribute:')
                rethrow(interpolation_error)
            end
            interpolations{k} = [];
            fprintf('%s at %s: %s\n',field_name(k), ...
                settings.SampleStage,interpolation_error.message);
        end
    end
end
has_surface = ~cellfun(@isempty,interpolations);
if ~any(has_surface); value_limits = [0 1]; end
if diff(value_limits) <= 0; value_limits = value_limits+[-0.5 0.5]; end

% Ribbon-shaped fields keep their equal aspect, so stacking them one per
% row is the only layout that does not leave most of the figure empty.
swap_axes = false(numel(setups),1);
for k = 1:numel(setups)
    swap_axes(k) = boundary_extent(setups(k).field_boundary_y_m) > ...
        boundary_extent(setups(k).field_boundary_x_m);
end
elongated = all(arrayfun(@(k) display_aspect(setups(k),swap_axes(k)) > 4, ...
    1:numel(setups)));
if elongated
    column_count = 1;
    row_count = numel(setups);
    figure_size = [1500 max(620,380*row_count+220)];
else
    column_count = ceil(sqrt(numel(setups)));
    row_count = ceil(numel(setups)/column_count);
    figure_size = [max(700,600*column_count) max(600,560*row_count)];
end
figure_handle = figure('Color','w','Position',[100 100 figure_size]);
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
    % Laying the long axis across the tile keeps the map readable instead
    % of a sliver a few pixels wide.
    if ~isempty(interpolations{k})
        draw_surface(axes_handle,interpolations{k},swap_axes(k));
    end
    if ~isempty(interpolations{k})
        cell_mode = 'mesh';
    elseif any(has_surface)
        cell_mode = 'blank';
    else
        cell_mode = 'overlap';
    end
    draw_grid_cells(axes_handle,topo,cell_mode,swap_axes(k));
    [boundary_x,boundary_y] = orient_pair(topo.field_boundary_x_m, ...
        topo.field_boundary_y_m,swap_axes(k));
    plot(axes_handle,boundary_x,boundary_y, ...
        '-','Color',[0.85 0.25 0.05],'LineWidth',2.5, ...
        'Tag',sprintf('FieldOutline%02d',k));
    if ~isempty(interpolations{k})
        draw_samples(axes_handle,interpolations{k},swap_axes(k));
    end
    axis(axes_handle,'equal','tight');
    box(axes_handle,'on');
    if any(has_surface)
        clim(axes_handle,value_limits);
    else
        clim(axes_handle,[0 1]);
    end
    if swap_axes(k)
        xlabel(axes_handle,'Rotated grid y [m]','Color','k');
        ylabel(axes_handle,'Rotated grid x [m]','Color','k');
    else
        xlabel(axes_handle,'Rotated grid x [m]','Color','k');
        ylabel(axes_handle,'Rotated grid y [m]','Color','k');
    end
    title(axes_handle,sprintf('%s: %d retained cells at %.2f m', ...
        field_name(k),topo.N_cells,topo.cell_size_m), ...
        'FontWeight','normal','Color','k');
    subtitle(axes_handle,tile_subtitle(topo,interpolations{k}), ...
        'Color','k','FontSize',9);
end
colormap(figure_handle,parula);
colorbar_handle = colorbar(axes_handle,'Location','eastoutside');
colorbar_handle.Layout.Tile = 'east';
colorbar_handle.Color = 'k';
colorbar_handle.Label.Color = 'k';
if any(has_surface)
    colorbar_handle.Label.String = settings.Attribute;
    header_lines = { ...
        sprintf('%s at %s over %d calibrated polygon grids', ...
        settings.Attribute,settings.SampleStage,numel(setups)), ...
        sprintf(['Ordinary kriging (%s); each polygon is simulated and ' ...
        'aggregated separately'],char(settings.VariogramModel))};
else
    colorbar_handle.Ticks = 0:0.2:1;
    colorbar_handle.Label.String = 'Cell area inside polygon [-]';
    header_lines = { ...
        sprintf('Independent retained full-cell setups for %d polygons', ...
        numel(setups)), ...
        'Each polygon is simulated and aggregated separately'};
end
% A tiledlayout title shares the top of the first tile with that tile's own
% title, and an equal-aspect ribbon pushes its title right to the edge. A
% reserved header strip keeps the two apart whatever the aspect turns out
% to be.
header_height = min(0.16,120/figure_size(2));
layout.OuterPosition = [0 0 1 1-header_height];
annotation(figure_handle,'textbox',[0 1-header_height 1 header_height], ...
    'String',header_lines,'EdgeColor','none','Color','k', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',12,'FitBoxToText','off','Tag','FigureHeader');

output_name = settings.OutputName;
if isempty(output_name)
    output_name = Field_Grid_Output_Name(shape_file,cell_size_m, ...
        feature_indices,settings.Attribute,settings.SampleStage);
end
figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
exportgraphics(figure_handle,fullfile(figure_directory, ...
    [output_name '.png']),'Resolution',300);
savefig(figure_handle,fullfile(figure_directory,[output_name '.fig']));
writetable(summary,fullfile(figure_directory,[output_name '_summary.csv']));

end

function draw_surface(axes_handle,interpolation,swap_axes)
[surface_x,surface_y] = meshgrid(interpolation.grid_x,interpolation.grid_y);
[surface_x,surface_y] = orient_pair(surface_x,surface_y,swap_axes);
surface(axes_handle,surface_x,surface_y,zeros(size(surface_x)), ...
    interpolation.prediction,'EdgeColor','none','FaceColor','flat', ...
    'Tag','KrigedSurface');
end

function draw_grid_cells(axes_handle,topo,cell_mode,swap_axes)
x = [topo.cell_x_lo topo.cell_x_hi topo.cell_x_hi topo.cell_x_lo];
y = [topo.cell_y_lo topo.cell_y_lo topo.cell_y_hi topo.cell_y_hi];
[x,y] = orient_pair(x,y,swap_axes);
vertices = [x(:) y(:)];
faces = [(1:topo.N_cells)' (topo.N_cells+1:2*topo.N_cells)' ...
    (2*topo.N_cells+1:3*topo.N_cells)' ...
    (3*topo.N_cells+1:4*topo.N_cells)'];
% The mesh is drawn at every cell count so the figure shows the simulation's
% true resolution; a fine grid is meant to read as fine, not be smoothed
% away or dropped. At tens of thousands of cells the line spacing is below
% one export pixel, so opaque lines would merge into a flat grey wash and
% bury the surface under it instead of showing it through a fine mesh. A
% partial edge alpha keeps every line real (it is still there in the
% vector .fig at full resolution) while letting density read as texture.
% The true field boundary and any cell straddling it stay fully opaque, so
% the one line that is never ambiguous is the one drawn boldest.
cell_edge_color = [0.35 0.35 0.35];
mesh_edge_alpha = min(1,max(0.06,4000/topo.N_cells));
switch cell_mode
    case 'mesh'
        patch(axes_handle,'Faces',faces,'Vertices',vertices, ...
            'FaceColor','none','EdgeColor',cell_edge_color, ...
            'EdgeAlpha',mesh_edge_alpha,'LineWidth',0.35,'Tag','GridCell');
        partial = find(topo.cell_overlap_fraction < 1-1e-12);
        if ~isempty(partial)
            patch(axes_handle,'Faces',faces(partial,:), ...
                'Vertices',vertices,'FaceColor','none', ...
                'EdgeColor',[0.15 0.15 0.15],'LineWidth',0.9, ...
                'Tag','BoundaryCell');
        end
    case 'blank'
        % The colour axis belongs to the attribute on the other tiles, so a
        % polygon with no data at this stage is drawn plain rather than
        % mapped onto a scale that does not describe it.
        patch(axes_handle,'Faces',faces,'Vertices',vertices, ...
            'FaceColor',[0.88 0.88 0.88],'EdgeColor',cell_edge_color, ...
            'EdgeAlpha',mesh_edge_alpha,'LineWidth',0.35,'Tag','GridCell');
    otherwise
        patch(axes_handle,'Faces',faces,'Vertices',vertices, ...
            'FaceVertexCData',topo.cell_overlap_fraction, ...
            'FaceColor','flat','EdgeColor',cell_edge_color, ...
            'EdgeAlpha',mesh_edge_alpha,'LineWidth',0.35,'Tag','GridCell');
end
end

function draw_samples(axes_handle,interpolation,swap_axes)
[sample_x,sample_y] = orient_pair(interpolation.samples.x_grid_m, ...
    interpolation.samples.y_grid_m,swap_axes);
scatter(axes_handle,sample_x,sample_y,45,interpolation.sample_values, ...
    'filled','MarkerEdgeColor','k','LineWidth',0.8,'Tag','SamplePoint');
end

function [a,b] = orient_pair(x,y,swap_axes)
if swap_axes
    a = y;
    b = x;
else
    a = x;
    b = y;
end
end

function aspect = display_aspect(topo,swap_axes)
x_extent = boundary_extent(topo.field_boundary_x_m);
y_extent = boundary_extent(topo.field_boundary_y_m);
if swap_axes
    aspect = y_extent/x_extent;
else
    aspect = x_extent/y_extent;
end
end

function extent = boundary_extent(coordinate)
coordinate = coordinate(isfinite(coordinate));
extent = max(coordinate)-min(coordinate);
end

function text = tile_subtitle(topo,interpolation)
text = sprintf('Field %.0f m^2, grid %.0f m^2, %.2f%% inside', ...
    topo.field_area_m2,topo.simulated_area_m2, ...
    100*topo.field_area_m2/topo.simulated_area_m2);
if isempty(interpolation)
    text = sprintf('%s, no usable samples at this stage',text);
else
    text = sprintf('%s, %d samples, %.3g to %.3g',text, ...
        interpolation.inside_sample_count, ...
        interpolation.value_limits(1),interpolation.value_limits(2));
end
end

