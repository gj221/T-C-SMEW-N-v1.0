function test_field_grid_plot()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Without samples the interpolation figure is the plain grid figure: one
%% subplot and one summary row per retained whole-site polygon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
addpath(code_directory,fullfile(v2_directory,'interpolation'));

shape_file = fullfile(v2_directory,'Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
[setups,summary,figure_handle] = ...
    Plot_Field_Grid_Interpolation(shape_file,[],[1 2]);
figure_cleanup = onCleanup(@() close(figure_handle));

assert(numel(setups) == 2);
assert(height(summary) == 2);
assert(isequal(summary.SourceFeatureIndices,["1";"2"]));
assert(isequal(string(summary.Field),["Control 2";"Control 1"]));
assert(isequal(summary.SourceVertexCount,[19;30]));
assert(isequal(summary.ProjectedPolygonVertexCount,[18;29]));
assert(all(summary.ProjectionRoundTripMaxM < 1e-6));
assert(all(summary.SourceToPolygonMaxM < 1e-9));
assert(all(summary.CellCount == reshape([setups.N_cells],[],1)));
assert(all(summary.GridSlotCount == reshape([setups.grid_slot_count],[],1)));
assert(all(summary.FieldFraction > 0.95));
assert(max(abs(summary.FieldAreaM2- ...
    reshape([setups.field_area_m2],[],1))) < 1e-9);
expected_outside = reshape([setups.simulated_area_m2],[],1)- ...
    reshape([setups.field_area_m2],[],1);
assert(max(abs(summary.OutsideAreaM2-expected_outside)) < 1e-9);

setup_axes = findobj(figure_handle,'Type','axes','Tag','FieldSetupAxes');
field_outlines = findobj(figure_handle,'Type','line', ...
    '-regexp','Tag','FieldOutline');
grid_cells = findobj(figure_handle,'Type','patch','Tag','GridCell');
assert(numel(setup_axes) == 2);
assert(numel(field_outlines) == 2);
assert(numel(grid_cells) == numel(setups));
assert(sum(arrayfun(@(x) size(x.Faces,1),grid_cells)) == ...
    sum([setups.N_cells]));
for k = 1:numel(setups)
    outline = findobj(figure_handle,'Type','line', ...
        'Tag',sprintf('FieldOutline%02d',k));
    polygon = polyshape(outline.XData,outline.YData,'Simplify',false);
    assert(abs(area(polygon)-setups(k).field_area_m2) < 1e-6);
end

figure_directory = fullfile(v2_directory,'figures');
base_name = ['Field_Grid_Interpolation_' ...
    'AdamWilbourne_Control_Fields_features1_2_geometry'];
assert(isfile(fullfile(figure_directory,[base_name '.png'])));
assert(isfile(fullfile(figure_directory,[base_name '.fig'])));
assert(isfile(fullfile(figure_directory,[base_name '_summary.csv'])));

fprintf('test_field_grid_plot: all checks passed\n');

clear figure_cleanup
end
