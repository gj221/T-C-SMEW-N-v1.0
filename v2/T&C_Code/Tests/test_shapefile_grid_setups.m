function test_shapefile_grid_setups()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Disconnected shapefile polygons use independent simulation grids.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
addpath(code_directory);

shape_file = fullfile(v2_directory,'Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
features = shaperead(shape_file);
assert(strcmp(features(3).Field,'Control 3'));

setups = Build_Shapefile_Grid_Setups(shape_file,[],[1 2]);
assert(numel(setups) == 2);
assert(all([setups.setup_index] == 1:2));
assert(all([setups.setup_count] == 2));
assert(isequal(setups(1).source_feature_indices,1));
assert(isequal(setups(2).source_feature_indices,2));
assert(abs(setups(1).field_area_m2-8410.290644) < 0.1);
assert(abs(setups(2).field_area_m2-15282.728509) < 0.1);
assert(abs(sum([setups.field_area_m2])-23693.019153) < 0.1);
assert(all([setups.field_region_count] == 1));
for k = 1:numel(setups)
    assert(abs(sum(setups(k).cell_overlap_area_m2)- ...
        setups(k).field_area_m2) < 1e-6);
    assert(setups(k).simulated_area_m2 < 1e6);
    assert(strcmp(setups(k).grid_calibration.mode,'polygon_geometry'));
end

expect_error(@() Build_Shapefile_Grid_Setups(shape_file,[],0), ...
    'Read_Shapefile_Polygons:BadFeatureIndices');
expect_error(@() Build_Shapefile_Grid_Setups(shape_file,[],[1 1]), ...
    'Read_Shapefile_Polygons:BadFeatureIndices');

fprintf('test_shapefile_grid_setups: all checks passed\n');

end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
