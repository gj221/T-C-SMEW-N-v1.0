function test_field_grid_topology()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exact polygon-grid intersections and Bruce Farm shapefile projection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
addpath(code_directory);

field_polygon = polyshape([0 150 150 50 50 0], ...
    [0 0 50 50 150 150]);
topo = Build_Polygon_Grid_Topology(field_polygon,100);

assert(topo.N_cells == 3);
assert(topo.main_grid_nx == 2 && topo.main_grid_ny == 2);
assert(all(abs(topo.cell_area_m2-10000) < 1e-12));
assert(max(abs(topo.cell_overlap_area_m2-[7500;2500;2500])) < 1e-9);
assert(max(abs(topo.cell_overlap_fraction-[0.75;0.25;0.25])) < 1e-12);
assert(abs(sum(topo.cell_overlap_area_m2)-area(field_polygon)) < 1e-9);
assert(topo.N_edges == 2);
assert(all(abs(topo.edge_face_width_m-100) < 1e-12));
assert(all(abs(topo.edge_face_distance_m-100) < 1e-12));
assert(abs(sum(topo.outer_face_width_m)-800) < 1e-12);
assert(all(strcmp(topo.cell_region,'main')));
assert(all(topo.cell_overlap_area_m2 > 0));

shape_file = fullfile(fileparts(code_directory), ...
    'Spatial inputs','Bruce Farm Virginia','applied_yield_shp.shp');
bruce = Build_Shapefile_Grid_Topology(shape_file,[],1);

assert(bruce.projected_crs_epsg == 32617);
assert(bruce.shapefile_feature_count == 12);
assert(bruce.field_region_count == 1);
assert(bruce.field_area_m2 > 0);
assert(abs(sum(bruce.cell_overlap_area_m2)-bruce.field_area_m2) < 1e-6);
assert(bruce.grid_calibration.field_fraction > 0.95);
assert(abs(bruce.simulated_area_m2- ...
    bruce.N_cells*bruce.cell_size_m^2) < 1e-9);
assert(all(bruce.cell_overlap_fraction >= 0));
assert(all(bruce.cell_overlap_fraction <= 1));
assert(all(bruce.cell_overlap_area_m2 > 0));

expect_error(@() Build_Polygon_Grid_Topology(field_polygon,0), ...
    'Build_Polygon_Grid_Topology:BadCellSize');
expect_error(@() Build_Polygon_Grid_Topology(polyshape(),100), ...
    'Build_Polygon_Grid_Topology:EmptyField');
expect_error(@() Grid_Cell_Count(100,0), ...
    'Grid_Cell_Count:BadCellSize');
expect_error(@() Grid_Cell_Count(-1,100), ...
    'Grid_Cell_Count:BadExtent');

fprintf('test_field_grid_topology: all checks passed\n');

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
