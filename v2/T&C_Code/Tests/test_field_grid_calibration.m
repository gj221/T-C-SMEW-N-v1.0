function test_field_grid_calibration()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell width follows the polygon envelope geometry without a cell cap.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
addpath(code_directory);

rectangle = polyshape([0 300 300 0],[0 0 100 100]);
[cell_size_m,calibration,oriented] = ...
    Calibrate_Field_Grid_Cell_Size(rectangle);
assert(abs(cell_size_m-75) < 1e-12);
assert(calibration.grid_nx == 4 && calibration.grid_ny == 2);
assert(calibration.grid_cell_count == 8);
assert(abs(calibration.envelope_area_m2-30000) < 1e-12);
assert(abs(calibration.envelope_perimeter_m-800) < 1e-12);
assert(abs(calibration.outside_grid_area_m2-15000) < 1e-9);
assert(abs(area(oriented)-area(rectangle)) < 1e-9);

rotated_rectangle = rotate(rectangle,31,[0 0]);
[rotated_size,rotated_calibration,rotated_polygon] = ...
    Calibrate_Field_Grid_Cell_Size(rotated_rectangle);
rotated_topo = Build_Polygon_Grid_Topology(rotated_polygon,rotated_size);
assert(rotated_topo.N_cells == 8);
assert(abs(rotated_calibration.rotation_degrees-31) < 1e-8);
assert(abs(rotated_calibration.outside_grid_area_m2-15000) < 1e-7);

shape_file = fullfile(fileparts(code_directory), ...
    'Spatial inputs','Bruce Farm Virginia','applied_yield_shp.shp');
bruce = Build_Shapefile_Grid_Topology(shape_file,[]);
assert(bruce.main_grid_nx == 2 && bruce.main_grid_ny == 13);
assert(bruce.N_cells == 26);
assert(bruce.cell_size_m > 49 && bruce.cell_size_m < 51);
assert(strcmp(bruce.grid_calibration.mode,'polygon_geometry'));
assert(abs(bruce.grid_calibration.outside_grid_area_m2-( ...
    bruce.simulated_area_m2-bruce.field_area_m2)) < 1e-8);
assert(isfinite(bruce.grid_rotation_degrees));
assert(all(isfinite(bruce.projected_grid_origin_m)));
[projected_x,projected_y] = Grid_Local_To_Projected([0;1],[0;0],bruce);
assert(abs(projected_x(1)-bruce.projected_grid_origin_m(1)) < 1e-9);
assert(abs(projected_y(1)-bruce.projected_grid_origin_m(2)) < 1e-9);
assert(abs(hypot(diff(projected_x),diff(projected_y))-1) < 1e-9);

prescribed = Build_Shapefile_Grid_Topology(shape_file,100);
assert(prescribed.N_cells < 21);
assert(strcmp(prescribed.grid_calibration.mode,'prescribed'));

slender = polyshape([0 1000 1000 0],[0 0 10 10]);
[~,slender_calibration] = Calibrate_Field_Grid_Cell_Size(slender);
assert(slender_calibration.grid_cell_count > 24);
assert(slender_calibration.grid_nx >= 2 && slender_calibration.grid_ny >= 2);

outer = polyshape([0 100 100 0],[0 0 100 100]);
inner = polyshape([25 75 75 25],[25 25 75 75]);
donut = subtract(outer,inner);
island = polyshape([140 160 160 140],[40 40 60 60]);
multipart = union(donut,island);
[multipart_size,~,multipart_oriented] = ...
    Calibrate_Field_Grid_Cell_Size(multipart);
multipart_topo = Build_Polygon_Grid_Topology( ...
    multipart_oriented,multipart_size);
assert(multipart_topo.field_region_count == 2);
assert(abs(sum(multipart_topo.cell_overlap_area_m2)-area(multipart)) < 1e-8);
assert(all(multipart_topo.cell_overlap_fraction >= 0));
assert(all(multipart_topo.cell_overlap_fraction <= 1));

rounded_case = polyshape( ...
    [234.5308 233.7740 207.1565 152.7239 10.5309 14.1687], ...
    [72.5813 55.4079 18.1021 18.0929 85.4820 117.6316]);
[~,rounded_calibration] = Calibrate_Field_Grid_Cell_Size(rounded_case);
assert(rounded_calibration.grid_cell_count == 10);

fprintf('test_field_grid_calibration: all checks passed\n');

end
