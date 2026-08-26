function topo = Build_Projected_Field_Topology( ...
        field_polygon,cell_size_m,metadata,source_feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibrate and annotate one connected projected field-grid setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(cell_size_m)
    [cell_size_m,grid_calibration,grid_polygon] = ...
        Calibrate_Field_Grid_Cell_Size(field_polygon);
else
    [grid_polygon,grid_calibration] = ...
        Orient_Field_Polygon_Grid(field_polygon,cell_size_m);
    grid_calibration.mode = 'prescribed';
end
topo = Build_Polygon_Grid_Topology(grid_polygon,cell_size_m);
grid_calibration.grid_nx = topo.main_grid_nx;
grid_calibration.grid_ny = topo.main_grid_ny;
grid_calibration.grid_cell_count = topo.N_cells;
grid_calibration.simulated_area_m2 = topo.simulated_area_m2;
grid_calibration.field_area_m2 = topo.field_area_m2;
grid_calibration.field_fraction = topo.field_area_m2/topo.simulated_area_m2;
topo.grid_calibration = grid_calibration;
topo.grid_rotation_degrees = grid_calibration.rotation_degrees;
topo.projected_grid_origin_m = grid_calibration.projected_grid_origin_m;
[topo.cell_xc_projected_m,topo.cell_yc_projected_m] = ...
    Grid_Local_To_Projected(topo.cell_xc,topo.cell_yc,topo);
[topo.field_boundary_projected_x_m,topo.field_boundary_projected_y_m] = ...
    Grid_Local_To_Projected( ...
    topo.field_boundary_x_m,topo.field_boundary_y_m,topo);
topo.source_shapefile = metadata.source_shapefile;
topo.source_crs_name = metadata.source_crs_name;
topo.projected_crs_name = metadata.projected_crs_name;
topo.projected_crs_epsg = metadata.projected_crs_epsg;
topo.shapefile_feature_count = metadata.shapefile_feature_count;
topo.shapefile_included_feature_count = numel(metadata.feature_indices);
topo.source_feature_indices = source_feature_indices;

target_crs = projcrs(metadata.projected_crs_epsg);
[latitude,longitude] = projinv(target_crs, ...
    topo.field_boundary_projected_x_m,topo.field_boundary_projected_y_m);
valid = isfinite(latitude) & isfinite(longitude);
topo.field_latitude_bounds = [min(latitude(valid)) max(latitude(valid))];
topo.field_longitude_bounds = [min(longitude(valid)) max(longitude(valid))];

end
