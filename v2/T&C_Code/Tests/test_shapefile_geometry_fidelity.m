function test_shapefile_geometry_fidelity()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Projected polygons preserve the supplied WGS84 boundary vertices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
addpath(code_directory);

shape_file = fullfile(v2_directory,'Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
features = shaperead(shape_file);
[feature_polygons,metadata] = ...
    Read_Shapefile_Polygons(shape_file,[1 2]);

assert(metadata.projected_crs_epsg == 32617);
assert(isequal(metadata.source_finite_vertex_count,[19 30]));
assert(isequal(metadata.projected_polygon_vertex_count,[18 29]));
assert(all(metadata.projection_roundtrip_max_error_m < 1e-6));
assert(all(metadata.source_to_polygon_max_error_m < 1e-9));

target_crs = projcrs(metadata.projected_crs_epsg);
for k = 1:2
    [x_m,y_m] = projfwd(target_crs,features(k).Y,features(k).X);
    valid = isfinite(x_m) & isfinite(y_m);
    vertices = feature_polygons{k}.Vertices;
    source_distance_m = zeros(sum(valid),1);
    source_index = find(valid);
    for j = 1:numel(source_index)
        i = source_index(j);
        source_distance_m(j) = min(hypot( ...
            vertices(:,1)-x_m(i),vertices(:,2)-y_m(i)));
    end
    assert(max(source_distance_m) < 1e-9);
end

fprintf('test_shapefile_geometry_fidelity: all checks passed\n');

end
