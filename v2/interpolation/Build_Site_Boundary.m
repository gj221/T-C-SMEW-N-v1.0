function site = Build_Site_Boundary(shape_file,feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One projected boundary covering every selected shapefile record.
%%
%% The notebook took the union of the selected fields as a single boundary,
%% fitted one kriging model to all the samples inside it and drew one
%% surface. This returns that union in the local UTM zone, shaped like a
%% grid setup so the same interpolation code serves both, but with no
%% rotation and no cell grid: coordinates stay true Easting and Northing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; feature_indices = []; end
script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(fullfile(v2_directory,'T&C_Code'));

[feature_polygons,metadata] = ...
    Read_Shapefile_Polygons(shape_file,feature_indices);
union_polygon = feature_polygons{1};
for k = 2:numel(feature_polygons)
    union_polygon = union(union_polygon,feature_polygons{k});
end
if isempty(union_polygon.Vertices) || area(union_polygon) <= 0
    error('Build_Site_Boundary:EmptyUnion', ...
        'The selected records must enclose positive area.');
end

% NaN separators keep the disjoint fields apart, so the vertices are passed
% through untouched rather than closed into one ring.
boundary_vertices = union_polygon.Vertices;
site.field_boundary_x_m = boundary_vertices(:,1);
site.field_boundary_y_m = boundary_vertices(:,2);
site.grid_rotation_degrees = 0;
site.projected_grid_origin_m = [0 0];
site.projected_crs_epsg = metadata.projected_crs_epsg;
site.projected_crs_name = metadata.projected_crs_name;
site.source_shapefile = shape_file;
site.source_feature_indices = metadata.feature_indices;
site.feature_polygons = feature_polygons;
site.union_polygon = union_polygon;
site.field_area_m2 = area(union_polygon);
site.region_count = numel(regions(union_polygon));

end
