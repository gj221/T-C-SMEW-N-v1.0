function topo = Build_Shapefile_Grid_Topology(shape_file,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Shapefile_Grid_Topology
%%
%% Read polygon records, merge overlaps and duplicates, project the union
%% to its local WGS84 UTM zone, then build the full equal-square grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; cell_size_m = []; end
if isstring(shape_file) && isscalar(shape_file); shape_file = char(shape_file); end
if ~(ischar(shape_file) && isrow(shape_file) && ~isempty(shape_file)) || ...
        ~isfile(shape_file)
    error('Build_Shapefile_Grid_Topology:MissingShapefile', ...
        'shape_file must name an existing .shp file.');
end
if ~strcmpi(file_extension(shape_file),'.shp')
    error('Build_Shapefile_Grid_Topology:BadFileType', ...
        'shape_file must have a .shp extension.');
end

info = shapeinfo(shape_file);
if ~startsWith(info.ShapeType,'Polygon','IgnoreCase',true)
    error('Build_Shapefile_Grid_Topology:BadGeometry', ...
        'The shapefile must contain polygon geometry.');
end
source_crs = info.CoordinateReferenceSystem;
if isempty(source_crs)
    error('Build_Shapefile_Grid_Topology:MissingCRS', ...
        'The shapefile must define its coordinate reference system.');
end

features = shaperead(shape_file);
if isempty(features)
    error('Build_Shapefile_Grid_Topology:EmptyShapefile', ...
        'The shapefile contains no polygon features.');
end

latitude = cell(numel(features),1);
longitude = cell(numel(features),1);
for k = 1:numel(features)
    if isa(source_crs,'geocrs')
        latitude{k} = features(k).Y;
        longitude{k} = features(k).X;
    elseif isa(source_crs,'projcrs')
        [latitude{k},longitude{k}] = projinv( ...
            source_crs,features(k).X,features(k).Y);
    else
        error('Build_Shapefile_Grid_Topology:UnsupportedCRS', ...
            'The shapefile CRS must be geographic or projected.');
    end
end

all_latitude = [latitude{:}];
all_longitude = [longitude{:}];
valid = isfinite(all_latitude) & isfinite(all_longitude);
all_latitude = all_latitude(valid);
all_longitude = all_longitude(valid);
if isempty(all_latitude) || any(all_latitude < -80 | all_latitude > 84) || ...
        range(all_longitude) >= 6
    error('Build_Shapefile_Grid_Topology:UnsupportedExtent', ...
        'The field must fit inside the latitude and width limits of one local UTM projection.');
end

latitude_center = mean(all_latitude);
longitude_center = mean(all_longitude);
utm_zone = min(60,max(1,floor((longitude_center+180)/6)+1));
if latitude_center >= 0
    projected_crs_epsg = 32600+utm_zone;
else
    projected_crs_epsg = 32700+utm_zone;
end
target_crs = projcrs(projected_crs_epsg);

field_polygon = polyshape();
for k = 1:numel(features)
    [x_m,y_m] = projfwd(target_crs,latitude{k},longitude{k});
    feature_polygon = polyshape(x_m,y_m,'Simplify',true);
    if area(feature_polygon) > 0
        field_polygon = union(field_polygon,feature_polygon);
    end
end
if isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Build_Shapefile_Grid_Topology:EmptyField', ...
        'The projected polygon union has no positive area.');
end

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
topo.source_shapefile = shape_file;
topo.source_crs_name = char(source_crs.Name);
topo.projected_crs_name = char(target_crs.Name);
topo.projected_crs_epsg = projected_crs_epsg;
topo.shapefile_feature_count = numel(features);
topo.field_latitude_bounds = [min(all_latitude) max(all_latitude)];
topo.field_longitude_bounds = [min(all_longitude) max(all_longitude)];

end

function extension = file_extension(file_path)
[~,~,extension] = fileparts(file_path);
end
