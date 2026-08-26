function setups = Build_Shapefile_Grid_Setups( ...
        shape_file,cell_size_m,feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build one independent full-grid setup per connected polygon region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; cell_size_m = []; end
if nargin < 3; feature_indices = []; end
[feature_polygons,metadata] = ...
    Read_Shapefile_Polygons(shape_file,feature_indices);

field_polygon = polyshape();
for k = 1:numel(feature_polygons)
    field_polygon = union(field_polygon,feature_polygons{k});
end
field_regions = regions(field_polygon);
if isempty(field_regions)
    error('Build_Shapefile_Grid_Setups:EmptyField', ...
        'The selected polygon union has no positive-area regions.');
end

setup_cell = cell(numel(field_regions),1);
sort_key = zeros(numel(field_regions),3);
for k = 1:numel(field_regions)
    region = field_regions(k);
    region_tol = max(1e-6,1e-9*area(region));
    member = false(size(feature_polygons));
    for j = 1:numel(feature_polygons)
        member(j) = area(intersect(region,feature_polygons{j})) > region_tol;
    end
    source_feature_indices = metadata.feature_indices(member);
    setup_cell{k} = Build_Projected_Field_Topology( ...
        region,cell_size_m,metadata,source_feature_indices);
    [center_x,center_y] = centroid(region);
    sort_key(k,:) = [min(source_feature_indices) center_x center_y];
end

[~,order] = sortrows(sort_key,[1 2 3]);
setups = vertcat(setup_cell{order});
for k = 1:numel(setups)
    setups(k).setup_index = k;
    setups(k).setup_count = numel(setups);
end

end
