function topo = Build_Polygon_Grid_Topology(field_polygon,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Polygon_Grid_Topology
%%
%% Construct a full rectangular grid of equal square cells over a metric
%% field polygon. Every square is retained for simulation. Polygon-cell
%% intersections are stored separately for field-scale accounting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Build_Polygon_Grid_Topology:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Build_Polygon_Grid_Topology:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end

vertices = field_polygon.Vertices;
finite_vertex = all(isfinite(vertices),2);
separator = all(isnan(vertices),2);
if any(~finite_vertex & ~separator)
    error('Build_Polygon_Grid_Topology:BadField', ...
        'field_polygon vertices must be finite metric coordinates.');
end
vertices = vertices(finite_vertex,:);

x_min = min(vertices(:,1));
x_max = max(vertices(:,1));
y_min = min(vertices(:,2));
y_max = max(vertices(:,2));
nx = Grid_Cell_Count(x_max-x_min,cell_size_m);
ny = Grid_Cell_Count(y_max-y_min,cell_size_m);
N_cells = nx*ny;

cell_id = (1:N_cells)';
cell_name = cell(N_cells,1);
cell_row = zeros(N_cells,1);
cell_col = zeros(N_cells,1);
cell_x_lo = zeros(N_cells,1);
cell_x_hi = zeros(N_cells,1);
cell_y_lo = zeros(N_cells,1);
cell_y_hi = zeros(N_cells,1);
cell_overlap_area_m2 = zeros(N_cells,1);

for iy = 1:ny
    for ix = 1:nx
        c = (iy-1)*nx+ix;
        cell_name{c} = sprintf('Cell_Grid_%04d',c);
        cell_row(c) = iy;
        cell_col(c) = ix;
        cell_x_lo(c) = x_min+(ix-1)*cell_size_m;
        cell_x_hi(c) = cell_x_lo(c)+cell_size_m;
        cell_y_lo(c) = y_min+(iy-1)*cell_size_m;
        cell_y_hi(c) = cell_y_lo(c)+cell_size_m;
        square = polyshape( ...
            [cell_x_lo(c) cell_x_hi(c) cell_x_hi(c) cell_x_lo(c)], ...
            [cell_y_lo(c) cell_y_lo(c) cell_y_hi(c) cell_y_hi(c)]);
        cell_overlap_area_m2(c) = area(intersect(field_polygon,square));
    end
end

cell_area_m2 = repmat(cell_size_m^2,N_cells,1);
area_tol = 1e-9*cell_size_m^2;
cell_overlap_area_m2(cell_overlap_area_m2 < area_tol) = 0;
full_cell = abs(cell_overlap_area_m2-cell_area_m2) < area_tol;
cell_overlap_area_m2(full_cell) = cell_area_m2(full_cell);
cell_overlap_fraction = cell_overlap_area_m2./cell_area_m2;

N_edges = (nx-1)*ny+nx*(ny-1);
edge_a = zeros(N_edges,1);
edge_b = zeros(N_edges,1);
e = 0;
for iy = 1:ny
    for ix = 1:nx-1
        e = e+1;
        edge_a(e) = (iy-1)*nx+ix;
        edge_b(e) = edge_a(e)+1;
    end
end
for iy = 1:ny-1
    for ix = 1:nx
        e = e+1;
        edge_a(e) = (iy-1)*nx+ix;
        edge_b(e) = edge_a(e)+nx;
    end
end

edge_face_width_m = repmat(cell_size_m,N_edges,1);
edge_face_distance_m = repmat(cell_size_m,N_edges,1);
edge_face_distance_a_m = repmat(cell_size_m/2,N_edges,1);
edge_face_distance_b_m = repmat(cell_size_m/2,N_edges,1);
edge_active = ones(N_edges,1);
adj = sparse([edge_a;edge_b],[edge_b;edge_a],1,N_cells,N_cells) > 0;

outer_face_width_m = zeros(N_cells,1);
for c = 1:N_cells
    if cell_col(c) == 1;  outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m; end
    if cell_col(c) == nx; outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m; end
    if cell_row(c) == 1;  outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m; end
    if cell_row(c) == ny; outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m; end
end

field_area_m2 = area(field_polygon);
overlap_tol = max(1e-6,1e-9*field_area_m2);
assert(abs(sum(cell_overlap_area_m2)-field_area_m2) < overlap_tol, ...
    'Build_Polygon_Grid_Topology:AreaMismatch', ...
    'Polygon-cell intersections do not recover the field area.');
assert(abs(sum(outer_face_width_m)-2*(nx+ny)*cell_size_m) < area_tol, ...
    'Build_Polygon_Grid_Topology:PerimeterMismatch', ...
    'Grid outer-face widths do not recover the rectangular perimeter.');

[field_boundary_x_m,field_boundary_y_m] = boundary(field_polygon);
topo.field_side_m = [];
topo.buffer_width_m = 0;
topo.cell_size_m = cell_size_m;
topo.main_grid_nx = nx;
topo.main_grid_ny = ny;
topo.N_cells = N_cells;
topo.cell_id = cell_id;
topo.cell_name = cell_name;
topo.cell_row = cell_row;
topo.cell_col = cell_col;
topo.cell_region = repmat({'main'},N_cells,1);
topo.cell_subtype = repmat({'main'},N_cells,1);
topo.cell_x_lo = cell_x_lo;
topo.cell_x_hi = cell_x_hi;
topo.cell_y_lo = cell_y_lo;
topo.cell_y_hi = cell_y_hi;
topo.cell_xc = 0.5*(cell_x_lo+cell_x_hi);
topo.cell_yc = 0.5*(cell_y_lo+cell_y_hi);
topo.cell_dx = repmat(cell_size_m,N_cells,1);
topo.cell_dy = repmat(cell_size_m,N_cells,1);
topo.cell_area_m2 = cell_area_m2;
topo.cell_overlap_area_m2 = cell_overlap_area_m2;
topo.cell_overlap_fraction = cell_overlap_fraction;
topo.field_area_m2 = field_area_m2;
topo.simulated_area_m2 = sum(cell_area_m2);
topo.field_boundary_x_m = field_boundary_x_m;
topo.field_boundary_y_m = field_boundary_y_m;
topo.field_region_count = numel(regions(field_polygon));
topo.N_edges = N_edges;
topo.edge_a = edge_a;
topo.edge_b = edge_b;
topo.edge_face_width_m = edge_face_width_m;
topo.edge_face_distance_m = edge_face_distance_m;
topo.edge_face_distance_a_m = edge_face_distance_a_m;
topo.edge_face_distance_b_m = edge_face_distance_b_m;
topo.edge_active = edge_active;
topo.adj = adj;
topo.outer_face_width_m = outer_face_width_m;

end
