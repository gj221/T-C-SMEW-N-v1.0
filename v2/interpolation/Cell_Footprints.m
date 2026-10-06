function [footprints,field] = Cell_Footprints(topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell_Footprints
%%
%% Each grid cell's footprint in the grid's local metric frame: its square
%% clipped to the field polygon, i.e. the ground its column represents.
%% Per-cell inputs drawn from maps (land cover, soil map units) are area
%% weights over these footprints.  field is the field polygon in the same
%% frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
need = {'N_cells','cell_x_lo','cell_x_hi','cell_y_lo','cell_y_hi', ...
    'field_boundary_x_m','field_boundary_y_m'};
if ~isstruct(topo) || ~all(isfield(topo,need))
    error('Cell_Footprints:BadTopology', ...
        'The topology must carry its cell squares and field boundary.');
end
warn = warning('off','MATLAB:polyshape:repairedBySimplify');
restore = onCleanup(@() warning(warn));
field = polyshape(topo.field_boundary_x_m,topo.field_boundary_y_m);
footprints = repmat(polyshape(),topo.N_cells,1);
for c = 1:topo.N_cells
    square = polyshape([topo.cell_x_lo(c) topo.cell_x_hi(c) topo.cell_x_hi(c) topo.cell_x_lo(c)], ...
        [topo.cell_y_lo(c) topo.cell_y_lo(c) topo.cell_y_hi(c) topo.cell_y_hi(c)]);
    footprints(c) = intersect(field,square);
end
end
