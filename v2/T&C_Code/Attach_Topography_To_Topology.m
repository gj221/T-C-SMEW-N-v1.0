function topo = Attach_Topography_To_Topology(topo, surface, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Attach_Topography_To_Topology
%%
%% Add surface elevation and derived topographic attributes to a
%% Build_*_Topology struct so the 2-D driver can resolve elevation-driven
%% lateral and vertical flow.  Without this call every topology is flat
%% (cell_z_m = 0, cell_slope = 0, closed perimeter), which is exactly the
%% behaviour of the code before topography was added.
%%
%% INPUTS
%%   topo      struct from Build_Cell_Topology / Build_Polygon_Grid_Topology /
%%             Build_Projected_Field_Topology / Build_Single_Cell_Topology.
%%   surface   one of
%%     (a) struct with .cell_z_m [N_cells x 1]  -- surface elevation at each
%%         cell centre [m].  Optional companion per-cell fields override the
%%         values that would otherwise be derived from a plane fit:
%%         .cell_slope [-], .cell_aspect_deg [deg from N], .cell_curv_prof
%%         [1/m], .cell_curv_plan [1/m], .cell_twi [-], .cell_svf [-].
%%     (b) struct with a gridded elevation model:
%%         .Z  [ny x nx] elevation [m], NaN outside the mapped area
%%         .x  [1 x nx]  projected easting  of column centres [m]
%%         .y  [1 x ny]  projected northing of row centres    [m]
%%         (.y descending is accepted and handled).  Slope, aspect,
%%         curvature, TWI and sky-view factor are computed on Z and then
%%         sampled at the cell centres.
%%   options   name/value pairs
%%     'OutletCells'   explicit [k x 1] cell ids to treat as domain outlets.
%%                     Default: derived from the elevation field (perimeter
%%                     cells with no lower graph neighbour).
%%     'MFDExponent'   Quinn (1991) multiple-flow-direction exponent for the
%%                     overland edge weights (default 1.1).
%%     'Source'        free-text provenance string stored on the topology.
%%
%% OUTPUT (fields added to topo)
%%   cell_z_m, cell_slope, cell_aspect_deg, cell_curv_prof, cell_curv_plan,
%%   cell_twi, cell_svf, cell_Asur = 1./cos(atan(cell_slope))
%%   edge_dz_m                    z(edge_a) - z(edge_b) [m]
%%   edge_slope                   |edge_dz_m| / edge_face_distance_m  [-]
%%   edge_face_distance_surface_m, _a_m, _b_m   slope (down-dip) distances [m]
%%   edge_overland_weight_ab, _ba  fraction of a cell's ponded outflow that
%%                                 leaves along this edge toward b (resp. a);
%%                                 the weights out of any one cell sum to <= 1,
%%                                 the deficit draining through the perimeter.
%%   cell_is_outlet               logical [N_cells x 1]
%%   topography_attached          true
%%   topography_source            provenance string
%%
%% Coupling note
%%   The elevation-head term this enables in Darcy_Lateral_Topology is
%%   layer-independent: layer k sits at the same depth below the surface in
%%   both cells of an edge, so the elevation-head difference is z_a - z_b for
%%   every layer.  A flat surface therefore reproduces the pre-topography
%%   result bit for bit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('OutletCells',[],'MFDExponent',1.1,'Source','');
if mod(numel(varargin),2) ~= 0
    error('Attach_Topography_To_Topology:BadOptions', ...
        'Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Attach_Topography_To_Topology:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~(isnumeric(opt.MFDExponent) && isscalar(opt.MFDExponent) && ...
        isfinite(opt.MFDExponent) && opt.MFDExponent > 0)
    error('Attach_Topography_To_Topology:BadOptions', ...
        'MFDExponent must be a positive finite scalar.');
end
if ~isstruct(topo) || ~isscalar(topo) || ~isfield(topo,'N_cells')
    error('Attach_Topography_To_Topology:BadTopology', ...
        'topo must be a scalar Build_*_Topology struct.');
end
required_topo = {'N_cells','N_edges','edge_a','edge_b','cell_xc','cell_yc', ...
    'edge_face_distance_m','outer_face_width_m','adj'};
missing = required_topo(~isfield(topo,required_topo));
if ~isempty(missing)
    error('Attach_Topography_To_Topology:BadTopology', ...
        'topo is missing fields: %s', strjoin(missing,', '));
end
N_cells = topo.N_cells;
N_edges = topo.N_edges;
if ~isstruct(surface) || ~isscalar(surface)
    error('Attach_Topography_To_Topology:BadSurface', ...
        'surface must be a scalar struct (see help).');
end

%% --- cell centre coordinates in the DEM/projected frame --------------
if isfield(topo,'cell_xc_projected_m') && isfield(topo,'cell_yc_projected_m')
    xc = topo.cell_xc_projected_m(:);
    yc = topo.cell_yc_projected_m(:);
    frame = 'projected';
else
    xc = topo.cell_xc(:);
    yc = topo.cell_yc(:);
    frame = 'local';
end

attr = struct();
if isfield(surface,'Z')
    %% (b) gridded elevation model -------------------------------------
    grid = Field_Topographic_Attributes(surface);
    x = grid.x; y = grid.y;
    sample = @(F) sample_grid(x,y,F,xc,yc);
    attr.cell_z_m        = sample(grid.Z);
    attr.cell_slope      = max(0,sample(grid.slope));
    attr.cell_aspect_deg = sample(grid.aspect_deg);
    attr.cell_curv_prof  = sample(grid.curv_prof);
    attr.cell_curv_plan  = sample(grid.curv_plan);
    attr.cell_twi        = sample(grid.twi);
    attr.cell_svf        = min(1,max(0,sample(grid.svf)));
    if any(~isfinite(attr.cell_z_m))
        error('Attach_Topography_To_Topology:CellOutsideDEM', ...
            'One or more cell centres fall outside the supplied elevation grid.');
    end
    source = ternary(isempty(opt.Source), ...
        sprintf('grid DEM %dx%d (%.2f m), %s frame', ...
        size(grid.Z,1),size(grid.Z,2),grid.cellsize,frame), opt.Source);
elseif isfield(surface,'cell_z_m')
    %% (a) direct per-cell injection ---------------------------------
    z = surface.cell_z_m(:);
    if numel(z) ~= N_cells || any(~isfinite(z))
        error('Attach_Topography_To_Topology:BadSurface', ...
            'surface.cell_z_m must hold one finite elevation per cell.');
    end
    attr.cell_z_m = z;
    attr = plane_fit_attributes(attr,topo);   % slope/aspect/curv from the graph
    passthrough = {'cell_slope','cell_aspect_deg','cell_curv_prof', ...
        'cell_curv_plan','cell_twi','cell_svf'};
    for f = passthrough
        if isfield(surface,f{1})
            v = surface.(f{1})(:);
            if numel(v) ~= N_cells || any(~isfinite(v))
                error('Attach_Topography_To_Topology:BadSurface', ...
                    'surface.%s must hold one finite value per cell.',f{1});
            end
            attr.(f{1}) = v;
        end
    end
    source = ternary(isempty(opt.Source),'per-cell elevation injection',opt.Source);
else
    error('Attach_Topography_To_Topology:BadSurface', ...
        'surface must contain either .Z (gridded) or .cell_z_m (per cell).');
end

%% --- defaults for any attribute not supplied ------------------------
attr = default_attr(attr,'cell_slope',zeros(N_cells,1));
attr = default_attr(attr,'cell_aspect_deg',zeros(N_cells,1));
attr = default_attr(attr,'cell_curv_prof',zeros(N_cells,1));
attr = default_attr(attr,'cell_curv_plan',zeros(N_cells,1));
attr = default_attr(attr,'cell_twi',nan(N_cells,1));
attr = default_attr(attr,'cell_svf',ones(N_cells,1));
attr.cell_slope = max(0,attr.cell_slope);
attr.cell_svf = min(1,max(0,attr.cell_svf));
attr.cell_Asur = 1./cos(atan(attr.cell_slope));

%% --- edge geometry -------------------------------------------------
z = attr.cell_z_m;
ea = topo.edge_a(:); eb = topo.edge_b(:);
edge_dz_m = z(ea) - z(eb);
plan_dist = topo.edge_face_distance_m(:);
edge_slope = abs(edge_dz_m)./max(plan_dist,eps);
surf_factor = 1./cos(atan(edge_slope));
edge_face_distance_surface_m = plan_dist.*surf_factor;
if isfield(topo,'edge_face_distance_a_m') && isfield(topo,'edge_face_distance_b_m')
    edge_face_distance_surface_a_m = topo.edge_face_distance_a_m(:).*surf_factor;
    edge_face_distance_surface_b_m = topo.edge_face_distance_b_m(:).*surf_factor;
else
    edge_face_distance_surface_a_m = 0.5*edge_face_distance_surface_m;
    edge_face_distance_surface_b_m = 0.5*edge_face_distance_surface_m;
end

%% --- overland multiple-flow-direction weights (Quinn 1991) ---------
[w_ab,w_ba] = overland_weights(topo,z,edge_face_distance_surface_m, ...
    opt.MFDExponent);

%% --- outlet cells -------------------------------------------------
if ~isempty(opt.OutletCells)
    ids = opt.OutletCells(:);
    if any(~ismember(ids,1:N_cells))
        error('Attach_Topography_To_Topology:BadOutlet', ...
            'OutletCells must be valid cell ids.');
    end
    is_outlet = false(N_cells,1); is_outlet(ids) = true;
else
    is_outlet = derive_outlets(topo,z);
end

%% --- pack --------------------------------------------------------
topo.cell_z_m           = attr.cell_z_m;
topo.cell_slope         = attr.cell_slope;
topo.cell_aspect_deg    = attr.cell_aspect_deg;
topo.cell_curv_prof     = attr.cell_curv_prof;
topo.cell_curv_plan     = attr.cell_curv_plan;
topo.cell_twi           = attr.cell_twi;
topo.cell_svf           = attr.cell_svf;
topo.cell_Asur          = attr.cell_Asur;
topo.edge_dz_m                    = edge_dz_m;
topo.edge_slope                   = edge_slope;
topo.edge_face_distance_surface_m = edge_face_distance_surface_m;
topo.edge_face_distance_surface_a_m = edge_face_distance_surface_a_m;
topo.edge_face_distance_surface_b_m = edge_face_distance_surface_b_m;
topo.edge_overland_weight_ab = w_ab;
topo.edge_overland_weight_ba = w_ba;
topo.cell_is_outlet     = is_outlet;
topo.topography_attached = true;
topo.topography_source   = source;
if exist('grid','var'); topo.topography_grid = grid; end

end

%% =======================================================================
function v = sample_grid(x,y,F,xq,yq)
G = griddedInterpolant({y,x},F,'linear','none');
v = G(yq(:),xq(:));
%% Fall back to nearest for query points a sliver outside the grid.
bad = ~isfinite(v);
if any(bad)
    Gn = griddedInterpolant({y,x},F,'nearest','nearest');
    v(bad) = Gn(yq(bad),xq(bad));
end
end

%% =======================================================================
function attr = plane_fit_attributes(attr,topo)
%% Least-squares plane through the cell centres -> a single slope/aspect,
%% assigned to every cell.  Used when only elevations are injected (tests,
%% idealised hillslopes); a real run supplies a gridded DEM instead.
x = topo.cell_xc(:); y = topo.cell_yc(:); z = attr.cell_z_m;
nz = numel(z);
if nz < 3
    attr.cell_slope = zeros(nz,1);
    attr.cell_aspect_deg = zeros(nz,1);
    attr.cell_curv_prof = zeros(nz,1);
    attr.cell_curv_plan = zeros(nz,1);
    return
end
xr = x-mean(x); yr = y-mean(y);
%% Only fit a gradient in coordinates that actually vary (a 1-D chain of
%% cells has no spread in one axis -> that column would be rank deficient).
cols = {}; picks = [];
if std(xr) > 1e-9*max(abs(xr)+1); cols{end+1} = xr; picks(end+1) = 1; end
if std(yr) > 1e-9*max(abs(yr)+1); cols{end+1} = yr; picks(end+1) = 2; end
gx = 0; gy = 0;
if ~isempty(cols)
    G = [cell2mat(cols), ones(nz,1)];
    c = G\(z-mean(z));
    for k = 1:numel(picks)
        if picks(k) == 1; gx = c(k); else; gy = c(k); end
    end
    resid = (z-mean(z)) - G*c;
else
    resid = z-mean(z);
end
slope = hypot(gx,gy);
aspect = mod(atan2d(-gx,-gy),360);          % steepest-descent direction, from N
attr.cell_slope = repmat(slope,nz,1);
attr.cell_aspect_deg = repmat(aspect,nz,1);
scale = max(hypot(xr,yr),eps);
attr.cell_curv_prof = -2*resid./scale.^2;    % convex-up positive
attr.cell_curv_plan = zeros(nz,1);
end

%% =======================================================================
function [w_ab,w_ba] = overland_weights(topo,z,surf_dist,p)
%% Multiple-flow-direction partition of ponded outflow on the cell graph.
%% w_ab(e): fraction of cell edge_a(e)'s ponded outflow that leaves along
%% edge e toward edge_b(e).  Weights out of a cell sum to <= 1; the rest
%% leaves through the perimeter (handled by the caller / outlet mask).
N_edges = topo.N_edges;
w_ab = zeros(N_edges,1);
w_ba = zeros(N_edges,1);
gain = zeros(topo.N_cells,1);               % sum of downslope (drop/dist)^p
for e = 1:N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    d = surf_dist(e);
    s_ab = (z(a)-z(b))/d;                   % positive: a drains to b
    if s_ab > 0
        w_ab(e) = s_ab^p;
        gain(a) = gain(a) + w_ab(e);
    else
        s_ba = -s_ab;
        w_ba(e) = s_ba^p;
        gain(b) = gain(b) + w_ba(e);
    end
end
for e = 1:N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    if w_ab(e) > 0 && gain(a) > 0; w_ab(e) = w_ab(e)/gain(a); end
    if w_ba(e) > 0 && gain(b) > 0; w_ba(e) = w_ba(e)/gain(b); end
end
end

%% =======================================================================
function is_outlet = derive_outlets(topo,z)
%% A cell is an outlet if it borders the domain perimeter and no graph
%% neighbour is lower -- ponded water there leaves through the perimeter.
N_cells = topo.N_cells;
is_outlet = false(N_cells,1);
on_perimeter = topo.outer_face_width_m(:) > 0;
for c = find(on_perimeter).'
    nb = [topo.edge_b(topo.edge_a==c); topo.edge_a(topo.edge_b==c)];
    if isempty(nb) || all(z(nb) >= z(c) - 1e-9)
        is_outlet(c) = true;
    end
end
if ~any(is_outlet) && any(on_perimeter)
    %% Degenerate flat perimeter: the single lowest perimeter cell drains.
    perim = find(on_perimeter);
    [~,lo] = min(z(perim));
    is_outlet(perim(lo)) = true;
end
end

%% =======================================================================
function attr = default_attr(attr,name,value)
if ~isfield(attr,name) || isempty(attr.(name))
    attr.(name) = value;
end
end

%% =======================================================================
function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end
