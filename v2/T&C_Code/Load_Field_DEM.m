function dem = Load_Field_DEM(source, topo_or_crs, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load_Field_DEM
%%
%% Read one or more projected elevation GeoTIFFs, mosaic them onto a common
%% grid, resample to a chosen spacing, clip to the field, and return a
%% plain gridded elevation model for Field_Topographic_Attributes /
%% Attach_Topography_To_Topology.
%%
%% The DEM must be supplied in a PROJECTED CRS (UTM).  It should match the
%% field topology's CRS; a same-zone datum difference (NAD83 vs WGS84,
%% ~1 m in CONUS) is tolerated with a warning.  Export the tile from the
%% DEM portal in UTM zone 17N for the Bruce Farm / Adam Wilbourne site.
%%
%% INPUTS
%%   source        a GeoTIFF path, a directory of GeoTIFF tiles
%%                 (*.tif / *.tiff), or a {cell} of paths.
%%   topo_or_crs   the field topology struct (uses .projected_crs_epsg and
%%                 the projected cell-centre bounds), an EPSG code, a
%%                 projcrs, or [] to keep the source extent and CRS.
%%
%% Name/value options
%%   'OutResolutionM'  output spacing [m] (default: native, snapped to
%%                     {0.5 1 2 3 5 10 30}).
%%   'BufferM'         margin kept around the field [m] (default 100).
%%   'BBoxLatLon'      [lonW lonE latS latN] clip window (WGS84 degrees).
%%   'FillMaxCells'    largest NaN void filled by interpolation (default 50).
%%
%% OUTPUT
%%   dem   struct: .Z [ny x nx] elevation [m]; .x [1 x nx], .y [1 x ny]
%%         ascending projected coordinates [m]; .cellsize_m;
%%         .native_resolution_m; .crs_name; .source_files.
%%
%% Requires the Mapping Toolbox (readgeoraster).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('OutResolutionM',[],'BufferM',100,'BBoxLatLon',[], ...
    'FillMaxCells',50);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if nargin < 2; topo_or_crs = []; end

files = resolve_source_files(source);
if isempty(files)
    error('Load_Field_DEM:NoFiles','No GeoTIFF found at %s.',string(source));
end

%% --- read tiles onto ascending-(x,y) grids -----------------------
tiles = cell(numel(files),1);
src_crs = [];
for k = 1:numel(files)
    [Z,R] = readgeoraster(files{k},'OutputType','double');
    Z(Z <= -9990 | Z >= 9000) = NaN;              % common NoData sentinels
    if ~(isprop(R,'ProjectedCRS') && ~isempty(R.ProjectedCRS))
        error('Load_Field_DEM:GeographicSource', ...
            ['%s is not in a projected CRS. Re-export the DEM from the ' ...
             'portal in a projected (UTM) CRS.'], files{k});
    end
    if isempty(src_crs); src_crs = R.ProjectedCRS; end
    nx = R.RasterSize(2); ny = R.RasterSize(1);
    %% Ascending pixel-centre coordinates.  Z rows run north->south when
    %% ColumnsStartFrom is 'north', so flip Z (not ys) to match ascending y.
    xs = R.XWorldLimits(1) + ((1:nx)-0.5)*R.CellExtentInWorldX;
    ys = R.YWorldLimits(1) + ((1:ny)-0.5)*R.CellExtentInWorldY;
    if ~isprop(R,'ColumnsStartFrom') || strcmp(R.ColumnsStartFrom,'north')
        Z = flipud(Z);
    end
    if isprop(R,'RowsStartFrom') && strcmp(R.RowsStartFrom,'east')
        Z = fliplr(Z);
    end
    tiles{k} = struct('Z',Z,'x',xs,'y',ys, ...
        'res',mean([R.CellExtentInWorldX R.CellExtentInWorldY]));
end
native_res = tiles{1}.res;
out_res = opt.OutResolutionM;
if isempty(out_res); out_res = snap_resolution(native_res); end

%% --- target CRS + clip window (in the TARGET CRS) --------------
[tgt_crs,tgt_epsg,bounds_xy] = resolve_target(topo_or_crs,src_crs,opt);
reproject = ~same_crs(src_crs,tgt_crs);

xt = (bounds_xy(1)+out_res/2) : out_res : bounds_xy(2);
yt = (bounds_xy(3)+out_res/2) : out_res : bounds_xy(4);
if numel(xt) < 3 || numel(yt) < 3
    error('Load_Field_DEM:TinyExtent', ...
        'Clip window is %g x %g m; check BBoxLatLon / topology bounds.', ...
        diff(bounds_xy(1:2)),diff(bounds_xy(3:4)));
end
[XT,YT] = meshgrid(xt,yt);

%% target grid points -> source-CRS coordinates for sampling
if reproject
    [latq,lonq] = projinv(tgt_crs,XT,YT);
    [XS,YS]     = projfwd(src_crs,latq,lonq);
else
    XS = XT; YS = YT;
end

%% --- sample tiles, average overlaps ----------------------------
acc = zeros(size(XT)); cnt = zeros(size(XT));
for k = 1:numel(tiles)
    G = griddedInterpolant({tiles{k}.y,tiles{k}.x},tiles{k}.Z,'linear','none');
    v = G(YS,XS);
    ok = isfinite(v);
    acc(ok) = acc(ok) + v(ok);
    cnt(ok) = cnt(ok) + 1;
end
Z = acc./cnt;
Z(cnt == 0) = NaN;
if all(isnan(Z(:)))
    error('Load_Field_DEM:NoOverlap', ...
        ['DEM tiles do not cover the requested area. Tile 1 spans ' ...
         'x=[%.0f %.0f] y=[%.0f %.0f] in its own CRS; the field/box maps ' ...
         'to source x=[%.0f %.0f] y=[%.0f %.0f]. If those disagree the DEM ' ...
         'does not cover this field -- fetch the tile(s) that do.'], ...
        min(tiles{1}.x),max(tiles{1}.x),min(tiles{1}.y),max(tiles{1}.y), ...
        min(XS(:)),max(XS(:)),min(YS(:)),max(YS(:)));
end
Z = fill_voids(Z,opt.FillMaxCells);

dem = struct('Z',Z,'x',xt,'y',yt,'cellsize_m',out_res, ...
    'native_resolution_m',native_res,'epsg',tgt_epsg, ...
    'crs_name',char(tgt_crs.Name),'source_crs_name',char(src_crs.Name), ...
    'reprojected',reproject,'source_files',{files});
fprintf(['Load_Field_DEM: %d tile(s) -> %d x %d at %.2f m, %.1f%% NaN ' ...
    '[%s%s]\n'], numel(files),size(Z,1),size(Z,2),out_res, ...
    100*mean(isnan(Z(:))),char(tgt_crs.Name), ...
    ternary(reproject,sprintf(' <- %s',char(src_crs.Name)),''));
end

%% =======================================================================
function [tgt_crs,tgt_epsg,bounds_xy] = resolve_target(topo_or_crs,src_crs,opt)
tgt_epsg = [];
if isempty(topo_or_crs)
    tgt_crs = src_crs;
elseif isa(topo_or_crs,'projcrs')
    tgt_crs = topo_or_crs;
elseif isnumeric(topo_or_crs) && isscalar(topo_or_crs)
    tgt_epsg = topo_or_crs; tgt_crs = projcrs(tgt_epsg);
elseif isstruct(topo_or_crs) && isfield(topo_or_crs,'projected_crs_epsg')
    tgt_epsg = topo_or_crs.projected_crs_epsg; tgt_crs = projcrs(tgt_epsg);
else
    error('Load_Field_DEM:BadTarget', ...
        'topo_or_crs must be a topology struct, EPSG code, projcrs or [].');
end
if isempty(tgt_epsg); tgt_epsg = epsg_of(tgt_crs); end

if ~isempty(opt.BBoxLatLon)
    bb = opt.BBoxLatLon;
    [x,y] = projfwd(tgt_crs,[bb(3) bb(3) bb(4) bb(4)],[bb(1) bb(2) bb(2) bb(1)]);
    bounds_xy = [min(x) max(x) min(y) max(y)] + opt.BufferM*[-1 1 -1 1];
elseif isstruct(topo_or_crs) && isfield(topo_or_crs,'cell_xc_projected_m')
    x = topo_or_crs.cell_xc_projected_m(:);
    y = topo_or_crs.cell_yc_projected_m(:);
    half = 0;
    if isfield(topo_or_crs,'cell_size_m') && ~isempty(topo_or_crs.cell_size_m)
        half = topo_or_crs.cell_size_m;
    end
    bounds_xy = [min(x) max(x) min(y) max(y)] + (opt.BufferM+half)*[-1 1 -1 1];
else
    error('Load_Field_DEM:NoExtent', ...
        'Pass a topology struct with projected cell centres or a BBoxLatLon.');
end
end

%% =======================================================================
function tf = same_crs(a,b)
try
    tf = strcmp(char(wktstring(a)),char(wktstring(b)));
catch
    tf = strcmp(char(a.Name),char(b.Name));
end
end

%% =======================================================================
function e = epsg_of(crs)
e = NaN;
try
    m = regexp(char(wktstring(crs)),'ID\["EPSG",(\d+)\]','tokens');
    if ~isempty(m); e = str2double(m{end}{1}); end
catch
end
end

%% =======================================================================
function files = resolve_source_files(source)
if iscell(source)
    files = source(:).';
elseif ischar(source) || isstring(source)
    source = char(source);
    if isfolder(source)
        d = [dir(fullfile(source,'*.tif')); dir(fullfile(source,'*.tiff'))];
        files = arrayfun(@(f) fullfile(f.folder,f.name),d,'UniformOutput',false).';
    elseif isfile(source)
        files = {source};
    else
        files = {};
    end
else
    files = {};
end
files = files(cellfun(@isfile,files));
end

%% =======================================================================
function r = snap_resolution(res)
ladder = [0.5 1 2 3 5 10 30];
[~,i] = min(abs(ladder-res));
r = ladder(i);
end

%% =======================================================================
function Zf = fill_voids(Z,max_cells)
Zf = Z;
bad = ~isfinite(Zf);
if ~any(bad(:)) || all(bad(:)); return; end
L = label_components(bad);
for lab = 1:max(L(:))
    m = L == lab;
    if nnz(m) > max_cells; continue; end
    [ii,jj] = find(m);
    for k = 1:numel(ii)
        i = ii(k); j = jj(k);
        w = Zf(max(1,i-2):min(end,i+2),max(1,j-2):min(end,j+2));
        w = w(isfinite(w));
        if ~isempty(w); Zf(i,j) = mean(w); end
    end
end
end

%% =======================================================================
function L = label_components(BW)
L = zeros(size(BW)); next = 0;
for j = 1:size(BW,2)
    for i = 1:size(BW,1)
        if BW(i,j) && L(i,j) == 0
            next = next + 1; L(i,j) = next; stack = [i j];
            while ~isempty(stack)
                p = stack(end,:); stack(end,:) = [];
                for d = [-1 0;1 0;0 -1;0 1].'
                    ni = p(1)+d(1); nj = p(2)+d(2);
                    if ni>=1 && ni<=size(BW,1) && nj>=1 && nj<=size(BW,2) ...
                            && BW(ni,nj) && L(ni,nj) == 0
                        L(ni,nj) = next; stack(end+1,:) = [ni nj]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

%% =======================================================================
function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end
