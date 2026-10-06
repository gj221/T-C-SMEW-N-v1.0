function [cover,diag] = Cell_Land_Cover(topo,mode,source,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell_Land_Cover
%%
%% Per-cell land cover for a field grid, as the cover_* columns of the
%% per-cell agronomic table (Configure_Cell_Agronomic_Inputs turns them
%% into each cell's Ccrown, Curb, Cwat, Cbare and Crock):
%%   cover_tree   crown 1: tree canopy, with the low vegetation beneath it
%%   cover_grass  crown 2: open low vegetation
%%   cover_urban, cover_water, cover_bare, cover_rock
%% Fractions of the cell's ground inside the field; each row sums to 1.
%%
%% Cell_Land_Cover(topo,'raster',raster_file,'ClassTable',T)
%%   A land-cover GeoTIFF (default classes: NLCD_Cover_Class_Table).  Each
%%   cell's footprint (Cell_Footprints) is intersected exactly with every
%%   raster pixel it touches, in the grid's local metric frame, so a cell
%%   straddling pixels takes each in proportion to the area it covers.
%%   Only the window around the field is read, so a continental raster is
%%   never loaded whole.  Pixels flagged missing must not touch the field.
%%
%% Cell_Land_Cover(topo,'manual',fractions)
%%   fractions: a table with cell_id and the six cover_* columns (one row
%%   per cell), or a one-row table / struct without cell_id applied to
%%   every cell.
%%
%% diag (raster mode): class codes, per-cell class areas and the area each
%% cell's pixels account for.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = {'cover_tree','cover_grass','cover_urban','cover_water','cover_bare','cover_rock'};
diag = struct();
switch lower(char(mode))
    case 'manual'
        F = manual_fractions(topo,source,names);
        label = "manual";
    case 'raster'
        p = inputParser;
        p.addParameter('ClassTable',NLCD_Cover_Class_Table());
        p.parse(varargin{:});
        classes = p.Results.ClassTable;
        if ~istable(classes); classes = readtable(classes); end
        [F,diag] = raster_fractions(topo,char(source),classes,names);
        [~,base,ext] = fileparts(char(source));
        label = "raster:" + string([base ext]);
    otherwise
        error('Cell_Land_Cover:BadMode','mode must be ''raster'' or ''manual''.');
end
F = close_fractions(F);
cover = array2table(F,'VariableNames',names);
cover = addvars(cover,topo.cell_id(:),'Before',1,'NewVariableNames','cell_id');
cover.cover_source = repmat(label,height(cover),1);
end

function F = manual_fractions(topo,source,names)
if isstruct(source); source = struct2table(source); end
if ~istable(source) || ~all(ismember(names,source.Properties.VariableNames))
    error('Cell_Land_Cover:BadManual', ...
        'Manual cover needs a table or struct with %s.',strjoin(names,', '));
end
if ~ismember('cell_id',source.Properties.VariableNames)
    if height(source) ~= 1
        error('Cell_Land_Cover:BadManual', ...
            'Without cell_id, manual cover must be a single row applied to every cell.');
    end
    F = repmat(source{1,names},topo.N_cells,1);
    return
end
[found,loc] = ismember(topo.cell_id(:),source.cell_id);
if ~all(found) || height(source) ~= topo.N_cells
    error('Cell_Land_Cover:CellSetMismatch', ...
        'Manual cover must give exactly one row for every topology cell.');
end
F = source{loc,names};
end

function [F,diag] = raster_fractions(topo,file,classes,names)
need = [{'class'},names];
if ~all(ismember(need,classes.Properties.VariableNames))
    error('Cell_Land_Cover:BadClassTable','The class table needs %s.',strjoin(need,', '));
end
rowsum = sum(classes{:,names},2);
if any(abs(rowsum-1) > 1e-9) || any(classes{:,names}(:) < 0)
    error('Cell_Land_Cover:BadClassTable', ...
        'Every class must divide among the cover components with shares summing to 1.');
end
I = georasterinfo(file);
R = I.RasterReference;
crs = I.CoordinateReferenceSystem;
nodata = I.MissingDataIndicator;
utm = projcrs(topo.projected_crs_epsg);

%% the window: the field's extent in the raster's grid, plus two pixels
[lat,lon] = projinv(utm,topo.field_boundary_projected_x_m,topo.field_boundary_projected_y_m);
ok = isfinite(lat) & isfinite(lon);
[r,c] = to_raster_index(R,crs,lat(ok),lon(ok));
r1 = max(1,min(r)-2); r2 = min(R.RasterSize(1),max(r)+2);
c1 = max(1,min(c)-2); c2 = min(R.RasterSize(2),max(c)+2);
A = imread(file,'PixelRegion',{[r1 r2],[c1 c2]});
A = double(A(:,:,1));

%% footprints and their bounding boxes in the grid's local frame
fp = Cell_Footprints(topo);
n = numel(fp);
bb = zeros(n,4);
fp_area = zeros(n,1);
for k = 1:n
    [xl,yl] = boundingbox(fp(k));
    bb(k,:) = [xl yl];
    fp_area(k) = area(fp(k));
end
codes = unique(A(:))';
class_area = zeros(n,numel(codes));
warn = warning('off','MATLAB:polyshape:repairedBySimplify');
restore = onCleanup(@() warning(warn));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        [px,py] = pixel_corners_local(R,crs,utm,topo,r1+i-1,c1+j-1);
        cand = find(bb(:,2) > min(px) & bb(:,1) < max(px) & ...
            bb(:,4) > min(py) & bb(:,3) < max(py));
        if isempty(cand); continue; end
        pix = polyshape(px,py);
        col = codes == A(i,j);
        for k = cand'
            a = area(intersect(fp(k),pix));
            if a > 0; class_area(k,col) = class_area(k,col) + a; end
        end
    end
end
accounted = sum(class_area,2);
if any(abs(accounted - fp_area) > 1e-6*max(1,fp_area))
    error('Cell_Land_Cover:IncompleteCoverage', ...
        'The raster window does not cover every cell footprint (worst %.3g m2 short).', ...
        max(fp_area - accounted));
end
if ~isempty(nodata) && any(codes == nodata) && any(class_area(:,codes == nodata) > 0)
    error('Cell_Land_Cover:NoDataInField', ...
        'Missing-data pixels (%g) touch the field.',nodata);
end
[known,loc] = ismember(codes,classes.class);
present = any(class_area > 0,1);
if any(present & ~known)
    error('Cell_Land_Cover:UnmappedClass', ...
        'Raster class(es) %s are not in the class table.',mat2str(codes(present & ~known)));
end
share = class_area(:,known) ./ accounted;
F = share * classes{loc(known),names};
diag = struct('class',codes,'class_area_m2',class_area,'footprint_area_m2',fp_area, ...
    'window_rows',[r1 r2],'window_cols',[c1 c2],'source',file);
end

function [r,c] = to_raster_index(R,crs,lat,lon)
if isa(crs,'projcrs')
    [x,y] = projfwd(crs,lat,lon);
    [r,c] = worldToDiscrete(R,x,y);
else
    [r,c] = geographicToDiscrete(R,lat,lon);
end
end

function [px,py] = pixel_corners_local(R,crs,utm,topo,row,col)
xi = col + [-0.5 0.5 0.5 -0.5];
yi = row + [-0.5 -0.5 0.5 0.5];
if isa(crs,'projcrs')
    [x,y] = intrinsicToWorld(R,xi,yi);
    [lat,lon] = projinv(crs,x,y);
else
    [lat,lon] = intrinsicToGeographic(R,xi,yi);
end
[ux,uy] = projfwd(utm,lat,lon);
[px,py] = Grid_Projected_To_Local(ux,uy,topo);
end

function F = close_fractions(F)
if any(~isfinite(F(:))) || any(F(:) < -1e-9) || any(F(:) > 1+1e-9)
    error('Cell_Land_Cover:BadFractions','Cover fractions must lie in [0, 1].');
end
s = sum(F,2);
if any(abs(s-1) > 1e-6)
    error('Cell_Land_Cover:BadFractions', ...
        'Cover fractions must sum to 1 in every cell (worst %.3g).',max(abs(s-1)));
end
%% Rows within rounding of one pass unchanged, so manual fractions are kept
%% exactly; a larger (still accepted) deviation is spread proportionally.
F = min(max(F,0),1);
off = abs(sum(F,2) - 1) > 1e-12;
F(off,:) = F(off,:) ./ sum(F(off,:),2);
end
