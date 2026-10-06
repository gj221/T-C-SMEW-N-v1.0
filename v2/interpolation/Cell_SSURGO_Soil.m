function soil = Cell_SSURGO_Soil(topo,mapunits,horizons,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell_SSURGO_Soil
%%
%% Per-cell soil from SSURGO for a field grid.
%%   mapunits  table: mukey, musym, wkt -- map-unit polygons as WGS84
%%             longitude/latitude WKT (POLYGON or MULTIPOLYGON), e.g. from
%%             USDA Soil Data Access
%%   horizons  table of the map units' component horizons (mukey, cokey,
%%             compname, comppct_r, hzname, hzdept_r, hzdepb_r,
%%             sandtotal_r, claytotal_r, om_r, dbthirdbar_r, ph1to1h2o_r,
%%             cec7_r, ecec_r)
%% Options
%%   'Restrictions'  table mukey, compname, resdept_r (depth to a
%%                   restrictive layer, cm)
%%   'Depth'         [top bottom] cm over which mineral horizons are
%%                   depth-weighted (default [0 30]); organic (O) horizons
%%                   and horizons without texture (bedrock) are skipped
%%   'Component'     '' (default): each cell takes the dominant major
%%                   component of the map unit covering most of its
%%                   footprint.  A component name ('Cheshire') gives every
%%                   cell that component's soil, still recording each
%%                   cell's own map unit.
%%
%% Each cell's map-unit shares are exact areas of its footprint
%% (Cell_Footprints) in the grid's local frame.  Within a complex map unit
%% SSURGO does not locate the components, so a cell takes one component
%% whole rather than a blend of soils that do not exist together.
%% CEC is CEC7 where SSURGO gives it for the component, ECEC otherwise
%% (CEC_basis says which).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
p.addParameter('Restrictions',table());
p.addParameter('Depth',[0 30]);
p.addParameter('Component','');
p.parse(varargin{:});
depth = p.Results.Depth;
restr = p.Results.Restrictions;
fixed = char(string(p.Results.Component));
if istable(mapunits) == 0; mapunits = readtable(mapunits,'TextType','string'); end
if istable(horizons) == 0; horizons = readtable(horizons,'TextType','string'); end
mapunits.mukey = string(mapunits.mukey); mapunits.musym = string(mapunits.musym);
horizons.mukey = string(horizons.mukey); horizons.compname = string(horizons.compname);

%% map-unit polygons in the grid's local frame, merged per map unit
utm = projcrs(topo.projected_crs_epsg);
keys = unique(mapunits.mukey,'stable');
warn = warning('off','MATLAB:polyshape:repairedBySimplify');
restore = onCleanup(@() warning(warn));
mu = repmat(polyshape(),numel(keys),1);
sym = strings(numel(keys),1);
for k = 1:numel(keys)
    rows = find(mapunits.mukey == keys(k));
    sym(k) = mapunits.musym(rows(1));
    for r = rows'
        [lon,lat] = parse_wkt(char(mapunits.wkt(r)));
        [ux,uy] = projfwd(utm,lat,lon);
        [lx,ly] = Grid_Projected_To_Local(ux,uy,topo);
        mu(k) = union(mu(k),polyshape(lx,ly));
    end
end

%% each cell's map-unit shares
fp = Cell_Footprints(topo);
n = numel(fp);
share = zeros(n,numel(keys));
for c = 1:n
    a = area(fp(c));
    for k = 1:numel(keys)
        share(c,k) = area(intersect(fp(c),mu(k)))/a;
    end
end
if any(abs(sum(share,2)-1) > 1e-3)
    error('Cell_SSURGO_Soil:IncompleteCoverage', ...
        'The map-unit polygons cover only %.4f of a cell footprint.',min(sum(share,2)));
end
[mshare,dom] = max(share,[],2);

%% component properties
soil = table(topo.cell_id(:),keys(dom),sym(dom),mshare, ...
    'VariableNames',{'cell_id','mukey','musym','mapunit_share'});
vars = {'component','comppct_r','sand_fraction','clay_fraction','OM_fraction', ...
    'bulk_density_g_cm3','pH_1to1_water','CEC_cmolc_kg','CEC_basis','restriction_depth_cm'};
props = cell(n,numel(vars));
cache = containers.Map();
for c = 1:n
    key = soil.mukey(c) + "|" + string(fixed);
    if ~isKey(cache,key)
        cache(key) = component_properties(horizons,restr,soil.mukey(c),fixed,depth);
    end
    props(c,:) = cache(key);
end
soil = [soil cell2table(props,'VariableNames',vars)];
for k = 1:numel(keys)
    soil.("share_" + matlab.lang.makeValidName(sym(k))) = share(:,k);
end
end

function out = component_properties(H,restr,mukey,fixed,depth)
if isempty(fixed)
    C = H(H.mukey == mukey,:);
    if ismember('majcompflag',C.Properties.VariableNames)
        C = C(C.majcompflag == "Yes",:);
    end
    [~,i] = max(C.comppct_r);
    name = C.compname(i);
    C = C(C.compname == name & C.cokey == C.cokey(i),:);
else
    C = H(H.compname == string(fixed),:);
    if isempty(C)
        error('Cell_SSURGO_Soil:UnknownComponent','No component named %s.',fixed);
    end
    C = C(C.cokey == C.cokey(1),:);
    name = string(fixed);
end
top = C.hzdept_r; bot = C.hzdepb_r;
w = max(0,min(bot,depth(2)) - max(top,depth(1)));
mineral = ~startsWith(string(C.hzname),"O") & isfinite(C.sandtotal_r) & ...
    isfinite(C.claytotal_r) & w > 0;
if ~any(mineral)
    error('Cell_SSURGO_Soil:NoMineralHorizon', ...
        '%s has no mineral horizon with texture in %g-%g cm.',name,depth(1),depth(2));
end
wm = @(v) sum(v(mineral & isfinite(v)).*w(mineral & isfinite(v))) / ...
    sum(w(mineral & isfinite(v)));
cec = wm(C.cec7_r); basis = "CEC7";
if ~isfinite(cec); cec = wm(C.ecec_r); basis = "ECEC"; end
rdepth = NaN;
if ~isempty(restr)
    m = string(restr.mukey) == mukey & string(restr.compname) == name;
    if any(m); rdepth = min(restr.resdept_r(m)); end
end
out = {name,C.comppct_r(1),wm(C.sandtotal_r)/100,wm(C.claytotal_r)/100, ...
    wm(C.om_r)/100,wm(C.dbthirdbar_r),wm(C.ph1to1h2o_r),cec,basis,rdepth};
end

function [x,y] = parse_wkt(wkt)
%% POLYGON / MULTIPOLYGON rings as NaN-separated vertex lists.
body = regexprep(wkt,'^\s*(MULTI)?POLYGON\s*','');
rings = regexp(body,'\(([^()]+)\)','tokens');
x = []; y = [];
for k = 1:numel(rings)
    v = sscanf(strrep(rings{k}{1},',',' '),'%f');
    v = reshape(v,2,[]);
    x = [x NaN v(1,:)]; %#ok<AGROW>
    y = [y NaN v(2,:)]; %#ok<AGROW>
end
x = x(2:end); y = y(2:end);
end
