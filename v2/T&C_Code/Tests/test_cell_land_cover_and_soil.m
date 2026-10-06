function test_cell_land_cover_and_soil()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_cell_land_cover_and_soil
%%
%% Per-cell land cover and SSURGO soil on the Hamden parcel grid (10 m,
%% diagnostic coverage rule):
%%   1  the cell footprints tile the field
%%   2  raster cover: each cell's class areas equal an independent
%%      intersection of its footprint with the pixels, done in the raster's
%%      own CRS, and the cover fractions follow the class table
%%   3  manual cover, uniform and per cell, and the refusals
%%   4  the cover reaches each cell's parameter row (Ccrown, Cbare, Curb,
%%      Cwat, Crock) and its input identity, and only when given
%%   5  SSURGO: map-unit shares per cell, the dominant major component,
%%      mineral-horizon depth weighting, the CEC basis and restriction
%%      depth, and one named component for every cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_cell_land_cover_and_soil ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code,fullfile(v2,'interpolation'),fullfile(v2,'spinup'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
    function ok = raises(fn,id)
        ok = false;
        try; fn(); catch ME
            ok = strcmp(ME.identifier,id);
            if ~ok; fprintf('    raised %s: %s\n',ME.identifier,ME.message); end
        end
    end
names = {'cover_tree','cover_grass','cover_urban','cover_water','cover_bare','cover_rock'};
shp = fullfile(v2,'Spatial inputs','Hamden Blue Trail','Hamden_Blue_Trail_Parcel.shp');
setups = Build_Shapefile_Grid_Setups(shp,10,[],0.5);
topo = setups(1);
tmp = tempname; mkdir(tmp); cleanup = onCleanup(@() rmdir(tmp,'s'));

%% 1  footprints
fp = Cell_Footprints(topo);
check('the cell footprints tile the field', ...
    abs(sum(area(fp)) - topo.field_area_m2) <= 1e-6*topo.field_area_m2, ...
    sprintf('%.6g vs %.6g m2',sum(area(fp)),topo.field_area_m2));

%% 2  raster cover from a synthetic 30 m raster in the grid's UTM zone
bx = [min(topo.field_boundary_projected_x_m) max(topo.field_boundary_projected_x_m)];
by = [min(topo.field_boundary_projected_y_m) max(topo.field_boundary_projected_y_m)];
xl = [floor(bx(1)/30)*30-60 ceil(bx(2)/30)*30+60];
yl = [floor(by(1)/30)*30-60 ceil(by(2)/30)*30+60];
nr = round(diff(yl)/30); nc = round(diff(xl)/30);
[J,I] = meshgrid(1:nc,1:nr);
pal = [41 21 90 31];
A = uint8(pal(1+mod(I+2*J,4)));
R = maprefcells(xl,yl,[nr nc],'ColumnsStartFrom','north');
tif = fullfile(tmp,'cover.tif');
geotiffwrite(tif,A,R,'CoordRefSysCode',topo.projected_crs_epsg);
[cover,diag] = Cell_Land_Cover(topo,'raster',tif);
expected = zeros(topo.N_cells,numel(diag.class));
for c = 1:topo.N_cells
    [vx,vy] = Grid_Local_To_Projected(fp(c).Vertices(:,1),fp(c).Vertices(:,2),topo);
    cell_utm = polyshape(vx,vy);
    for i = 1:nr
        for j = 1:nc
            x0 = xl(1)+(j-1)*30; y1 = yl(2)-(i-1)*30;
            a = area(intersect(cell_utm,polyshape([x0 x0+30 x0+30 x0],[y1 y1 y1-30 y1-30])));
            if a > 0; k = diag.class == A(i,j); expected(c,k) = expected(c,k) + a; end
        end
    end
end
check('each cell''s class areas equal the intersection done in the raster''s CRS', ...
    max(abs(diag.class_area_m2(:)-expected(:))) <= 1e-6*max(expected(:)), ...
    sprintf('worst %.3g m2',max(abs(diag.class_area_m2(:)-expected(:)))));
s = diag.class_area_m2 ./ sum(diag.class_area_m2,2);
k = @(code) s(:,diag.class == code);
check('the cover fractions follow the NLCD class table', ...
    max(abs(cover.cover_tree - (k(41)+k(90)))) < 1e-12 && ...
    max(abs(cover.cover_grass - 0.9*k(21))) < 1e-12 && ...
    max(abs(cover.cover_urban - 0.1*k(21))) < 1e-12 && ...
    max(abs(cover.cover_bare - k(31))) < 1e-12);
check('every cell''s cover sums to one', max(abs(sum(cover{:,names},2)-1)) < 1e-12);
bad_classes = NLCD_Cover_Class_Table(); bad_classes.cover_tree(1) = 0.5;
check('a class table whose shares do not sum to one is refused', ...
    raises(@() Cell_Land_Cover(topo,'raster',tif,'ClassTable',bad_classes), ...
    'Cell_Land_Cover:BadClassTable'));
missing_class = NLCD_Cover_Class_Table(); missing_class(missing_class.class == 31,:) = [];
check('a raster class missing from the table is refused', ...
    raises(@() Cell_Land_Cover(topo,'raster',tif,'ClassTable',missing_class), ...
    'Cell_Land_Cover:UnmappedClass'));

%% 3  manual cover
U = Cell_Land_Cover(topo,'manual',struct('cover_tree',0.6,'cover_grass',0.3, ...
    'cover_urban',0.1,'cover_water',0,'cover_bare',0,'cover_rock',0));
check('one manual row applies to every cell', ...
    all(U.cover_tree == 0.6 & U.cover_grass == 0.3 & U.cover_urban == 0.1));
P = U(:,[{'cell_id'} names]); P.cover_tree(1) = 0.2; P.cover_grass(1) = 0.7;
P = P(end:-1:1,:);
U2 = Cell_Land_Cover(topo,'manual',P);
check('per-cell manual cover is matched by cell_id', ...
    U2.cover_tree(1) == 0.2 && U2.cover_grass(1) == 0.7 && all(U2.cover_tree(2:end) == 0.6));
P.cover_tree(2) = 0.9;
check('manual cover that does not sum to one is refused', ...
    raises(@() Cell_Land_Cover(topo,'manual',P),'Cell_Land_Cover:BadFractions'));

%% 4  into the cell's parameter row
n = topo.N_cells;
T = table(topo.cell_id(:),5.3*ones(n,1),nan(n,1),false(n,1),0.6*ones(n,1),0.1*ones(n,1), ...
    nan(n,1),nan(n,1),false(n,1),0.03*ones(n,1),nan(n,1),false(n,1),2*ones(n,1),nan(n,1),false(n,1), ...
    'VariableNames',{'cell_id','baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction','measured_sand_fraction', ...
    'measured_clay_fraction','use_measured_texture','baseline_OM_fraction', ...
    'measured_OM_fraction','use_measured_OM','baseline_CEC_cmolc_kg', ...
    'measured_CEC_cmolc_kg','use_measured_CEC'});
MOD = table(5.3,0.6,0.1,0.03,300,1295,1,"[0 0]",'VariableNames', ...
    {'PH','Psan','Pcla','Porg','Zbio','rsd','Cbare','Ccrown'});
ERW = table(2,'VariableNames',{'CEC_tot_in'});
plain = Configure_Cell_Agronomic_Inputs(MOD,ERW,Resolve_Cell_Agronomic_Inputs(topo,T));
check('without cover the row keeps its own cover and identity', ...
    plain(1).MOD_PARAM.Ccrown == "[0 0]" && ~ismember('Curb',plain(1).MOD_PARAM.Properties.VariableNames) && ...
    ~contains(plain(1).agronomic_identity,'cover='));
Tc = [T cover(:,names)];
withc = Configure_Cell_Agronomic_Inputs(MOD,ERW,Resolve_Cell_Agronomic_Inputs(topo,Tc));
ok = true;
for c = 1:n
    M = withc(c).MOD_PARAM;
    ok = ok && M.Ccrown == string(sprintf('[%.12g %.12g]',cover.cover_tree(c),cover.cover_grass(c))) && ...
        M.Cbare == cover.cover_bare(c) && M.Curb == cover.cover_urban(c) && ...
        M.Cwat == cover.cover_water(c) && M.Crock == cover.cover_rock(c);
end
check('each cell''s row carries its crowns, bare, urban, water and rock', ok);
check('the cover enters the cell''s input identity', contains(withc(1).agronomic_identity,'cover='));
check('a partial set of cover columns is refused', raises(@() Resolve_Cell_Agronomic_Inputs(topo, ...
    [T cover(:,{'cover_tree','cover_grass'})]),'Resolve_Cell_Agronomic_Inputs:MissingColumns'));

%% 5  SSURGO: two synthetic map units splitting the field west / east
lonb = topo.field_longitude_bounds; latb = topo.field_latitude_bounds;
mid = mean(lonb); m = 0.002;
wkt = @(a,b) sprintf('POLYGON ((%.9f %.9f, %.9f %.9f, %.9f %.9f, %.9f %.9f, %.9f %.9f))', ...
    a,latb(1)-m,b,latb(1)-m,b,latb(2)+m,a,latb(2)+m,a,latb(1)-m);
mapunits = table(["1";"2"],["W";"E"],string({wkt(lonb(1)-m,mid);wkt(mid,lonb(2)+m)}), ...
    'VariableNames',{'mukey','musym','wkt'});
hz = cell2table({ ...
    "1",11,"Westsoil",80,"Oe",0,5,NaN,NaN,90,0.4,NaN,NaN,NaN
    "1",11,"Westsoil",80,"A",5,20,60,10,4,1.2,5.0,NaN,2
    "1",11,"Westsoil",80,"B",20,60,40,20,1,1.5,5.5,NaN,3
    "1",12,"Minor",20,"A",0,30,90,2,1,1.4,6.0,4,NaN
    "2",21,"Eastsoil",60,"Ap",0,25,30,15,5,1.0,6.0,8,NaN
    "2",21,"Eastsoil",60,"Bw",25,50,35,12,2,1.4,6.2,6,NaN}, ...
    'VariableNames',{'mukey','cokey','compname','comppct_r','hzname','hzdept_r','hzdepb_r', ...
    'sandtotal_r','claytotal_r','om_r','dbthirdbar_r','ph1to1h2o_r','cec7_r','ecec_r'});
restr = table("2","Eastsoil",50,'VariableNames',{'mukey','compname','resdept_r'});
soil = Cell_SSURGO_Soil(topo,mapunits,hz,'Restrictions',restr,'Depth',[0 30]);
utm = projcrs(topo.projected_crs_epsg);
[ex,ey] = projfwd(utm,[latb(1)-m latb(2)+m],[mid mid]);
share_w = zeros(n,1);
for c = 1:n
    [vx,vy] = Grid_Local_To_Projected(fp(c).Vertices(:,1),fp(c).Vertices(:,2),topo);
    cell_utm = polyshape(vx,vy);
    west = polyshape([ex(1)-1e4 ex(1) ex(2) ex(2)-1e4],[ey(1) ey(1) ey(2) ey(2)]);
    share_w(c) = area(intersect(cell_utm,west))/area(cell_utm);
end
check('each cell''s map-unit shares match an independent split', ...
    max(abs(soil.share_W - share_w)) < 1e-6 && max(abs(soil.share_W + soil.share_E - 1)) < 1e-6, ...
    sprintf('worst %.3g',max(abs(soil.share_W - share_w))));
check('each cell takes the map unit covering most of it', ...
    all((soil.musym == "W") == (soil.share_W > 0.5)));
w = soil.musym == "W"; e = ~w;
check('a cell takes the dominant major component, depth-weighted over mineral horizons', ...
    all(soil.component(w) == "Westsoil") && all(abs(soil.sand_fraction(w) - 0.52) < 1e-12) && ...
    all(abs(soil.OM_fraction(w) - (4*15+1*10)/25/100) < 1e-12) && ...
    all(soil.component(e) == "Eastsoil") && all(abs(soil.sand_fraction(e) - (30*25+35*5)/30/100) < 1e-12));
check('CEC is CEC7 where given, ECEC otherwise', ...
    all(soil.CEC_basis(w) == "ECEC") && all(abs(soil.CEC_cmolc_kg(w) - 2.4) < 1e-12) && ...
    all(soil.CEC_basis(e) == "CEC7") && all(abs(soil.CEC_cmolc_kg(e) - (8*25+6*5)/30) < 1e-12));
check('the restriction depth is the component''s', ...
    all(soil.restriction_depth_cm(e) == 50) && all(isnan(soil.restriction_depth_cm(w))));
one = Cell_SSURGO_Soil(topo,mapunits,hz,'Component','Eastsoil');
check('a named component gives every cell its soil and keeps each cell''s map unit', ...
    all(one.component == "Eastsoil") && isequal(one.musym,soil.musym));

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_cell_land_cover_and_soil: %d check(s) failed',n_fail); end
end
