function files = Run_Hamden_Field_Inputs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The per-cell inputs of the Hamden (145 Blue Trail, CT) field spin-up on
%% its calibrated grid (95 % coverage rule), written to
%% T&C_Code/Inputs/Field_Spinup/Hamden_calibrated_grid:
%%   soil       Cheshire in every cell (SSURGO, 0-30 cm, depth-weighted
%%              mineral horizons): the sand, clay and organic-matter
%%              baselines;
%%   chemistry  the site's representative bulk pH, CEC and base saturation
%%              (Inputs/Hamden/Soils/Representative_soil_chemistry.csv) as
%%              every cell's baselines, until sub-plot measurements fill the
%%              measured_* columns (Krige_Cell_Agronomic_Values,
%%              Apply_Kriged_Cell_Values) and a warm start re-spins the
%%              field from this one (WarmStartStates);
%%   cover      NLCD 2020 shares per cell by exact area overlap
%%              (Cell_Land_Cover).
%% Also writes a review table of every cell's dominant SSURGO map unit and
%% component and its NLCD class areas.
%%
%% Name/value
%%   'Chemistry'  another one-row representative set (same columns)
%%   'Output'     another output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
code = fullfile(v2,'T&C_Code');
addpath(code,here,fullfile(v2,'spinup'));
site = fullfile(v2,'Spatial inputs','Hamden Blue Trail');
soil_dir = fullfile(code,'Inputs','Hamden','Soils');
opt = struct('Chemistry',fullfile(soil_dir,'Representative_soil_chemistry.csv'), ...
    'Output',fullfile(code,'Inputs','Field_Spinup','Hamden_calibrated_grid'));
for k = 1:2:numel(varargin)
    if ~isfield(opt,char(varargin{k}))
        error('Run_Hamden_Field_Inputs:BadOptions','Unknown option %s.',char(varargin{k}));
    end
    opt.(char(varargin{k})) = varargin{k+1};
end

%% the representative set first: nothing is written without it
chem = readtable(opt.Chemistry,'TextType','string','VariableNamingRule','preserve');
need = {'bulk_pH_1to1_water','CEC_cmolc_kg','bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction'};
if height(chem) ~= 1 || ~all(ismember(need,chem.Properties.VariableNames))
    error('Run_Hamden_Field_Inputs:BadChemistry', ...
        '%s must hold one row with %s.',opt.Chemistry,strjoin(need,', '));
end
v = cellfun(@(n) double(chem.(n)),need);
if any(~isfinite(v))
    error('Run_Hamden_Field_Inputs:ChemistryNotSet', ...
        ['%s: %s not set. The representative bulk pH, CEC and base saturation ' ...
         'are the site''s own choice (Inputs/Hamden/README.md, Soils).'], ...
        opt.Chemistry,strjoin(need(~isfinite(v)),', '));
end

shp = fullfile(site,'Hamden_Blue_Trail_Parcel.shp');
tif = fullfile(site,'Land cover','NLCD_2020_Hamden_Blue_Trail_window.tif');
files = Write_Field_Agronomic_Templates(shp,opt.Output,[],[],[]);
T = readtable(files{1});
setups = Build_Shapefile_Grid_Setups(shp,[],[],[]); topo = setups(1);
if ~isequal(T.cell_id(:),topo.cell_id(:))
    error('Run_Hamden_Field_Inputs:GridMismatch','The template and the grid disagree.');
end
mapunits = readtable(fullfile(soil_dir,'SSURGO_mapunit_polygons_parcel_buffer.csv'), ...
    'TextType','string','Delimiter',',');
hz = readtable(fullfile(soil_dir,'SSURGO_parcel_horizons.csv'),'TextType','string');
rs = readtable(fullfile(soil_dir,'SSURGO_parcel_restrictions.csv'),'TextType','string');
cheshire = Cell_SSURGO_Soil(topo,mapunits,hz,'Restrictions',rs,'Component','Cheshire');
dominant = Cell_SSURGO_Soil(topo,mapunits,hz,'Restrictions',rs);
[cover,diag] = Cell_Land_Cover(topo,'raster',tif);

N = height(T);
T.baseline_sand_fraction = cheshire.sand_fraction;
T.baseline_clay_fraction = cheshire.clay_fraction;
T.baseline_OM_fraction = cheshire.OM_fraction;
T.baseline_bulk_pH = repmat(v(1),N,1);
T.baseline_CEC_cmolc_kg = repmat(v(2),N,1);
T.baseline_bs_Ca_fraction = repmat(v(3),N,1);
T.baseline_bs_Mg_fraction = repmat(v(4),N,1);
T.baseline_bs_K_fraction = repmat(v(5),N,1);
names = {'cover_tree','cover_grass','cover_urban','cover_water','cover_bare','cover_rock'};
T = [T cover(:,names)];
Resolve_Cell_Agronomic_Inputs(topo,T);   % the file must resolve before it is written
writetable(T,files{1});

R = dominant(:,{'cell_id','mukey','musym','mapunit_share','component','sand_fraction', ...
    'clay_fraction','OM_fraction','pH_1to1_water','CEC_cmolc_kg','CEC_basis','restriction_depth_cm'});
R = renamevars(R,R.Properties.VariableNames(2:end),strcat('ssurgo_',R.Properties.VariableNames(2:end)));
R.footprint_area_m2 = diag.footprint_area_m2;
for k = 1:numel(diag.class)
    R.(sprintf('nlcd_%d_area_m2',diag.class(k))) = diag.class_area_m2(:,k);
end
R = [R cover(:,names)];
review = fullfile(opt.Output,'Cell_Review_SSURGO_and_NLCD.csv');
writetable(R,review);
files{end+1} = review;
a = diag.footprint_area_m2;
fprintf(['Run_Hamden_Field_Inputs: %d cells of %.4g m -> %s\n' ...
    '  Cheshire sand %.3f clay %.3f OM %.4f; representative pH %.2f, CEC %.3g cmolc/kg, ' ...
    'base saturation Ca %.3f Mg %.3f K %.3f\n' ...
    '  cover (area-weighted): tree %.4f grass %.4f urban %.4f\n'], ...
    N,topo.cell_size_m,files{1},T.baseline_sand_fraction(1),T.baseline_clay_fraction(1), ...
    T.baseline_OM_fraction(1),v(1),v(2),v(3),v(4),v(5), ...
    sum(cover.cover_tree.*a)/sum(a),sum(cover.cover_grass.*a)/sum(a),sum(cover.cover_urban.*a)/sum(a));
end
