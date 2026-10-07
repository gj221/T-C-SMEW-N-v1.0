function files = Run_Hamden_Field_Inputs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The per-cell inputs of the Hamden (145 Blue Trail, CT) field spin-up on
%% its process grid -- 5 m squares, the coarsest width at which the
%% parcel's terrain statistics have converged (Run_Topographic_Grid_Convergence),
%% kept by the majority rule on the lattice that minimises the
%% misrepresented area (Build_Projected_Field_Topology) -- written to
%% T&C_Code/Inputs/Field_Spinup/Hamden:
%%   texture    Cheshire in every cell (SSURGO, 0-30 cm, depth-weighted
%%              mineral horizons): the sand and clay baselines;
%%   pH, OM     each cell's own SSURGO map unit, its dominant component
%%              over the same depth: Cheshire on the Cheshire-Holyoke
%%              units, Wilbraham on the wet Wilbraham-Menlo unit;
%%   exchanger  the site's representative effective CEC and base
%%              saturation (Representative_Soil_Chemistry): pooled from the
%%              laboratory layers of the Cheshire pedons near the site
%%              (Inputs/Hamden/Soils/KSSL_Cheshire_pedons_0-30cm.csv) whose
%%              pH lies within 0.3 of the cells' area-mean pH, with
%%              SoilGrids' CEC at pH 7 as the cross-check; written to
%%              Inputs/Hamden/Soils/Representative_soil_chemistry.csv;
%%   cover      NLCD 2020 shares per cell by exact area overlap
%%              (Cell_Land_Cover).
%% All of these are baselines: sub-plot measurements fill the measured_*
%% columns (Krige_Cell_Agronomic_Values, Apply_Kriged_Cell_Values) and a
%% warm start re-spins the field from this one (WarmStartStates).  Also
%% writes a review table of every cell's dominant SSURGO map unit and
%% component and its NLCD class areas.
%%
%% Name/value
%%   'Chemistry'     a one-row representative set to use instead of the
%%                   derived one (columns CEC_cmolc_kg, bs_Ca_fraction,
%%                   bs_Mg_fraction, bs_K_fraction)
%%   'CellSizeM', 'CoverageRule'
%%                   another grid (default 5 and 'majority'; the spin-up's
%%                   TC_CELL_M and TC_MIN_COVERAGE must match)
%%   'Output'        another output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
code = fullfile(v2,'T&C_Code');
addpath(code,here,fullfile(v2,'spinup'));
site = fullfile(v2,'Spatial inputs','Hamden Blue Trail');
soil_dir = fullfile(code,'Inputs','Hamden','Soils');
opt = struct('Chemistry','','CellSizeM',5,'CoverageRule','majority', ...
    'Output',fullfile(code,'Inputs','Field_Spinup','Hamden'));
for k = 1:2:numel(varargin)
    if ~isfield(opt,char(varargin{k}))
        error('Run_Hamden_Field_Inputs:BadOptions','Unknown option %s.',char(varargin{k}));
    end
    opt.(char(varargin{k})) = varargin{k+1};
end

shp = fullfile(site,'Hamden_Blue_Trail_Parcel.shp');
tif = fullfile(site,'Land cover','NLCD_2020_Hamden_Blue_Trail_window.tif');
setups = Build_Shapefile_Grid_Setups(shp,opt.CellSizeM,[],opt.CoverageRule);
topo = setups(1);
mapunits = readtable(fullfile(soil_dir,'SSURGO_mapunit_polygons_parcel_buffer.csv'), ...
    'TextType','string','Delimiter',',');
hz = readtable(fullfile(soil_dir,'SSURGO_parcel_horizons.csv'),'TextType','string');
rs = readtable(fullfile(soil_dir,'SSURGO_parcel_restrictions.csv'),'TextType','string');
cheshire = Cell_SSURGO_Soil(topo,mapunits,hz,'Restrictions',rs,'Component','Cheshire');
dominant = Cell_SSURGO_Soil(topo,mapunits,hz,'Restrictions',rs);
[cover,diag] = Cell_Land_Cover(topo,'raster',tif);
a = diag.footprint_area_m2;
pH = dominant.pH_1to1_water; OM = dominant.OM_fraction;
if any(~isfinite(pH)) || any(~isfinite(OM))
    error('Run_Hamden_Field_Inputs:NoCellSoil', ...
        '%d cell(s) have no SSURGO pH or OM in 0-30 cm.',nnz(~isfinite(pH) | ~isfinite(OM)));
end

%% the representative exchanger, at the cells' own pH
target_pH = sum(a.*pH)/sum(a);
if isempty(opt.Chemistry)
    chem = Representative_Soil_Chemistry(fullfile(soil_dir,'KSSL_Cheshire_pedons_0-30cm.csv'), ...
        target_pH,'SoilGrids',fullfile(soil_dir,'SoilGrids_parcel_centre.csv'), ...
        'Output',fullfile(soil_dir,'Representative_soil_chemistry.csv'));
else
    chem = readtable(opt.Chemistry,'TextType','string','VariableNamingRule','preserve');
end
need = {'CEC_cmolc_kg','bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction'};
if height(chem) ~= 1 || ~all(ismember(need,chem.Properties.VariableNames)) || ...
        any(~isfinite(cellfun(@(n) double(chem.(n)),need)))
    error('Run_Hamden_Field_Inputs:BadChemistry', ...
        'The representative set must be one row with finite %s.',strjoin(need,', '));
end
v = cellfun(@(n) double(chem.(n)),need);

%% the blank template (it never overwrites, so it is written aside), filled
%% here and written over this builder's previous output
blank = tempname; cleanup = onCleanup(@() rmdir(blank,'s'));
template = Write_Field_Agronomic_Templates(shp,blank,opt.CellSizeM,[], ...
    opt.CoverageRule);
T = readtable(template{1});
if ~isequal(T.cell_id(:),topo.cell_id(:))
    error('Run_Hamden_Field_Inputs:GridMismatch','The template and the grid disagree.');
end
if ~isfolder(opt.Output); mkdir(opt.Output); end
[~,name,ext] = fileparts(template{1});
files = {fullfile(opt.Output,[name ext])};
N = height(T);
T.baseline_sand_fraction = cheshire.sand_fraction;
T.baseline_clay_fraction = cheshire.clay_fraction;
T.baseline_bulk_pH = pH;
T.baseline_OM_fraction = OM;
T.baseline_CEC_cmolc_kg = repmat(v(1),N,1);
T.baseline_bs_Ca_fraction = repmat(v(2),N,1);
T.baseline_bs_Mg_fraction = repmat(v(3),N,1);
T.baseline_bs_K_fraction = repmat(v(4),N,1);
names = {'cover_tree','cover_grass','cover_urban','cover_water','cover_bare','cover_rock'};
T = [T cover(:,names)];
Resolve_Cell_Agronomic_Inputs(topo,T);   % the file must resolve before it is written
writetable(T,files{1});

R = dominant(:,{'cell_id','mukey','musym','mapunit_share','component','sand_fraction', ...
    'clay_fraction','OM_fraction','pH_1to1_water','CEC_cmolc_kg','CEC_basis','restriction_depth_cm'});
R = renamevars(R,R.Properties.VariableNames(2:end),strcat('ssurgo_',R.Properties.VariableNames(2:end)));
R.footprint_area_m2 = a;
for k = 1:numel(diag.class)
    R.(sprintf('nlcd_%d_area_m2',diag.class(k))) = diag.class_area_m2(:,k);
end
R = [R cover(:,names)];
review = fullfile(opt.Output,'Cell_Review_SSURGO_and_NLCD.csv');
writetable(R,review);
files{end+1} = review;
[u,~,g] = unique(dominant.component);
parts = arrayfun(@(k) sprintf('%s %d cells (pH %.2f, OM %.4f)',u(k),nnz(g == k), ...
    pH(find(g == k,1)),OM(find(g == k,1))),1:numel(u),'UniformOutput',false);
fprintf(['Run_Hamden_Field_Inputs: %d cells of %.4g m -> %s\n' ...
    '  texture Cheshire (sand %.3f, clay %.3f); pH and OM by map unit: %s\n' ...
    '  exchanger at the cells'' mean pH %.2f: CEC %.3g cmolc/kg, base saturation Ca %.3f Mg %.3f K %.3f\n' ...
    '  cover (area-weighted): tree %.4f grass %.4f urban %.4f\n'], ...
    N,topo.cell_size_m,files{1},T.baseline_sand_fraction(1),T.baseline_clay_fraction(1), ...
    strjoin(parts,'; '),target_pH,v(1),v(2),v(3),v(4), ...
    sum(cover.cover_tree.*a)/sum(a),sum(cover.cover_grass.*a)/sum(a),sum(cover.cover_urban.*a)/sum(a));
end
