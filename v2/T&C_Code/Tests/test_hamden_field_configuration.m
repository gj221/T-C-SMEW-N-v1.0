function test_hamden_field_configuration()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_hamden_field_configuration
%%
%% The Hamden (145 Blue Trail, CT) field spin-up as the cluster configures
%% it, without spinning anything up:
%%   1  the site's MOD and ERW column files read as the workbook rows hold
%%      them (numbers, NaN, text); they leave the borrowed row its soil
%%      layers and set no rain chemistry
%%   2  a per-cell file without base saturation is refused: the site seeds
%%      the exchanger from it (f_CEC_inventory), and the Kantola row's 80 %
%%      would otherwise stand in
%%   0  the process grid: 5 m squares kept by the majority rule on the
%%      lattice of least misrepresented area, as the cluster builds it
%%      (the environment example's TC_CELL_M and TC_MIN_COVERAGE)
%%   3  the site's own per-cell file: every cell runs the site's settings
%%      over the borrowed case (location, its default 20-layer column with
%%      free drainage, vegetation file, CO2 and deposition records, the 2020
%%      handoff); its exchanger seeded from the representative base
%%      saturation (Representative_Soil_Chemistry) on the representative
%%      CEC; pH by a 1:1 extract; no rock; Cheshire texture, its own map
%%      unit's pH and OM, and its own land cover
%%   4  a column override may not set a cell's own soil or cover, nor an
%%      unknown column
%%   5  a cell with open grass assembles the site's parameters: the forcing
%%      reference height 2 m above the 24 m canopy, the default column, its
%%      own cover, and the vegetation, management and deposition of the
%%      Hamden vegetation file: event deposition (CT15 wet on the day's
%%      precipitation, TDep dry every day) and no fertiliser, nor any
%%      fertiliser year for the plot-scale frame to refill; a cold start
%%      takes the base file's 20-layer moisture profile
%%   6  with the parcel's DEM window (not versioned; skipped without it) the
%%      flux-tracking plan the cluster runs first: one class (every cell
%%      shares soil and vegetation), a topographically stratified ~10 %
%%      tracked, the rest to the forest, every cell once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_hamden_field_configuration ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code,fullfile(code,'Parameters'),fullfile(code,'Inputs'),fullfile(v2,'spinup'), ...
    fullfile(v2,'interpolation'));
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

grid_dir = fullfile(code,'Inputs','Field_Spinup','Hamden');
shp = fullfile(v2,'Spatial inputs','Hamden Blue Trail','Hamden_Blue_Trail_Parcel.shp');
agro = fullfile(grid_dir,'Cell_Agronomic_Hamden_Blue_Trail_Parcel_polygon_01.csv');
tmp = tempname; mkdir(tmp); cleanup = onCleanup(@() rmdir(tmp,'s'));
grid = {'GridCellSizeM',5,'MinCoverage','majority'};

%% --- 0 --------------------------------------------------------------------
env = fileread(fullfile(v2,'cluster','field_spinup_hamden.env.example'));
check('the cluster builds the same grid (TC_CELL_M=5, TC_MIN_COVERAGE=majority)', ...
    ~isempty(regexp(env,'(?m)^TC_CELL_M=5\s','once')) && ...
    ~isempty(regexp(env,'(?m)^TC_MIN_COVERAGE=majority\s','once')) && ...
    contains(env,'Inputs/Field_Spinup/Hamden/Cell_Agronomic_'));
setups = Build_Shapefile_Grid_Setups(shp,5,[],'majority'); topo = setups(1);
g = topo.grid_calibration;
check('the grid: every kept square more than half in the parcel, 2.3 % misrepresented', ...
    topo.N_cells == 936 && all(topo.cell_overlap_fraction > 0.5) && ...
    abs(g.misrepresented_area_m2-547.8) < 0.5 && abs(g.missed_field_area_m2-273.3) < 0.5 && ...
    abs(sum(topo.cell_overlap_area_m2)+topo.field_missed_area_m2-topo.field_area_m2) < 1e-6, ...
    sprintf('%d cells, misrepresented %.1f m2',topo.N_cells,g.misrepresented_area_m2));

%% --- 1 --------------------------------------------------------------------
M = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_MOD_columns.csv'));
E = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_ERW_columns.csv'));
check('the MOD columns read as the workbook holds them', ...
    M.Lat == 41.437 && M.Lon == -72.866 && isnan(M.Kbot) && M.Zbio == 300 && ...
    M.NumberOfDays == 2464 && M.Location == "Hamden" && ...
    ~isfield(M,'Zs') && ~isfield(M,'rain_chem_file') && ...
    M.meteo_start == "2015-01-01" && M.meteo_end == "2019-12-31" && ...
    M.transient_start == "2020-01-01" && M.transient_end == "2026-09-29");
check('the ERW columns read as the workbook holds them', ...
    E.fX_conc_convention == "f_CEC_inventory" && E.pH_DilutionFactor == 1 && ...
    E.M_rock_in == 0 && E.M_rock_event == 0 && E.Location == "Hamden");
%% The exchange selectivities are tabulated by USDA texture class
%% (ERW_MOD_PARAM_excel); Cheshire (55-62 % sand, 9.5 % clay) is a sandy loam.
erw_src = fileread(fullfile(code,'Parameters','ERW_MOD_PARAM_excel.m'));
known = regexp(erw_src,'ismember\(soil, \{([^}]*)\}\)','tokens');
known = regexp(strjoin(cellfun(@(t) t{1},known,'UniformOutput',false),','),'''([^'']+)''','tokens');
known = cellfun(@(t) t{1},known,'UniformOutput',false);
check('the soil class is one the ERW exchange constants know: a sandy loam', ...
    E.soil == "sandy loam" && ismember(char(E.soil),known), char(E.soil));
check('a column file without a value is refused', raises(@() ...
    Read_Parameter_Columns(write_text(fullfile(tmp,'blank.csv'),"a,b"+newline+"1,")), ...
    'Read_Parameter_Columns:BadFile'));

%% --- 2 --------------------------------------------------------------------
dry = [grid {'Method','per_cell','MODParamColumns',M,'ERWParamColumns',E, ...
    'OutputFile',fullfile(tmp,'CellBGStates.mat'),'DryRun',true}];
T = readtable(agro);
no_bs = removevars(T,{'baseline_bs_Ca_fraction','baseline_bs_Mg_fraction','baseline_bs_K_fraction'});
agro_no_bs = fullfile(tmp,'agronomic_no_bs.csv'); writetable(no_bs,agro_no_bs);
check('cells without a base saturation are refused before any work', ...
    raises(@() Run_Field_Cell_Spinup(shp,agro_no_bs,dry{:}),'Run_Field_Cell_Spinup:NoBaseSaturation'));

%% --- 3 --------------------------------------------------------------------
chem = readtable(fullfile(code,'Inputs','Hamden','Soils','Representative_soil_chemistry.csv'));
review = readtable(fullfile(grid_dir,'Cell_Review_SSURGO_and_NLCD.csv'),'TextType','string');
agro_bs = agro;
[~,CI] = Run_Field_Cell_Spinup(shp,agro_bs,dry{:});
N = numel(CI);
mod_ok = true; erw_ok = true; own_ok = true;
for c = 1:N
    m = CI(c).MOD_PARAM; e = CI(c).ERW_PARAM;
    r = T(T.cell_id == CI(c).cell_id,:);
    mod_ok = mod_ok && m.Location == "Hamden" && m.Lat == 41.437 && isnan(m.Kbot) && ...
        m.pft_override == "PFT_Hamden_Deciduous_Grass" && ...
        m.co2_file == "Ca_Data_NOAA_MLO_daily.mat" && m.forcing_file == "Data_Hamden.mat" && ...
        (~ismember('rain_chem_file',m.Properties.VariableNames) || ...
         strlength(strtrim(string(m.rain_chem_file))) == 0) && ...
        m.transient_start == "2020-01-01";
    f = str2num(char(e.f_CEC_in)); %#ok<ST2NM>
    fs = [chem.bs_Ca_fraction chem.bs_Mg_fraction chem.bs_K_fraction];
    erw_ok = erw_ok && e.fX_conc_convention == "f_CEC_inventory" && e.pH_DilutionFactor == 1 && ...
        e.M_rock_in == 0 && abs(e.CEC_tot_in - chem.CEC_cmolc_kg) < 1e-9 && ...
        max(abs(f - [fs 0 0 1-sum(fs)])) < 1e-9;
    q = review(review.cell_id == CI(c).cell_id,:);
    own_ok = own_ok && m.PH == q.ssurgo_pH_1to1_water && m.Porg == q.ssurgo_OM_fraction && ...
        m.Psan == r.baseline_sand_fraction && ...   % the row holds the cover to 12 digits
        max(abs(str2num(char(m.Ccrown)) - [r.cover_tree r.cover_grass])) < 1e-11; %#ok<ST2NM>
end
check(sprintf('all %d cells run the site''s MOD settings',N),N == topo.N_cells && mod_ok);
check('every cell''s exchanger is the representative one, pH by a 1:1 extract, no rock',erw_ok);
check('every cell keeps its own map unit''s pH and OM and its own land cover',own_ok);
pHs = arrayfun(@(ci) ci.MOD_PARAM.PH,CI);
check('the cells carry both map units'' pH (Cheshire and Wilbraham)', ...
    isequal(unique(pHs(:)).',[5.3 5.5]) && chem.n_layers >= 1 && ...
    sum(fs) < 1 && abs(chem.target_pH - mean(pHs)) < 0.05);
check('the baseline base saturation is part of each cell''s identity', ...
    contains(CI(1).agronomic_identity,':baseline'));

%% --- 4 --------------------------------------------------------------------
M_pH = M; M_pH.PH = 6;
check('a MOD override of a cell''s own soil is refused', raises(@() ...
    Run_Field_Cell_Spinup(shp,agro_bs,dry{:},'MODParamColumns',M_pH),'Run_Field_Cell_Spinup:CellColumn'));
M_cov = M; M_cov.Ccrown = "[1 0]";
check('a MOD override of a cell''s own cover is refused', raises(@() ...
    Run_Field_Cell_Spinup(shp,agro_bs,dry{:},'MODParamColumns',M_cov),'Run_Field_Cell_Spinup:CellColumn'));
E_cec = E; E_cec.CEC_tot_in = 10;
check('an ERW override of a cell''s own CEC is refused', raises(@() ...
    Run_Field_Cell_Spinup(shp,agro_bs,dry{:},'ERWParamColumns',E_cec),'Run_Field_Cell_Spinup:CellColumn'));
E_bad = E; E_bad.not_a_column = 1;
check('an unknown ERW column is refused', raises(@() ...
    Run_Field_Cell_Spinup(shp,agro_bs,dry{:},'ERWParamColumns',E_bad), ...
    'Run_Field_Cell_Spinup:UnknownERWParamColumn'));

%% --- 5 --------------------------------------------------------------------
cover_grass = arrayfun(@(ci) str2num(char(ci.MOD_PARAM.Ccrown)),CI,'UniformOutput',false); %#ok<ST2NM>
cover_grass = cellfun(@(x) x(2),cover_grass);
[~,c] = max(cover_grass);
Zs = eval(char(CI(c).MOD_PARAM.Zs));
dn = (datenum(2015,1,1):1/24:datenum(2015,12,31,23,0,0))';
[YE,MO,DA,HO] = datevec(dn);
K = load(fullfile(code,'Inputs','Data_Kantola_BG.mat'));
ctx = struct('MOD_PARAM_input',CI(c).MOD_PARAM,'ERW_PARAM_input',CI(c).ERW_PARAM, ...
    'NN',numel(dn),'NNd',numel(dn)/24,'ms',numel(Zs)-1,'Datam',[YE MO DA HO], ...
    'Lat',41.437,'Lon',-72.866,'DeltaGMT',-5,'t_bef',-0.67,'t_aft',1.67, ...
    'Ta',10*ones(numel(dn),1),'Ca',415*ones(numel(dn),1), ...
    'OPT_SoilBiogeochemistry',1,'NN_case_study_ERW_MOD',5, ...
    'B_out',K.B_out(end,:),'O_out',0.3*ones(1,numel(Zs)-1),'NavlI_out',K.NavlI_out(end,:), ...
    'RexmyI_out',K.RexmyI_out(end,:),'NupI_L_out',K.NupI_L_out);
P = Load_Cell_Parameters(fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(CI(c).MOD_PARAM.pft_override) '.m']),ctx);
check('the forcing reference height is 2 m above the 24 m canopy', ...
    P.zatm == 26 && max(P.hc_H(1,:)) == 24, sprintf('zatm %g',P.zatm));
check('the soil column is the default: 20 layers to 5 m, with free drainage', ...
    P.Zs(end) == 5000 && P.ms == 20 && isnan(P.Kbot));
check('the cell''s own cover reaches the parameters', ...
    cover_grass(c) > 0 && max(abs(P.Ccrown - str2num(char(CI(c).MOD_PARAM.Ccrown)))) < 1e-12 && ... %#ok<ST2NM>
    abs(sum(P.Ccrown) + P.Cbare + P.Curb + P.Cwat + P.Crock - 1) < 1e-9);
R = load(fullfile(code,'Inputs','Deposition_Hamden.mat'),'Deposition');
ev = R.Deposition.events; row = find(ev.year == 2021);
[~,rain_day] = Deposition_Day(P.B_IO,[2021 7 1 1 0 0],10);
check('the Hamden vegetation, mowing and event deposition are assembled', ...
    P.ZR95_H(1) == 1000 && P.ZR95_L(1) == 250 && ...
    any(abs(P.Mpar_L(1).jDay_cut - datenum(2015,9,1)) < 1e-9) && P.Mpar_L(1).fract_left == 0.1 && ...
    isequal(P.B_IO.Dep_events,ev) && P.B_IO.DepN == R.Deposition.mean_g_m2_d(1) && ...
    P.B_IO.DepCa == R.Deposition.mean_g_m2_d(3) && ...
    max(abs(rain_day([1 3 4 5])-(ev.dry_g_m2_d(row,:)+10*ev.wet_g_m2_per_mm(row,:)))) < 1e-15);
check('no fertiliser, and no fertiliser year for the plot-scale frame to refill', ...
    sum(P.B_IO.FertN) == 0 && isempty(P.B_IO.FertN_years) && ...
    ~any(ismember(2014:2026,P.B_IO.FertN_years)));
ctx_cold = ctx; ctx_cold.O_out = zeros(1,numel(Zs)-1);
P_cold = Load_Cell_Parameters(fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(CI(c).MOD_PARAM.pft_override) '.m']),ctx_cold);
check('a cold start takes the base file''s 20-layer moisture profile', ...
    isequal(size(P_cold.O(1,:)),[1 20]) && all(P_cold.O(1,:) > 0));

%% --- 6 --------------------------------------------------------------------
dem = fullfile(v2,'Spatial inputs','Hamden Blue Trail','DEM');
if isfolder(dem) && ~isempty(dir(fullfile(dem,'*.tif')))
    plan_dir = fullfile(tmp,'plan');
    plan = readtable(Run_Field_Cell_Spinup(shp,agro,grid{:},'MODParamColumns',M, ...
        'ERWParamColumns',E,'FieldDEM',dem,'OutputFile',fullfile(plan_dir,'CellBGStates.mat'), ...
        'Stage','plan'),'TextType','string');
    tracked = logical(plan.is_tracked);
    check('the flux-tracking plan: one class, ~10 % tracked, the rest to the forest', ...
        height(plan) == N && all(plan.class_index == 1) && nnz(plan.is_class_seed) == 1 && ...
        nnz(tracked) >= 0.08*N && nnz(tracked) <= 0.15*N && ...
        all(plan.method(tracked) == "replay") && all(plan.method(~tracked) == "rf"), ...
        sprintf('%d tracked of %d',nnz(tracked),N));
else
    fprintf('  SKIP  flux-tracking plan: no DEM in %s\n',dem);
end

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_hamden_field_configuration: %d check(s) failed',n_fail); end
end

function f = write_text(f,txt)
fid = fopen(f,'w'); fprintf(fid,'%s\n',txt); fclose(fid);
end
