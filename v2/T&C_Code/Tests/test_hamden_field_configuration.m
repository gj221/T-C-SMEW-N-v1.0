function test_hamden_field_configuration()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_hamden_field_configuration
%%
%% The Hamden (145 Blue Trail, CT) field spin-up as the cluster configures
%% it, without spinning anything up:
%%   1  the site's MOD and ERW column files read as the workbook rows hold
%%      them (numbers, NaN, text); they leave the borrowed row its soil
%%      layers and name CT15's rain record
%%   2  the exchanger emerges from the spin-up (spinup_equilibrium): a
%%      per-cell file without base saturation runs, the representative one
%%      is not imposed, and a cell given a measured one is seeded from it
%%   0  the process grid: 5 m squares kept by the majority rule on the
%%      lattice of least misrepresented area, as the cluster builds it
%%      (the environment example's TC_CELL_M and TC_MIN_COVERAGE)
%%   3  the site's own per-cell file: every cell runs the site's settings
%%      over the borrowed case (location, its default 20-layer column with
%%      free drainage, vegetation file, CO2, deposition and rain records, the
%%      2020 handoff); its exchanger left to emerge on the representative
%%      CEC at pH 7 (Representative_Soil_Chemistry); pH by a 1:1 extract; no rock;
%%      Cheshire texture, its own map
%%      unit's pH and OM, and its own land cover
%%   4  a column override may not set a cell's own soil or cover, nor an
%%      unknown column
%%   5  a cell with trees and open grass assembles the site's parameters: the
%%      forcing reference height 2 m above the 24 m canopy, the default column, its
%%      own cover, and the vegetation, management and deposition of the
%%      Hamden vegetation file: event deposition (CT15 wet on the day's
%%      precipitation, TDep dry every day) and no fertiliser, nor any
%%      fertiliser year for the plot-scale frame to refill; a cold start
%%      takes the base file's 20-layer moisture profile; trees and grass
%%      are perennial, a run that follows another starting from its stand,
%%      and grow without T&C's sink limitation; the site's weathering keeps
%%      the uplift's Ca + Mg charge, split to the measured exchanger, a cold
%%      start's primary minerals at the steady state of their supply; the
%%      site's own exchange constants hold the measured exchanger on the
%%      class seed's soil solution; a crown without cover carries no
%%      vegetation (a wooded cell no open grass, a lawn no trees), and a
%%      crown with cover that the stand handed in leaves empty -- a lawn
%%      seeded from the wooded class seed -- starts from the file's stand
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
    ~isfield(M,'Zs') && M.rain_chem_file == "Rain_Chem_Hamden.mat" && ...
    M.meteo_start == "2015-01-01" && M.meteo_end == "2019-12-31" && ...
    M.transient_start == "2020-01-01" && M.transient_end == "2026-09-29");
check('the ERW columns read as the workbook holds them', ...
    E.fX_conc_convention == "spinup_equilibrium" && E.pH_DilutionFactor == 1 && ...
    E.M_rock_in == 0 && E.M_rock_event == 0 && E.Location == "Hamden");
check('Na and the strong anions come from deposition alone, not the steady closure', ...
    E.keyword_add == 0);
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
[~,CI0] = Run_Field_Cell_Spinup(shp,agro_no_bs,dry{:});
check('cells without a base saturation run: their exchanger emerges from the spin-up', ...
    all(arrayfun(@(ci) ci.ERW_PARAM.fX_conc_convention == "spinup_equilibrium",CI0)));
one = T;
one.measured_bs_Ca_fraction(1) = 0.40; one.measured_bs_Mg_fraction(1) = 0.10;
one.measured_bs_K_fraction(1) = 0.02; one.use_measured_base_saturation(1) = 1;
agro_one = fullfile(tmp,'agronomic_one_measured.csv'); writetable(one,agro_one);
[~,CI1] = Run_Field_Cell_Spinup(shp,agro_one,dry{:});
k = find([CI1.cell_id] == one.cell_id(1));
rest = CI1(setdiff(1:numel(CI1),k));
check('a cell given a measured base saturation is seeded from it; the rest emerge', ...
    isscalar(k) && CI1(k).ERW_PARAM.fX_conc_convention == "f_CEC_inventory" && ...
    max(abs(str2num(char(CI1(k).ERW_PARAM.f_CEC_in)) - [0.40 0.10 0.02 0 0 0.48])) < 1e-9 && ... %#ok<ST2NM>
    all(arrayfun(@(ci) ci.ERW_PARAM.fX_conc_convention == "spinup_equilibrium",rest)));

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
        m.rain_chem_file == "Rain_Chem_Hamden.mat" && ...
        m.transient_start == "2020-01-01";
    erw_ok = erw_ok && e.fX_conc_convention == "spinup_equilibrium" && e.pH_DilutionFactor == 1 && ...
        e.M_rock_in == 0 && abs(e.CEC_tot_in - chem.CEC_cmolc_kg) < 1e-9 && ...
        all(isnan(str2double(string(e.f_CEC_in))));
    q = review(review.cell_id == CI(c).cell_id,:);
    own_ok = own_ok && m.PH == q.ssurgo_pH_1to1_water && m.Porg == q.ssurgo_OM_fraction && ...
        m.Psan == r.baseline_sand_fraction && ...   % the row holds the cover to 12 digits
        max(abs(str2num(char(m.Ccrown)) - [r.cover_tree r.cover_grass])) < 1e-11; %#ok<ST2NM>
end
check(sprintf('all %d cells run the site''s MOD settings',N),N == topo.N_cells && mod_ok);
check('every cell''s exchanger emerges on the representative CEC, pH by a 1:1 extract, no rock',erw_ok);
check('every cell keeps its own map unit''s pH and OM and its own land cover',own_ok);
pHs = arrayfun(@(ci) ci.MOD_PARAM.PH,CI);
check('the cells carry both map units'' pH (Cheshire and Wilbraham)', ...
    isequal(unique(pHs(:)).',[5.3 5.5]) && chem.n_layers >= 1 && ...
    chem.bs_Ca_fraction + chem.bs_Mg_fraction + chem.bs_K_fraction < 1 && ...
    abs(chem.target_pH - mean(pHs)) < 0.05);
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
covers = arrayfun(@(ci) str2num(char(ci.MOD_PARAM.Ccrown)),CI,'UniformOutput',false); %#ok<ST2NM>
cover_tree = cellfun(@(x) x(1),covers);
cover_grass = cellfun(@(x) x(2),covers);
%% the grassiest cell that also has trees
mixed = find(cover_tree > 0 & cover_grass > 0);
[~,k] = max(cover_grass(mixed));
c = mixed(k);
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
[~,rain_day,rain_erw] = Deposition_Day(P.B_IO,[2021 7 1 1 0 0],10);
check('the Hamden vegetation, mowing and event deposition are assembled', ...
    P.ZR95_H(1) == 1000 && P.ZR95_L(1) == 250 && ...
    any(abs(P.Mpar_L(1).jDay_cut - datenum(2015,9,1)) < 1e-9) && P.Mpar_L(1).fract_left == 0.1 && ...
    isequal(P.B_IO.Dep_events,ev) && P.B_IO.DepN == R.Deposition.mean_g_m2_d(1) && ...
    P.B_IO.DepCa == R.Deposition.mean_g_m2_d(3) && ...
    max(abs(rain_day([1 3 4 5])-(ev.dry_g_m2_d(row,:)+10*ev.wet_g_m2_per_mm(row,:)))) < 1e-15);
check('the record also brings the ERW its Na and strong anions, wet and dry', ...
    max(abs(rain_erw-(ev.erw_dry_umol_m2_d(row,:)+10*ev.erw_wet_umol_m2_per_mm(row,:)))) < 1e-12 && ...
    all(rain_erw > 0));
check('no fertiliser, manure or N fixation, and no fertiliser year for the plot-scale frame to refill', ...
    sum(P.B_IO.FertN) == 0 && sum(P.B_IO.FertP) == 0 && sum(P.B_IO.FertK) == 0 && ...
    sum(P.B_IO.FertCa) == 0 && sum(P.B_IO.FertMg) == 0 && sum(P.B_IO.ManF) == 0 && ...
    isempty(P.B_IO.FertN_years) && ~any(ismember(2014:2026,P.B_IO.FertN_years)) && ...
    ~any([P.ParEx_H.bfix P.ParEx_L.bfix]));
ctx_cold = ctx; ctx_cold.O_out = zeros(1,numel(Zs)-1);
P_cold = Load_Cell_Parameters(fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(CI(c).MOD_PARAM.pft_override) '.m']),ctx_cold);
check('a cold start takes the base file''s 20-layer moisture profile', ...
    isequal(size(P_cold.O(1,:)),[1 20]) && all(P_cold.O(1,:) > 0));
%% Hamden's weathering: the uplift's Ca + Mg charge, split so that with the
%% deposition it supplies the dissolved Ca/Mg in which the measured exchanger
%% is at equilibrium under SCEPTER's selectivity; a cold start's
%% primary-mineral pools at the steady state of that supply, a run that
%% follows another keeping its own.
chem_r = (chem.bs_Ca_fraction/chem.bs_Mg_fraction)/10^(10.786-10.47);
mol = @(ca,mg) [ca/40.078 mg/24.305];
w_mol = mol(P.B_IO.Tup_Ca,P.B_IO.Tup_Mg);
d_mol = mol(R.Deposition.mean_g_m2_d(3),R.Deposition.mean_g_m2_d(4));
BGP = Biogeochemistry_Parameter(1);
ctx_seed = rmfield(ctx,{'B_out','O_out','NavlI_out','RexmyI_out','NupI_L_out'});
ctx_seed.O = zeros(1,numel(Zs)-1); ctx_seed.P = zeros(1,79);   % as Phase 2 starts a cold spin-up
P_seed = Load_Cell_Parameters(fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(CI(c).MOD_PARAM.pft_override) '.m']),ctx_seed);
check('weathering keeps the uplift''s Ca + Mg charge, split to the measured exchanger''s dissolved Ca/Mg', ...
    abs(sum(w_mol)/sum(mol(0.003,0.002)) - 1) < 1e-5 && ...
    abs((w_mol(1)+d_mol(1))/(w_mol(2)+d_mol(2))/chem_r - 1) < 1e-4 && ...
    P.B_IO.Tup_Si == 0.0001, ...
    sprintf('Tup_Ca %.7f Tup_Mg %.7f, supply Ca/Mg %.4f against %.4f',P.B_IO.Tup_Ca,P.B_IO.Tup_Mg, ...
    (w_mol(1)+d_mol(1))/(w_mol(2)+d_mol(2)),chem_r));
check('a cold start''s primary K, Ca and Mg begin at the steady state of their supply; a run that follows another keeps its own', ...
    max(abs(P_seed.P(1,[55 63 71]) - [P_seed.B_IO.Tup_K/BGP.K3k P_seed.B_IO.Tup_Ca/BGP.K3ca ...
    P_seed.B_IO.Tup_Mg/BGP.K3mg])) < 1e-9 && isequal(P.P(1,:),K.B_out(end,:)), ...
    sprintf('primary %s',mat2str(P_seed.P(1,[55 63 71]),6)));
%% Trees and grass are perennial: a cold start begins from the file's stand,
%% a run that follows another from that run's end state, as the 2-D driver
%% loads it for each cell (Cell_BG_State_Context -> Perennial_Vegetation_Apply);
%% the stand grows without T&C's sink limitation, as the legacy forest runs do.
[veg_names,veg_depth] = Perennial_Vegetation_States();
depth = struct('B',8,'NupI',6,'PARI',3,'NuLit',6,'ISOIL',27,'Rexmy',3);
Veg = struct();
for k = 1:numel(veg_names)
    if veg_depth(k)
        Veg.(veg_names{k}) = rand(1,2,depth.(regexprep(veg_names{k},'_[HL]$','')));
    else
        Veg.(veg_names{k}) = rand(1,2);
    end
end
Veg.hc_H = [24 0]; Veg.hc_L = [0.3 0.3];   % a carried canopy keeps the reference height
ctx_veg = ctx; ctx_veg.Veg_out = Veg;
P_veg = Load_Cell_Parameters(fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(CI(c).MOD_PARAM.pft_override) '.m']),ctx_veg);
carried = true;
for k = 1:numel(veg_names)
    carried = carried && isequal(P_veg.(veg_names{k})(1,:,:),Veg.(veg_names{k}));
end
check('the stand is perennial: it starts from a run''s end state handed in, else from the file''s mature stand over an established sward', ...
    isequal(P.PFT_perennial,true) && carried && isequal(squeeze(P.B_H(1,1,:)).',[0 3000 300 800 0 8000 0 0]) && ...
    isequal(squeeze(P.B_L(1,2,:)).',[0 0 140 93 0 0 0 0]) && P_veg.zatm == 26);
check('the stand grows as the legacy coupled forest runs: no sink limitation, heartwood turning over in 40 years', ...
    P.OPT_EnvLimitGrowth == 0 && P_veg.OPT_EnvLimitGrowth == 0 && ...
    abs(P.Wm_H(1)*40*365 - 1) < 1e-12 && isequal(P.Wm_L,[0 0]));
%% Crowns as the cell has them.  A crown without cover carries no vegetation:
%% T&C would grow it with no share of the soil water, it would starve and
%% its leftover reserve trip Nutrients_Available's e_relN guard every day.
%% Every cell is seeded from the wooded class seed, so a lawn is handed
%% trees and an empty open-grass crown: it keeps no trees and starts its
%% sward from the file's.  A stand that died keeps its remains.
wood = find(cover_tree == 1 & cover_grass == 0,1);
lawn = find(cover_tree == 0 & cover_grass > 0,1);
cell_params = @(ci,base) Load_Cell_Parameters( ...
    fullfile(code,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(code,'Parameters',[char(ci.MOD_PARAM.pft_override) '.m']), ...
    cell_context(base,ci));
P_wood = cell_params(CI(wood),ctx);
P_lawn = cell_params(CI(lawn),ctx);
sward = [0 0 140 93 0 0 0 0];
check('a wooded cell has no open grass: trees and the grass beneath them only', ...
    isequal(P_wood.ZR95_H,[1000 0]) && isequal(P_wood.ZR95_L,[250 0]) && ...
    ~any(P_wood.B_L(1,2,:)) && P_wood.Nreserve_L(1,2) == 0 && ...
    isequal(squeeze(P_wood.B_L(1,1,:)).',sward) && P_wood.zatm == 26);
check('a lawn has no trees: its open grass only, the reference height 2 m above the sward', ...
    isequal(P_lawn.ZR95_H,[0 0]) && isequal(P_lawn.ZR95_L,[0 250]) && ...
    ~any(P_lawn.B_H(1,1,:)) && ~any(P_lawn.B_L(1,1,:)) && P_lawn.Nreserve_H(1,1) == 0 && ...
    isequal(squeeze(P_lawn.B_L(1,2,:)).',sward) && abs(P_lawn.zatm - 2.05) < 1e-12, ...
    sprintf('zatm %g',P_lawn.zatm));
Seed = Veg;                                 % the wooded class seed's stand:
for k = 1:numel(veg_names)                  % its open-grass crown is empty
    Seed.(veg_names{k})(1,2,:) = 0;
end
Seed.B_H(1,1,:) = [0 1035 607 723 13 12618 1 0]; Seed.hc_H = [24 0];
Seed.B_L(1,1,:) = [1e-178 0 5e-265 7e-109 0 0 3.2 0];   % its shaded-out understory
ctx_seed_veg = ctx; ctx_seed_veg.Veg_out = Seed;
P_lawn_seeded = cell_params(CI(lawn),ctx_seed_veg);
P_wood_seeded = cell_params(CI(wood),ctx_seed_veg);
check('a lawn handed the wooded class seed''s stand keeps no trees and starts its sward from the file''s', ...
    ~any(P_lawn_seeded.B_H(1,1,:)) && P_lawn_seeded.Nreserve_H(1,1) == 0 && ...
    P_lawn_seeded.hc_H(1,1) == 0 && ~any(P_lawn_seeded.B_L(1,1,:)) && ...
    isequal(squeeze(P_lawn_seeded.B_L(1,2,:)).',sward) && ...
    P_lawn_seeded.Nreserve_L(1,2) == P_lawn.Nreserve_L(1,2) && P_lawn_seeded.Nreserve_L(1,2) > 0 && ...
    P_lawn_seeded.FNC_L(1,2) == P_lawn.FNC_L(1,2) && abs(P_lawn_seeded.zatm - 2.05) < 1e-12, ...
    sprintf('zatm %g',P_lawn_seeded.zatm));
check('the wooded cell carries its own stand on, the dead understory with its remains, not sown again', ...
    isequal(P_wood_seeded.B_H(1,1,:),Seed.B_H(1,1,:)) && isequal(P_wood_seeded.B_L(1,1,:),Seed.B_L(1,1,:)) && ...
    ~any(P_wood_seeded.B_L(1,2,:)) && P_wood_seeded.zatm == 26);
%% Hamden's own exchange constants (Fit_Exchange_Affinity): on the class
%% seed's ERW seed inputs, the exchanger that emerges at the representative
%% layers' pH holds their base saturation and KCl Al.
X = load(fullfile(code,'Inputs','Hamden','Calibration','ERW_seed_inputs_classseed_01.mat'));
Cx = X.ERW_TC_const; Cx.PH_in = chem.layers_pH;
Cx.bulk_pH_params.appelo = P.bulk_pH_params.appelo;
Cx.bulk_pH_params.Al_pin_logKso = P.bulk_pH_params.Al_pin_logKso;
[Bx,~,Wx,ECx] = ERW_MOD_PARAM_excel(X.t,X.TT_m,Cx,X.ERW_biogeochem,X.ERW_PARAM_input, ...
    X.opt_ERW_preprocessing);
Dx = ERW_Handoff_Diagnostic(ERW_Initialize_TC_Cation_Pools(X.B_out,Bx,ECx,X.TT_m.s(1,:)), ...
    Bx,Wx,ECx,X.cov_moll);
bs_x = Bx.f_Ca(1) + Bx.f_Mg(1) + Bx.f_K(1) + Bx.f_Na(1);
bs_m = chem.bs_Ca_fraction + chem.bs_Mg_fraction + chem.bs_K_fraction + chem.bs_Na_fraction;
delta = P.bulk_pH_params.appelo.logKH(1) - 5.9;
check('the site''s exchange constants hold the measured base saturation and KCl Al on the class seed''s soil solution', ...
    P.bulk_pH_params.appelo.alpha == 3.4 && isequal(P.bulk_pH_params.Al_solubility_pin,true) && ...
    max(abs(P.bulk_pH_params.appelo.logKH - ([5.9 4.8 10.47 10.786 16.47] + [1 1 2 2 3]*delta))) < 1e-12 && ...
    abs(X.ERW_PARAM_input.CEC_tot_in - chem.CEC_cmolc_kg) < 1e-9 && ...
    abs(bs_x - bs_m) < 1e-4 && abs(Bx.f_Al(1) - chem.Al_fraction) < 1e-4 && ...
    string(ERW_Handoff_Readiness(Dx).status) == "pass", ...
    sprintf('base saturation %.4f (measured %.4f), Al %.4f (%.4f)',bs_x,bs_m,Bx.f_Al(1),chem.Al_fraction));

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

function ctx = cell_context(ctx,ci)
%% A load context for cell ci: its own rows over the shared one.
ctx.MOD_PARAM_input = ci.MOD_PARAM;
ctx.ERW_PARAM_input = ci.ERW_PARAM;
end
