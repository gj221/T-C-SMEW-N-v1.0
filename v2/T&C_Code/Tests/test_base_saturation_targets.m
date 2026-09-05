function test_base_saturation_targets()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_base_saturation_targets
%%
%% Measured exchangeable base saturation enters the model as spin-up
%% observation targets for the exchangeable pools B(61) Ca, B(69) Mg,
%% B(53) K.  Checks the chain end to end at the unit level:
%%   1  BG_Observation_Targets converts cmol_c kg-1 -> g m-2 with the
%%      equivalent weights and the biogeochemical soil mass
%%   2  Initialize_BG_State_From_Observations writes exactly those pools
%%      and touches nothing else
%%   3  Resolve_Cell_Agronomic_Inputs turns fractions x CEC into the
%%      exchangeable stocks, is optional, validates the fractions, and
%%      appends to the identity only where used
%%   4  Configure_Cell_Agronomic_Inputs carries them into the per-cell
%%      MOD_PARAM row that BG_Observation_Targets reads
%%   5  the template carries the optional columns and Apply_Kriged_Cell_
%%      Values treats the three fractions as one jointly-accepted toggle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_base_saturation_targets ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code, fullfile(v2,'spinup'), fullfile(v2,'interpolation'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% --- 1. unit conversion ----------------------------------------------
MOD = table(1295,300,2.0,0.5,0.3, ...
    'VariableNames',{'rsd','Zbio','InitExchCa_cmolc_kg','InitExchMg_cmolc_kg','InitExchK_cmolc_kg'});
tg = BG_Observation_Targets(MOD);
soil_mass = 1295*300/1000;                       % kg m-2
check('Ca: cmol_c/kg -> g/m2 uses 40.078/2 g per mol_c', ...
    abs(tg.exch_Ca_gCam2 - 2.0*0.01*(40.078/2)*soil_mass) < 1e-9, ...
    sprintf('%.6g',tg.exch_Ca_gCam2));
check('Mg: uses 24.305/2', abs(tg.exch_Mg_gMgm2 - 0.5*0.01*(24.305/2)*soil_mass) < 1e-9);
check('K: monovalent, uses 39.098/1', abs(tg.exch_K_gKm2 - 0.3*0.01*39.098*soil_mass) < 1e-9);
MOD0 = MOD(:,{'rsd','Zbio'});
tg0 = BG_Observation_Targets(MOD0);
check('absent columns produce no target (spin-up exchanger stands)', ...
    ~isfield(tg0,'exch_Ca_gCam2') && ~isfield(tg0,'exch_K_gKm2'));

%% --- 2. pool mapping ---------------------------------------------------
B = 0.5 + rand(1,79);
Bn = Initialize_BG_State_From_Observations(B, ...
    struct('exch_Ca_gCam2',77.7,'exch_Mg_gMgm2',23.6,'exch_K_gKm2',45.5));
check('targets land on B(61) Ca, B(69) Mg, B(53) K exactly', ...
    Bn(61) == 77.7 && Bn(69) == 23.6 && Bn(53) == 45.5);
untouched = setdiff(1:79,[61 69 53]);
check('no other pool is changed by an exchangeable-base target', ...
    isequal(Bn(untouched),B(untouched)));

%% --- 3. resolve -----------------------------------------------------
topo = struct('N_cells',3,'cell_id',(1:3)','cell_xc_projected_m',[0;10;20], ...
    'cell_yc_projected_m',[0;0;0]);
T = base_table(3);
T.measured_bs_Ca_fraction = [0.55;0.60;NaN];
T.measured_bs_Mg_fraction = [0.15;0.12;NaN];
T.measured_bs_K_fraction  = [0.03;0.02;NaN];
T.use_measured_base_saturation = [true;true;false];
R = Resolve_Cell_Agronomic_Inputs(topo,T);
check('exchangeable stock = fraction x resolved CEC [cmol_c/kg]', ...
    abs(R.exch_Ca_cmolc_kg(1) - 0.55*8.6) < 1e-12 && ...
    abs(R.exch_K_cmolc_kg(2) - 0.02*8.6) < 1e-12);
check('a cell without the toggle carries NaN and source "none"', ...
    isnan(R.exch_Ca_cmolc_kg(3)) && R.base_saturation_source(3) == "none");
check('identity appended only where used', ...
    contains(R.input_identity(1),'exch=') && ~contains(R.input_identity(3),'exch='));
R0 = Resolve_Cell_Agronomic_Inputs(topo,base_table(3));
check('a template without the columns still resolves (all NaN, "none")', ...
    all(isnan(R0.exch_Ca_cmolc_kg)) && all(R0.base_saturation_source == "none"));
check('the identity of an unmeasured cell is unchanged by the new columns', ...
    R0.input_identity(3) == R.input_identity(3));
Tbad = T; Tbad.measured_bs_Ca_fraction(1) = 0.9;      % 0.9+0.15+0.03 > 1
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,Tbad), ...
    'Resolve_Cell_Agronomic_Inputs:InvalidBaseSaturation');
check('fractions that leave no room for acid saturation are rejected', true);

%% --- 4. configure ---------------------------------------------------
MODc = table(6.1,0.10,0.35,0.025,300,1295,'VariableNames', ...
    {'PH','Psan','Pcla','Porg','Zbio','rsd'});
ERWc = table(8.6,'VariableNames',{'CEC_tot_in'});
CI = Configure_Cell_Agronomic_Inputs(MODc,ERWc,R);
check('measured cells get the three InitExch* MOD_PARAM columns', ...
    ismember('InitExchCa_cmolc_kg',CI(1).MOD_PARAM.Properties.VariableNames) && ...
    abs(CI(1).MOD_PARAM.InitExchCa_cmolc_kg - 0.55*8.6) < 1e-12);
check('an unmeasured cell gets none (targets absent, spin-up exchanger stands)', ...
    ~ismember('InitExchCa_cmolc_kg',CI(3).MOD_PARAM.Properties.VariableNames));
tg3 = BG_Observation_Targets(CI(1).MOD_PARAM);
check('the per-cell row feeds BG_Observation_Targets directly', ...
    isfield(tg3,'exch_Ca_gCam2') && tg3.exch_Ca_gCam2 > 0);

%% --- 5. template + merge -----------------------------------------------
tmp = tempname; mkdir(tmp);
shp = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp');
if isfile(shp)
    files = Write_Field_Agronomic_Templates(shp,tmp,[],1);
    Tt = readtable(files{1},'VariableNamingRule','preserve');
    check('the template carries the optional base-saturation columns', ...
        all(ismember({'measured_bs_Ca_fraction','measured_bs_Mg_fraction', ...
        'measured_bs_K_fraction','use_measured_base_saturation'}, ...
        Tt.Properties.VariableNames)));
    n = height(Tt);
    Tt.baseline_bulk_pH(:) = 6; Tt.baseline_sand_fraction(:) = 0.1;
    Tt.baseline_clay_fraction(:) = 0.3; Tt.baseline_OM_fraction(:) = 0.02;
    Tt.baseline_CEC_cmolc_kg(:) = 8;
    tfile = fullfile(tmp,'template.csv'); writetable(Tt,tfile);
    ids = Tt.cell_id(1:4);
    K = table(ids,[0.5;0.5;0.5;0.5],[0.01;0.01;0.01;0.01],true(4,1), ...
        [0.1;0.1;0.1;0.1],[0.01;0.01;0.01;3],true(4,1), ...
        [0.02;0.02;0.02;0.02],[0.001;0.001;0.001;0.001],[true;true;false;true], ...
        'VariableNames',{'cell_id','bs_Ca_fraction_predicted','bs_Ca_fraction_std', ...
        'bs_Ca_fraction_in_sample_hull','bs_Mg_fraction_predicted','bs_Mg_fraction_std', ...
        'bs_Mg_fraction_in_sample_hull','bs_K_fraction_predicted','bs_K_fraction_std', ...
        'bs_K_fraction_in_sample_hull'});
    [out_file,~] = Apply_Kriged_Cell_Values(tfile,K,'MaxRelativeStd',0.35, ...
        'OutputFile',fullfile(tmp,'merged.csv'));
    M = readtable(out_file,'VariableNamingRule','preserve');
    r = @(i) find(M.cell_id == ids(i));
    check('all three fractions accepted -> the joint toggle is on', ...
        M.use_measured_base_saturation(r(1)) && M.use_measured_base_saturation(r(2)));
    check('one fraction out of hull -> the joint toggle stays off', ...
        ~M.use_measured_base_saturation(r(3)));
    check('one fraction over the std gate -> the joint toggle stays off', ...
        ~M.use_measured_base_saturation(r(4)));
    check('measured values are written whether or not the toggle turns on', ...
        abs(M.measured_bs_Ca_fraction(r(3)) - 0.5) < 1e-12);
else
    fprintf('  SKIP  template/merge checks (control-field shapefile not present)\n');
end

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_base_saturation_targets: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function T = base_table(n)
T = table((1:n)',[0;10;20],zeros(n,1), ...
    repmat(6.1,n,1),nan(n,1),false(n,1), ...
    repmat(0.10,n,1),repmat(0.35,n,1),nan(n,1),nan(n,1),false(n,1), ...
    repmat(0.025,n,1),nan(n,1),false(n,1), ...
    repmat(8.6,n,1),nan(n,1),false(n,1), ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction', ...
    'measured_sand_fraction','measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC'});
end

%% =======================================================================
function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end
