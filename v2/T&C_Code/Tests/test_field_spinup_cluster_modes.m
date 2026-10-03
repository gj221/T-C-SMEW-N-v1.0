function test_field_spinup_cluster_modes()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_field_spinup_cluster_modes
%%
%% The cell-parallel modes of Run_Field_Cell_Spinup and the validation
%% experiment, without spinning anything up:
%%   1  PlanOnly writes one plan row per cell and SLURM array-index lists
%%   2  CellIndex is validated before any work
%%   3  AssembleOnly refuses a missing or mis-labelled cell state and
%%      assembles a complete set into the file the transient loader reads
%%   4  Run_RF_Validation_Experiment scores a forest against fabricated
%%      truth states across fractions and seeds, reports base-saturation
%%      attainment on the measured cells, and scores an RF-path directory
%% Uses the Control 2 diagnostic grid (15 m, 50 % coverage rule) and the
%% real DEM tiles for the forest's covariates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_field_spinup_cluster_modes ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code, fullfile(v2,'spinup'), fullfile(v2,'interpolation'), fullfile(v2,'run'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
    function ok = raises(fn,id)
        ok = false;
        try; fn(); catch ME; ok = strcmp(ME.identifier,id); if ~ok; disp(ME.message); end; end
    end

shp = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp');
dem = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','DEM');
h = 15; cov = 0.5; feats = [1 2 3]; poly = 1;
setups = Build_Shapefile_Grid_Setups(shp,h,feats,cov); topo = setups(poly);
N = topo.N_cells;
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp,'s')); %#ok<NASGU>

%% template: baseline constants, two cells with measured base saturation
T = table(topo.cell_id,topo.cell_xc_projected_m,topo.cell_yc_projected_m, ...
    repmat(6.1,N,1),nan(N,1),false(N,1), ...
    repmat(0.10,N,1),repmat(0.35,N,1),nan(N,1),nan(N,1),false(N,1), ...
    repmat(0.025,N,1),nan(N,1),false(N,1), ...
    repmat(8.6,N,1),nan(N,1),false(N,1), ...
    nan(N,1),nan(N,1),nan(N,1),false(N,1), ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction', ...
    'measured_sand_fraction','measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC', ...
    'measured_bs_Ca_fraction','measured_bs_Mg_fraction','measured_bs_K_fraction', ...
    'use_measured_base_saturation'});
T.measured_bs_Ca_fraction(1:2) = [0.55; 0.60];
T.measured_bs_Mg_fraction(1:2) = [0.15; 0.12];
T.measured_bs_K_fraction(1:2)  = [0.03; 0.02];
T.use_measured_base_saturation(1:2) = true;
agro = fullfile(tmp,'agronomic.csv'); writetable(T,agro);
outdir = fullfile(tmp,'truth'); mkdir(outdir);
state_file = fullfile(outdir,'CellBGStates.mat');
common = {'FieldFeatureIndices',feats,'FieldPolygonIndex',poly,'GridCellSizeM',h, ...
    'MinCoverage',cov,'OutputFile',state_file};

%% --- 1. plan ------------------------------------------------------------
plan_file = Run_Field_Cell_Spinup(shp,agro,common{:},'PlanOnly',true);
plan = readtable(plan_file);
check('plan has one row per cell with the cell ids in grid order', ...
    height(plan) == N && isequal(plan.cell_id,topo.cell_id(:)) && all(plan.method == "full"));
all_idx = strtrim(fileread(fullfile(outdir,'array_indices_all.txt')));
check('the all-cells array list is the SLURM range 1-N', strcmp(all_idx,sprintf('1-%d',N)));
sub_idx = strtrim(fileread(fullfile(outdir,'array_indices_subsample.txt')));
check('per_cell plans no subsample', isempty(sub_idx));

%% --- 2. CellIndex validation ----------------------------------------------
check('CellIndex outside 1..N is refused before any work', ...
    raises(@() Run_Field_Cell_Spinup(shp,agro,common{:},'CellIndex',N+1), ...
    'Run_Field_Cell_Spinup:BadCellIndex'));

%% --- 3. assembly ----------------------------------------------------------
[~,CellInputs] = Run_Field_Cell_Spinup(shp,agro,common{:},'DryRun',true);
soil_mass = 1295*300/1000;
    function st = mock_state(c)
        st = struct();
        st.cell_id = topo.cell_id(c);
        st.cell_agronomic_identity = string(CellInputs(c).agronomic_identity);
        rng(1000+c);
        B = 1 + 50*rand(1,79);
        fr = [T.measured_bs_Ca_fraction(c) T.measured_bs_Mg_fraction(c) T.measured_bs_K_fraction(c)];
        if any(isnan(fr)); fr = [0.5 0.15 0.03]; end     % unmeasured cells: any plausible exchanger
        B(61) = fr(1)*8.6*soil_mass*0.01*(40.078/2);   % exchangeable Ca at this cell's measured target
        B(69) = fr(2)*8.6*soil_mass*0.01*(24.305/2);
        B(53) = fr(3)*8.6*soil_mass*0.01*39.098;
        %% the handoff exchanger seed the transient starts from
        st.spinup_erw_handoff_diagnostic = struct('f_CEC',[fr 0 0 1-sum(fr)]);
        st.B_out = B;                                   % one 79-pool row, as a real state
        st.O_out = 0.3*ones(1,20); st.NavlI_out = zeros(1,6);
        st.NupI_L_out = zeros(1,6); st.RexmyI_out = zeros(1,3);
        st.spinup_converged = true; st.spinup_complete = true; st.spinup_final_pass = true;
        st.spinup_erw_handoff_evaluated = true; st.spinup_erw_handoff_ready = true;
        st.spinup_nutrient_budget_vers = 12;
        b = struct('closed',true,'residual',0,'tolerance',1);
        st.spinup_nutrient_budgets = struct('N',b,'P',b,'K',b,'Ca',b,'Mg',b,'Si',b);
        st.spinup_biogeochemical_soil_mass_kg_m2 = soil_mass;
        st.spinup_water_loss_mean_l_m2_d = 0.5;
    end
    function write_cellstate(dirname,c,method)
        st = mock_state(c); cell_id = topo.cell_id(c); cell_index = c; %#ok<NASGU>
        save(fullfile(dirname,sprintf('cellstate_%04d.mat',topo.cell_id(c))), ...
            'st','method','cell_id','cell_index','-v7.3');
    end
for c = 1:N-1; write_cellstate(outdir,c,"full"); end
check('assembly refuses an incomplete set and names the missing cell', ...
    raises(@() Run_Field_Cell_Spinup(shp,agro,common{:},'AssembleOnly',true), ...
    'Run_Field_Cell_Spinup:CellStatesMissing'));
write_cellstate(outdir,N,"full");
bad = load(fullfile(outdir,sprintf('cellstate_%04d.mat',topo.cell_id(3))));
bad.st.cell_id = topo.cell_id(4);
save(fullfile(outdir,sprintf('cellstate_%04d.mat',topo.cell_id(3))),'-struct','bad');
check('assembly refuses a state labelled with another cell', ...
    raises(@() Run_Field_Cell_Spinup(shp,agro,common{:},'AssembleOnly',true), ...
    'Run_Field_Cell_Spinup:CellStateIdentity'));
write_cellstate(outdir,3,"full");
out_file = Run_Field_Cell_Spinup(shp,agro,common{:},'AssembleOnly',true);
S = load(out_file,'CellBGStates','rf_report');
check('a complete set assembles into one struct array with the cells in order', ...
    isfile(out_file) && numel(S.CellBGStates) == N && ...
    isequal([S.CellBGStates.cell_id].',topo.cell_id(:)) && S.rf_report.n_full == N);
states = Load_Field_Spinup_States(out_file,topo,string({CellInputs.agronomic_identity}));
check('the transient loader accepts the assembled file', numel(states) == N);

%% --- 4. validation experiment ----------------------------------------------
rfdir = fullfile(tmp,'rf'); mkdir(rfdir);
for c = 1:5; write_cellstate(rfdir,c,"rf_seed"); end
val = Run_RF_Validation_Experiment(shp,agro,'FieldFeatureIndices',feats, ...
    'FieldPolygonIndex',poly,'GridCellSizeM',h,'MinCoverage',cov,'FieldDEM',dem, ...
    'TruthStateFile',out_file,'Fractions',[0.3 0.5],'Seeds',[1 2],'KFold',3, ...
    'RFStatesDirectory',rfdir,'OutputDirectory',tmp);
R = val.regression;
check('one regression row per fraction x seed, all trained', ...
    height(R) == 4 && all(R.trained) && all(R.n_subsample + R.n_heldout == N));
check('held-out errors are finite for every conserved total', ...
    all(isfinite(R.heldout_rel_err_POC_gCm2)) && all(isfinite(R.heldout_rel_err_total_P_gPm2)));
check('base-saturation attainment is exact on the fabricated truth', ...
    height(val.base_saturation) == 2 && max(abs(val.base_saturation.rel_err_Ca_Mg_K(:))) < 1e-9);
check('the RF-path table scores the cells that have an RF state', ...
    istable(val.rf_path) && height(val.rf_path) == 5 && all(val.rf_path.handoff_ready) && ...
    all(abs(val.rf_path.rel_err_POC_gCm2) < 1e-12));
check('CSV and figure written', isfile(val.files.regression) && isfile(val.files.figure) && ...
    isfile(val.files.base_saturation) && isfile(val.files.rf_path));

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_field_spinup_cluster_modes: %d check(s) failed',n_fail); end
end
