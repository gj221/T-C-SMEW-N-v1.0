function test_spinup_handoff_only()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_spinup_handoff_only
%%
%% The ERW charge-balance handoff of a spun-up state, and the handoff a
%% forest-reconstructed cell of the flux-tracking field spin-up runs
%% (Run_Automatic_Site_Spinup with SPINUP_HANDOFF_ONLY):
%%   1  the handoff diagnostic balances the species the state holds: v0.9
%%      speciates the pore water with activity-conditioned constants, and
%%      re-speciating its Al with the thermodynamic constants reported the
%%      correction as an imbalance (0.7 umolc/L on Control 2 cell 1, 75 on
%%      cell 17) and failed every acidic Bruce Farm cell; a held species
%%      that does not close still fails the gate
%%   2  on a real Control 2 diagnostic cell with the Kantola case and
%%      forcing, no spin-up phase runs; the ERW seed is built from the
%%      cell's own inputs and passes the gate; the state keeps its pools and
%%      soil moisture, carries the template's convergence record and
%%      reconstruction provenance, hands off at the transient start and
%%      validates for the transient
%%   3  a state without plant reserves or without a converged record is
%%      refused
%% The template is the archived converged Kantola state with a fabricated
%% convergence record, as Reconstruct_Cell_Pools_RF hands a replayed one on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_spinup_handoff_only ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code,fullfile(code,'Parameters'),fullfile(code,'Inputs'),fullfile(v2,'spinup'));
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
old_dir = pwd; restore_dir = onCleanup(@() cd(old_dir));

%% 1  the balance over the held species (umol/L; Al species on the
%%    Al_CEC_scale basis, distributed as conditioned constants put them)
E = struct('H',2.45,'Ca',530,'Mg',310,'Na',48,'K',96,'NH4',21, ...
    'Al',3159.48,'AlOH',6497.28,'AlOH2',14002.1,'AlOH3',1000,'AlOH4',297.957, ...
    'Al_w',24956.817,'k_w',1e-2,'HCO3',74.1438,'CO3',9.75e-4,'NO3',402,'PO4',1.2);
C = struct('Al_CEC_scale',1000,'K1',10^-4.97*1e6,'K2',10^-5.1*1e6, ...
    'K3',10^-6.7*1e6,'K4',10^-6.2*1e6);
E.An = E.H + 2*E.Ca + 2*E.Mg + E.Na + E.K + E.NH4 + ...
    (3*E.Al + 2*E.AlOH + E.AlOH2)/C.Al_CEC_scale - ...
    (E.k_w/E.H + E.HCO3 + 2*E.CO3 + E.NO3 + E.PO4 + E.AlOH4/C.Al_CEC_scale);
charge_flagged = @(D) any(contains(ERW_Handoff_Readiness(D).warnings,"charge residual"));
D0 = ERW_Handoff_Diagnostic(zeros(1,79),E,struct(),C,1);
check('the handoff gate balances the species the state holds', ...
    abs(D0.charge_residual_umolc_L) <= 1e-9 && ~charge_flagged(D0), ...
    sprintf('%.3g umolc/L',D0.charge_residual_umolc_L));
E_off = E; E_off.NO3 = E_off.NO3 + 1;
D1 = ERW_Handoff_Diagnostic(zeros(1,79),E_off,struct(),C,1);
check('a held species that does not close fails the gate', ...
    abs(D1.charge_residual_umolc_L + 1) <= 1e-9 && charge_flagged(D1), ...
    sprintf('%.3g umolc/L',D1.charge_residual_umolc_L));

%% 2  the handoff-only path on Control 2 cell 1
shp = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','Whole farm','AdamWilbourne_Control_Fields.shp');
agro = fullfile(code,'Inputs','Field_Spinup','diagnostic_Control2_15m','Cell_Agronomic_merged.csv');
cols = struct('forcing_file',"Data_Kantola.mat",'meteo_file',"Data_Kantola.mat", ...
    'meteo_start',"2009-01-01",'meteo_end',"2015-12-31", ...
    'transient_start',"2016-01-01",'transient_end',"2020-12-30");
tmp = tempname; mkdir(tmp); cleanup = onCleanup(@() rmdir(tmp,'s'));
[~,CI] = Run_Field_Cell_Spinup(shp,agro,'FieldFeatureIndices',[1 2 3],'FieldPolygonIndex',1, ...
    'GridCellSizeM',15,'MinCoverage',0.5,'Method','per_cell','MODParamColumns',cols, ...
    'OutputFile',fullfile(tmp,'CellBGStates.mat'),'DryRun',true);
c = 1; cell = CI(c);

K = load(fullfile(code,'Inputs','Data_Kantola_BG.mat'));
st = struct('B_out',K.B_out(end,:),'O_out',0.97*K.O_out(end,:),'NavlI_out',K.NavlI_out(end,:), ...
    'NupI_L_out',K.NupI_L_out,'RexmyI_out',K.RexmyI_out(end,:));
for nut = {'N','P','K','Ca','Mg','Si'}
    st.([nut{1} 'reserve_L_out']) = zeros(1,cell.MOD_PARAM.cc);
    st.([nut{1} 'reserve_H_out']) = zeros(1,cell.MOD_PARAM.cc);
end
st.spinup_converged = true; st.spinup_complete = true; st.spinup_final_pass = true;
st.spinup_cycles_completed = 42; st.spinup_ring_streak_achieved = 1;
st.spinup_nutrient_budget_vers = 12;
b = struct('closed',true,'residual',0,'tolerance',1);
st.spinup_nutrient_budgets = struct('N',b,'P',b,'K',b,'Ca',b,'Mg',b,'Si',b);
st.spinup_biogeochemical_soil_mass_kg_m2 = 1295*300/1000;
st.reconstruction_method = "rf"; st.reconstruction_template_cell_id = 7;
st.reconstruction_rf_cv_rel_error = 0.05*ones(1,8);
st.cell_id = c; st.cell_agronomic_identity = cell.agronomic_identity;

    function f = seed_file(name,S)
        f = fullfile(tmp,[name '.mat']); save(f,'-struct','S','-v7.3');
    end
    function settings = handoff_settings(seed,out)
        settings = struct('NN_case_study_ERW_MOD',5,'RUN_MODE','spinup','RUN_CHEMISTRY','erw', ...
            'RUN_DOMAIN','single','SAVE_RESULTS',0,'MOD_PARAM_override',cell.MOD_PARAM, ...
            'ERW_PARAM_override',cell.ERW_PARAM,'CELL_ID',c, ...
            'CELL_AGRONOMIC_IDENTITY',cell.agronomic_identity, ...
            'BG_OUTPUT_FILE',fullfile(tmp,out),'SPINUP_SEED_FILE',seed, ...
            'SPINUP_HANDOFF_ONLY',true,'SPINUP_PROVENANCE','test-handoff');
    end

settings = handoff_settings(seed_file('seed',st),'handoff_BG.mat');
log_text = evalc('Run_Automatic_Site_Spinup(settings);');
check('no spin-up phase runs', ~contains(log_text,'=== Phase 1') && ...
    ~contains(log_text,'=== Phase 2') && ~contains(log_text,'Phase 3 loop') && ...
    contains(log_text,'=== Phase 3: handoff only, no convergence loop'));
check('the ERW seed of the cell''s own inputs passes the handoff gate', ...
    contains(log_text,'ERW handoff readiness: pass'), ...
    char(regexp(log_text,'ERW handoff readiness:[^\n]*','match','once')));
R = load_bg_state(settings.BG_OUTPUT_FILE);
ok_valid = true;
try; Validate_BG_State_For_Transient(R,settings.BG_OUTPUT_FILE); catch; ok_valid = false; end
check('the handed-on state validates for the transient', ok_valid);
check('its charge balance closes', R.spinup_erw_handoff_ready && ...
    abs(R.spinup_erw_handoff_charge_residual_umolc_L) <= 1e-6, ...
    sprintf('%.3g umolc/L',R.spinup_erw_handoff_charge_residual_umolc_L));
check('pools and soil moisture are the reconstructed ones', ...
    isequal(R.O_out,st.O_out) && max(abs(R.B_out - st.B_out)) <= 1e-9*max(abs(st.B_out)));
check('the template''s record and the reconstruction provenance are kept', ...
    R.spinup_converged && R.spinup_cycles_completed == 42 && ...
    R.reconstruction_method == "rf" && R.reconstruction_template_cell_id == 7);
check('the exchanger seed of the cell''s own inputs is recorded', ...
    isfield(R.spinup_erw_handoff_diagnostic,'f_CEC') && ...
    abs(sum(R.spinup_erw_handoff_diagnostic.f_CEC) - 1) < 1e-3);
check('it hands off at the transient start', abs(R.spinup_handoff_dn - datenum(2016,1,1)) < 1/48);

%% 3  refusals
no_reserve = rmfield(st,'Nreserve_L_out');
check('a state without plant reserves is refused', raises(@() Run_Automatic_Site_Spinup( ...
    handoff_settings(seed_file('no_reserve',no_reserve),'x_BG.mat')), ...
    'Run_BG_Spinup:HandoffOnlyNeedsState'));
unconverged = st; unconverged.spinup_converged = false;
check('a state without a converged record is refused', raises(@() Run_Automatic_Site_Spinup( ...
    handoff_settings(seed_file('unconverged',unconverged),'y_BG.mat')), ...
    'Run_BG_Spinup:HandoffOnlyUnconverged'));

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_spinup_handoff_only: %d check(s) failed',n_fail); end
end
