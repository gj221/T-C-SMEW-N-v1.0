% =========================================================================
% Run_BG_Spinup.m  |  3-phase legacy TC-SMEW spin-up
%
%   Phase 1: T&C hydrology + vegetation (OPT_SoilBiogeochemistry=0, OPT_pH=0)
%            Run each spin-up cycle year and extract annual-average fluxes.
%
%   Phase 2: Standalone BIOGEO_UNIT daily loop with constant Phase-1 fluxes.
%            Iterate until all BG pools converge.
%
%   Phase 3: Coupled MAIN_FRAME (OPT_SoilBiogeochemistry=1, OPT_ERW=0) looped
%            over the transient period until dynamic equilibrium.
%            One final cycle with OPT_pH forced on for ERW charge-balance check.
%
%   Checkpoint files allow automatic resume at any phase boundary.
%
%   Field spin-up (Lian et al. 2025 flux tracking, Run_Field_Cell_Spinup):
%     LATERAL_REPLAY_FILE  a tracked cell replays, in every phase, the
%                          lateral exchange recorded for it in the 2-D field
%                          run (water in Phases 1 and 3, the mean daily
%                          solute map in Phase 2, both in Phase 3 and the
%                          handoff extension), with its terrain radiation.
%     SPINUP_HANDOFF_ONLY  a reconstructed state (random forest) runs only
%                          the observation constraint and the ERW handoff.
% =========================================================================

%% Checkpoint file names
p1_chk        = spinup_phase_checkpoint_filename(bg_file, 'p1_avgs');
p2_inprog_chk = spinup_phase_checkpoint_filename(bg_file, 'p2_inprog');
p2_chk        = spinup_phase_checkpoint_filename(bg_file, 'p2_standalone');
p3_chk        = spinup_phase_checkpoint_filename(bg_file, 'p3_loop');
%%% The run these checkpoints belong to (Spinup_Provenance: cell, inputs,
%%% configuration and seed).  Every checkpoint and the final state carry it
%%% and a checkpoint of another run is not resumed, so an interrupted task
%%% resumes its own work and any change of inputs or seed starts clean.
%%% Empty for a plain site spin-up, which keeps the earlier behaviour.
spinup_provenance = '';
if exist('RUN_SETTINGS','var') && isstruct(RUN_SETTINGS) && ...
        isfield(RUN_SETTINGS,'SPINUP_PROVENANCE')
    spinup_provenance = char(string(RUN_SETTINGS.SPINUP_PROVENANCE));
end
%%% The CO2 record (Run_Automatic_Site_Spinup sets the site's; v0.9's
%%% historical Ca_Data.mat otherwise).
if ~exist('co2_record_file','var') || isempty(co2_record_file)
    co2_record_file = 'Ca_Data.mat';
end
%%% A tracked cell's lateral record (Run_Field_Lateral_Tracking).  The 1-D
%%% driver replays it through Lateral_Replay_*; without one nothing changes.
LATERAL_REPLAY = struct([]);
if exist('RUN_SETTINGS','var') && isstruct(RUN_SETTINGS) && ...
        isfield(RUN_SETTINGS,'LATERAL_REPLAY_FILE') && ...
        strlength(strtrim(string(RUN_SETTINGS.LATERAL_REPLAY_FILE))) > 0
    lateral_record_file = char(string(RUN_SETTINGS.LATERAL_REPLAY_FILE));
    tmp = load(lateral_record_file,'record');
    if ~isfield(tmp,'record') || ~isstruct(tmp.record) || ~isscalar(tmp.record)
        error('Run_BG_Spinup:BadLateralRecord','%s holds no lateral record.', ...
            lateral_record_file);
    end
    LATERAL_REPLAY = tmp.record;
    fprintf('Lateral replay: %s (%d steps, %d days)\n',lateral_record_file, ...
        numel(LATERAL_REPLAY.hour_dn),numel(LATERAL_REPLAY.day_dn));
    clear tmp lateral_record_file
end
lateral_replay_spinup = ~isempty(fieldnames(LATERAL_REPLAY));
%%% The handoff alone: the state handed in (SPINUP_SEED_FILE) already stands
%%% at the transient start -- the random-forest reconstruction of a cell the
%%% flux-tracking spin-up did not replay -- so only the observation
%%% constraint and the ERW charge-balance handoff run, under the cell's own
%%% inputs, and the state keeps the convergence record of its template.
spinup_handoff_only = exist('RUN_SETTINGS','var') && isstruct(RUN_SETTINGS) && ...
    isfield(RUN_SETTINGS,'SPINUP_HANDOFF_ONLY') && ...
    logical(RUN_SETTINGS.SPINUP_HANDOFF_ONLY);
if spinup_handoff_only && ~(exist('bg_state','var') && isstruct(bg_state) && ...
        isfield(bg_state,'B_out'))
    error('Run_BG_Spinup:HandoffOnlyNeedsState', ...
        'SPINUP_HANDOFF_ONLY needs the state to hand off (SPINUP_SEED_FILE).');
end
if spinup_handoff_only && lateral_replay_spinup
    error('Run_BG_Spinup:HandoffOnlyReplay', ...
        'A handoff-only state has no lateral replay to run.');
end
%%% Convergence is one drift limit stated per complete meteorological ring.
%%% Each phase applies it over its own step: Phase 3 advances a whole ring and
%%% uses the limit directly; Phase 2 advances one mean-forcing year and uses the
%%% same limit divided by the ring length.  Both therefore enforce an identical
%%% per-year drift whatever number of forcing years a scenario supplies.  The
%%% limit is relative, taken against a pool scale floored at the resolved mass
%%% below.
P3_RING_REL_TOL = 0.01; % [-] fractional pool drift per complete meteorological ring
P2_ANNUAL_REL_TOL = P3_RING_REL_TOL/n_spinup_years; % [-] the same limit over one year
P3_SLOW_RING_REL_TOL = P3_RING_REL_TOL; % slow pools share the dynamic-pool limit
%%% A pool holding less than this is below the mass the spin-up resolves, so
%%% it is judged against the floor rather than against itself.  This is what
%%% keeps a near-zero pool from failing on a negligible change, and it is
%%% ring-length independent, so it needs no per-phase scaling.
P3_POOL_SCALE_FLOOR = 10; % [g element m-2] resolved mass scale
P2_POOL_SCALE_FLOOR = P3_POOL_SCALE_FLOOR;
P3_POOL_NEG_TOL = 1e-6; % [g element m-2], below the resolved convergence scale
P2_POOL_NEG_TOL = P3_POOL_NEG_TOL;
P2_STREAK_REQ = n_spinup_years; % one complete meteorological ring of consecutive passing mean-forcing years
P3_RING_STREAK_REQ = 1; % one complete meteorological ring
%%% Pools the Phase-3 ring gate does not hold to stationarity.  Occluded P
%%% (46) has no return flux: dB46/dt = K3p*B45.  Secondary P (45) turns over
%%% in 1/K2p = 37 years toward mineral P and 1/K3p ~ 82,000 years toward
%%% occluded P, and a field's P surplus accumulates in it ring after ring
%%% (+3.8 g P m-2 per ring on a Control 2 cell after 234 coupled years), so
%%% it cannot settle on a spin-up's timescale.  Phase 2 still converges it
%%% under the mean forcing, cheaply, which hands Phase 3 a near-equilibrium
%%% split between mineral and secondary P; mineral P (43) stays gated.
P3_UNGATED_POOLS = [45 46];
%%% A cycle year of a vegetated case whose canopy never reaches this leaf
%%% area is bare soil: Phases 2 and 3 would converge the soil of a bare
%%% field (the pre-merge reserve limiter left LAI ~0.003 and passed).
P1_MIN_PEAK_LAI = 0.5; % [m2 m-2]
P3_NUTRIENT_BUDGET_ABS_TOL = 0.02; % [g element m-2] per meteorological year
P3_NUTRIENT_BUDGET_REL_TOL = 1e-3; % [-]
STANDING_NUTRIENT_STATE_VERS = 7; % v7: v0.9 package u (legacy plant-reserve dilution; Phase 2 forced by the mean anaerobic fraction)
P3_NUTRIENT_BUDGET_VERS = 12;
%%% v12: maize leaf Ca/Mg targets and grain Ca corrected (see
%%% Veg_Stoichiometric_Parameter.m).
%%% Rejects Phase 1/2 v5 and Phase 3 v11 checkpoints.
%%% The observed MAP dose remains split equally between spin-up maize years;
%%% the transient retains the full observed 2022 dose.
%%% Phase 3 applies no total-P boundary. A boundary would hide a management
%%% imbalance by creating or removing P. Non-occluded P other than secondary P
%%% (P3_UNGATED_POOLS) must settle across the complete rotation; full total P
%%% remains diagnostic because pool 46 has no return flux. Total P is screened
%%% in Validate_BG_Pool_Physicality.

%%% Rebuilt from the workbook row here, before any parameter file has run, so a
%%% checkpoint carrying a differently constrained pool state can be rejected
%%% rather than silently resumed on.  MOD_PARAM re-derives the identical struct
%%% from the same row.  Phase 1 carries no biogeochemistry and so is unaffected
%%% by the constraint; only the Phase 2 and Phase 3 checkpoints are gated on it.
BG_SPINUP_CONSTRAINT = BG_Spinup_Constraint(MOD_PARAM_input);
%%% The same measured slow stocks are a one-time state constraint at the real
%%% transient start. They are not reapplied inside the recurrent ring.
BG_TRANSIENT_INITIAL_CONSTRAINT = BG_Spinup_Constraint(MOD_PARAM_input, true);
%%% v1.0: the water the dissolved pools leach with (Leaching_Water_Basis), as
%%% the frame reads it from the same row; Phase 2 leaches as the frame does,
%%% also when Phase 1 is taken from its checkpoint and the frame has not run.
leach_layer = strcmp(Leaching_Water_Basis(MOD_PARAM_input),'layer');

[p1_done, p1_avgs] = load_compatible_p1_checkpoint(p1_chk, spinup_cycle_years, STANDING_NUTRIENT_STATE_VERS);
[p2_done, p2_bg_state] = load_compatible_p2_checkpoint(p2_chk, spinup_cycle_years, BG_SPINUP_CONSTRAINT, STANDING_NUTRIENT_STATE_VERS);
[p3_done, p3_resume_state] = load_compatible_bg_checkpoint(p3_chk, spinup_cycle_years, BG_SPINUP_CONSTRAINT, 'p3_observation_constraint');
if p1_done && ~spinup_provenance_matches(p1_avgs, spinup_provenance)
    fprintf('  Ignoring a Phase 1 checkpoint written for another run.\n');
    p1_done = false;
    p1_avgs = struct();
end
if p2_done && ~spinup_provenance_matches(p2_bg_state, spinup_provenance)
    fprintf('  Ignoring a Phase 2 checkpoint written for another run.\n');
    p2_done = false;
    p2_bg_state = struct();
end
if p3_done && ~spinup_provenance_matches(p3_resume_state, spinup_provenance)
    fprintf('  Ignoring a Phase 3 checkpoint written for another run.\n');
    p3_done = false;
    p3_resume_state = struct();
end

if p3_done && isfield(p3_resume_state, 'p3_P_boundary_constraint')
    fprintf('  Ignoring Phase 3 checkpoint: contains deprecated P boundary constraint.\n');
    p3_resume_state = struct();
    p3_done = false;
end
if p3_done && (~isfield(p3_resume_state, 'p3_transient_observation_constraint') || ...
        ~isequaln(p3_resume_state.p3_transient_observation_constraint, BG_TRANSIENT_INITIAL_CONSTRAINT))
    fprintf('  Ignoring Phase 3 checkpoint: transient observation constraint changed.\n');
    p3_resume_state = struct();
    p3_done = false;
end
if p3_done && (~isfield(p3_resume_state, 'p3_nutrient_budget_vers') || ...
        ~isnumeric(p3_resume_state.p3_nutrient_budget_vers) || ...
        ~isscalar(p3_resume_state.p3_nutrient_budget_vers) || ...
        ~isfinite(p3_resume_state.p3_nutrient_budget_vers) || ...
        p3_resume_state.p3_nutrient_budget_vers ~= P3_NUTRIENT_BUDGET_VERS)
    fprintf('  Ignoring Phase 3 checkpoint: nutrient budget schema changed; restarting from Phase 2 state.\n');
    p3_resume_state = struct();
    p3_done = false;
end

if p3_done
    p3_done = isfield(p3_resume_state, 'p3_loop') && ...
        isscalar(p3_resume_state.p3_loop) && isfinite(p3_resume_state.p3_loop) && ...
        p3_resume_state.p3_loop >= 0 && ...
        isfield(p3_resume_state, 'p3_pass_by_ring') && ...
        numel(p3_resume_state.p3_pass_by_ring) == n_spinup_years && ...
        isfield(p3_resume_state, 'p3_pool_ring_scale') && ...
        (isempty(p3_resume_state.p3_pool_ring_scale) || ...
        (isnumeric(p3_resume_state.p3_pool_ring_scale) && ...
        isreal(p3_resume_state.p3_pool_ring_scale) && ...
        numel(p3_resume_state.p3_pool_ring_scale) == 79 && ...
        all(isfinite(p3_resume_state.p3_pool_ring_scale(:))) && ...
        all(p3_resume_state.p3_pool_ring_scale(:) >= P3_POOL_SCALE_FLOOR))) && ...
        isfield(p3_resume_state, 'p3_ring_streak') && ...
        isnumeric(p3_resume_state.p3_ring_streak) && ...
        isscalar(p3_resume_state.p3_ring_streak) && ...
        isfinite(p3_resume_state.p3_ring_streak) && ...
        p3_resume_state.p3_ring_streak >= 0;
    if p3_done
        p3_done = spinup_nutrient_budget_resumed(p3_resume_state);
        if ~p3_done
            fprintf('  Ignoring Phase 3 checkpoint: nutrient budget incomplete or not closed.\n');
            p3_resume_state = struct();
        end
    end
    if ~p3_done
        fprintf('  Ignoring incomplete Phase 3 checkpoint; restarting from an earlier phase.\n');
        p3_resume_state = struct();
    end
end
if spinup_handoff_only
    p1_done = false; p2_done = false; p3_done = false;
end
%%% Warm start (SPINUP_WARM_START): Phase 3 begins from the converged state
%%% of an earlier spin-up of this cell (SPINUP_SEED_FILE) under the cell's
%%% present inputs -- the cold field spin-up, say, once site observations
%%% replace its provisional soil inputs.  Phases 1 and 2, which only bring a
%%% cold soil near equilibrium, are skipped; Phase 3 converges its ring as
%%% always and the handoff is built from the new inputs.
spinup_warm_start = exist('RUN_SETTINGS','var') && isstruct(RUN_SETTINGS) && ...
    isfield(RUN_SETTINGS,'SPINUP_WARM_START') && logical(RUN_SETTINGS.SPINUP_WARM_START);
if spinup_warm_start
    if spinup_handoff_only
        error('Run_BG_Spinup:WarmHandoffOnly', ...
            'A warm start runs Phase 3; it cannot also be handoff-only.');
    end
    warm_fields = {'B_out','O_out','NavlI_out','RexmyI_out','NupI_L_out'};
    if ~(exist('bg_state','var') && isstruct(bg_state) && all(isfield(bg_state,warm_fields)))
        error('Run_BG_Spinup:WarmStartNeedsState', ...
            'SPINUP_WARM_START needs the earlier state (SPINUP_SEED_FILE) with %s.', ...
            strjoin(warm_fields,', '));
    end
    if ~p3_done
        p1_done = true; p2_done = true;
        p2_bg_state = bg_state;
        for warm_name = {'B_out','O_out','NavlI_out','RexmyI_out'}
            p2_bg_state.(warm_name{1}) = bg_state.(warm_name{1})(end,:);
        end
        fprintf('Warm start: Phase 3 from the earlier state handed in (Phases 1-2 skipped)\n');
    end
    clear warm_fields warm_name
end

%% =========================================================================
%% Phase 1: T&C-only run; extract annual-average fluxes
%% =========================================================================
if ~spinup_handoff_only && ~p1_done && ~p2_done && ~p3_done
    fprintf('\n=== Phase 1: T&C-only run (%d years) ===\n', n_spinup_years);
    %%% Every Phase-1 year lies in a tracked cell's record.
    if lateral_replay_spinup
        LATERAL_REPLAY(1).required = true;
    end

    p1_n          = 0;
    p1_Ta_acc     = 0;
    p1_Tdp_acc    = 0;
    p1_Se_acc     = 0;
    p1_Psi_acc    = 0;
    p1_VSUM_acc   = 0;
    p1_VTSUM_acc  = 0;
    p1_WFPS_acc   = 0;
    p1_fanaer_acc = 0;
    p1_peak_LAI   = zeros(1,0);
    p1_Lk_acc     = 0;
    p1_T_H_acc      = zeros(1, cc);
    p1_T_L_acc      = zeros(1, cc);
    p1_LAI_H_acc    = zeros(1, cc);
    p1_LAI_L_acc    = zeros(1, cc);
    p1_B_H_root_acc = zeros(1, cc);
    p1_B_L_root_acc = zeros(1, cc);
    p1_SupN_H_acc   = zeros(1, cc);
    p1_SupP_H_acc   = zeros(1, cc);
    p1_SupK_H_acc   = zeros(1, cc);
    p1_SupCa_H_acc  = zeros(1, cc);
    p1_SupMg_H_acc  = zeros(1, cc);
    p1_SupSi_H_acc  = zeros(1, cc);
    p1_SupN_L_acc   = zeros(1, cc);
    p1_SupP_L_acc   = zeros(1, cc);
    p1_SupK_L_acc   = zeros(1, cc);
    p1_SupCa_L_acc  = zeros(1, cc);
    p1_SupMg_L_acc  = zeros(1, cc);
    p1_SupSi_L_acc  = zeros(1, cc);
    p1_IS_acc       = [];
    p1_REXMY_acc    = [];

    p1_keep = [who; { ...
        'p1_keep'; 'p1_yi'; 'p1_n'; 'p1_chk'; 'p2_inprog_chk'; 'p2_chk'; 'p3_chk'; ...
        'p1_Ta_acc'; 'p1_Tdp_acc'; 'p1_Se_acc'; 'p1_Psi_acc'; ...
        'p1_VSUM_acc'; 'p1_VTSUM_acc'; 'p1_WFPS_acc'; 'p1_fanaer_acc'; 'p1_peak_LAI'; 'p1_Lk_acc'; ...
        'p1_T_H_acc'; 'p1_T_L_acc'; 'p1_LAI_H_acc'; 'p1_LAI_L_acc'; ...
        'p1_B_H_root_acc'; 'p1_B_L_root_acc'; ...
        'p1_SupN_H_acc'; 'p1_SupP_H_acc'; 'p1_SupK_H_acc'; ...
        'p1_SupCa_H_acc'; 'p1_SupMg_H_acc'; 'p1_SupSi_H_acc'; ...
        'p1_SupN_L_acc'; 'p1_SupP_L_acc'; 'p1_SupK_L_acc'; ...
        'p1_SupCa_L_acc'; 'p1_SupMg_L_acc'; 'p1_SupSi_L_acc'; ...
        'p1_IS_acc'; 'p1_REXMY_acc'; ...
        'OPT_SoilBiogeochemistry_saved'; 'OPT_pH_saved'; 'OPT_ERW_saved'; ...
        'id_location'; 't_bef'; 't_aft'; 'DeltaGMT'; 'Lon'; 'Lat'; 'O_out'}];

    OPT_SoilBiogeochemistry_saved = OPT_SoilBiogeochemistry;
    OPT_pH_saved                  = OPT_pH;
    OPT_ERW_saved                 = OPT_ERW;
    OPT_SoilBiogeochemistry = 0;
    OPT_pH                  = 0;
    %% ERW chemistry is seeded from the biogeochemical pools (B_out) and
    %% cannot run in this T&C-only phase; it returns with them in Phase 2.
    OPT_ERW                 = 0;

    for p1_yi = 1:n_spinup_years
        if p1_yi > 1
            clearvars('-except', p1_keep{:});
        end
        cycle_year = spinup_cycle_years(p1_yi);
        fprintf('  Phase 1: year %d/%d (meteo %d)\n', p1_yi, n_spinup_years, cycle_year);

        yr_mask = (spinup_all_yr == cycle_year);
        Date = spinup_all_met.Date(yr_mask);
        Pr   = spinup_all_met.Pr(yr_mask);
        Ta   = spinup_all_met.Ta(yr_mask);
        Ws   = spinup_all_met.Ws(yr_mask);
        ea   = spinup_all_met.ea(yr_mask);
        Pre  = spinup_all_met.Pre(yr_mask);
        SAD1 = spinup_all_met.SAD1(yr_mask);
        SAD2 = spinup_all_met.SAD2(yr_mask);
        SAB1 = spinup_all_met.SAB1(yr_mask);
        SAB2 = spinup_all_met.SAB2(yr_mask);
        N    = spinup_all_met.N(yr_mask);
        Tdew = spinup_all_met.Tdew(yr_mask);
        esat = spinup_all_met.esat(yr_mask);
        PARB = spinup_all_met.PARB(yr_mask);
        PARD = spinup_all_met.PARD(yr_mask);
        NN   = numel(Date);
        Ws(Ws <= 0) = 0.01;
        Ds   = max(esat - ea, 0);
        [YE, MO, DA, HO, ~, ~] = datevec(Date);
        Datam = [YE, MO, DA, HO];
        Ca = Load_CO2_Record(Date, Input_Directory, co2_record_file);
        clear YE MO DA HO

        PARAM_IC     = fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']);
        ERW_PARAM_IC = fullfile(Parameter_Directory, 'ERW_MOD_PARAM_excel.m');

        %% Series the driver creates only when absent (irrigation, wetland
        %% level): one left from another year has that year's length.
        clear IrD Wlev
        cd(Directory);
        MAIN_FRAME_BG_ERW_10mins;
        cd(current_directory);

        % Extract averages; skip sentinel row 1
        NNd_y = ceil(NN/24) + 1;
        jrows = 2:NNd_y;

        % Carry ending soil moisture into next cycle year (MOD_PARAM_*.m
        % restores O(1,:) from O_out independently of OPT_SoilBiogeochemistry)
        % so Phase 1 is a continuous multi-year run rather than resetting to
        % field capacity every year.
        O_out = O(NNd_y, :);
        %% A perennial stand is not carried through Phase 1: without the
        %% biogeochemistry the plants take up no nutrients, yet growth still
        %% draws on their reserves, so a carried stand would starve (Hamden's
        %% grass within a year, its trees in four).  Each Phase-1 year starts
        %% from the vegetation file, as in v0.9; the stand is carried from
        %% Phase 3 on.

        p1_n       = p1_n + 1;
        p1_Ta_acc  = p1_Ta_acc  + mean(Ta,                    'omitnan');
        p1_Tdp_acc = p1_Tdp_acc + mean(Tdp_bio_out(jrows),   'omitnan');
        p1_Se_acc  = p1_Se_acc  + mean(Se_bio_out(jrows),    'omitnan');
        p1_Psi_acc = p1_Psi_acc + mean(Psi_bio_out(jrows),   'omitnan');
        %%% v1.0, leaching_water 'layer' (the frame's leach_layer): Phase 2
        %%% leaches with the layer's own drainage over the layer's water, so
        %%% these two carry the frame's daily records of those instead.
        if leach_layer
            p1_VSUM_acc = p1_VSUM_acc + mean(V_leach_out(jrows),'omitnan');
        else
            p1_VSUM_acc  = p1_VSUM_acc  + mean(VSUM_out(jrows),  'omitnan');
        end
        p1_VTSUM_acc = p1_VTSUM_acc + mean(VTSUM_out(jrows),  'omitnan');
        p1_WFPS_acc  = p1_WFPS_acc  + mean(WFPS_bio_out(jrows),'omitnan');
        p1_fanaer_acc = p1_fanaer_acc + mean(f_anaer_bio_out(jrows),'omitnan');
        %% Peak leaf area of the year over every crown of both layers.
        p1_peak_LAI(end+1) = max([reshape(LAI_H(jrows,:),1,[]), ...
            reshape(LAI_L(jrows,:),1,[]), 0]);
        if leach_layer
            p1_Lk_acc = p1_Lk_acc + mean(Lk_leach_out(jrows),'omitnan');
        else
            p1_Lk_acc    = p1_Lk_acc    + sum(Lk(1:NN)) / (NNd_y - 1);
        end

        p1_T_H_acc      = p1_T_H_acc      + sum(T_H(1:NN,:), 1) / (NNd_y - 1);
        p1_T_L_acc      = p1_T_L_acc      + sum(T_L(1:NN,:), 1) / (NNd_y - 1);
        p1_LAI_H_acc    = p1_LAI_H_acc    + mean(LAI_H(jrows,:),    1, 'omitnan');
        p1_LAI_L_acc    = p1_LAI_L_acc    + mean(LAI_L(jrows,:),    1, 'omitnan');
        p1_B_H_root_acc = p1_B_H_root_acc + mean(B_H(jrows,:,3),   1, 'omitnan');
        p1_B_L_root_acc = p1_B_L_root_acc + mean(B_L(jrows,:,3),   1, 'omitnan');
        p1_SupN_H_acc   = p1_SupN_H_acc   + mean(SupN_H(jrows,:),  1, 'omitnan');
        p1_SupP_H_acc   = p1_SupP_H_acc   + mean(SupP_H(jrows,:),  1, 'omitnan');
        p1_SupK_H_acc   = p1_SupK_H_acc   + mean(SupK_H(jrows,:),  1, 'omitnan');
        p1_SupCa_H_acc  = p1_SupCa_H_acc  + mean(SupCa_H(jrows,:), 1, 'omitnan');
        p1_SupMg_H_acc  = p1_SupMg_H_acc  + mean(SupMg_H(jrows,:), 1, 'omitnan');
        p1_SupSi_H_acc  = p1_SupSi_H_acc  + mean(SupSi_H(jrows,:), 1, 'omitnan');
        p1_SupN_L_acc   = p1_SupN_L_acc   + mean(SupN_L(jrows,:),  1, 'omitnan');
        p1_SupP_L_acc   = p1_SupP_L_acc   + mean(SupP_L(jrows,:),  1, 'omitnan');
        p1_SupK_L_acc   = p1_SupK_L_acc   + mean(SupK_L(jrows,:),  1, 'omitnan');
        p1_SupCa_L_acc  = p1_SupCa_L_acc  + mean(SupCa_L(jrows,:), 1, 'omitnan');
        p1_SupMg_L_acc  = p1_SupMg_L_acc  + mean(SupMg_L(jrows,:), 1, 'omitnan');
        p1_SupSi_L_acc  = p1_SupSi_L_acc  + mean(SupSi_L(jrows,:), 1, 'omitnan');

        % IS: Ccrown-weighted sum of H and L radiation interception (1x27)
        ISOIL_H_avg = reshape(squeeze(mean(ISOIL_H(jrows,:,:), 1)), cc, []);
        ISOIL_L_avg = reshape(squeeze(mean(ISOIL_L(jrows,:,:), 1)), cc, []);
        IS_y = reshape(Ccrown, 1, []) * ISOIL_H_avg + reshape(Ccrown, 1, []) * ISOIL_L_avg;
        % REXMY (1x3); zeros when OPT_SoilBiogeochemistry=0 (Rexmy arrays stay 0)
        Rexmy_H_avg = reshape(squeeze(mean(Rexmy_H(jrows,:,:), 1)), cc, []);
        Rexmy_L_avg = reshape(squeeze(mean(Rexmy_L(jrows,:,:), 1)), cc, []);
        REXMY_y = reshape(Ccrown, 1, []) * Rexmy_H_avg + reshape(Ccrown, 1, []) * Rexmy_L_avg;

        if isempty(p1_IS_acc)
            p1_IS_acc    = IS_y;
            p1_REXMY_acc = REXMY_y;
        else
            p1_IS_acc    = p1_IS_acc    + IS_y;
            p1_REXMY_acc = p1_REXMY_acc + REXMY_y;
        end
    end

    OPT_SoilBiogeochemistry = OPT_SoilBiogeochemistry_saved;
    OPT_pH                  = OPT_pH_saved;
    OPT_ERW                 = OPT_ERW_saved;

    % Fail fast rather than propagate NaN/Inf averages into Phase 2 if Phase 1
    % accumulated no timesteps (e.g. an empty or misconfigured spin-up window).
    if ~(isscalar(p1_n) && isfinite(p1_n) && p1_n > 0)
        error('Run_BG_Spinup:NoPhase1Data', ...
            ['Phase 1 accumulated no timesteps (p1_n = %s); check the spin-up ' ...
            'meteo window (meteo_start/meteo_end) and cycle years.'], mat2str(p1_n));
    end
    p1_O = reshape(O_out, 1, []);
    p1_Ohy = reshape(Ohy, 1, []);
    p1_Osat = reshape(Osat, 1, []);
    if mod(numel(p1_O), numel(p1_Ohy)) ~= 0 || ...
            mod(numel(p1_O), numel(p1_Osat)) ~= 0
        error('Run_BG_Spinup:InvalidPhase1State', ...
            'Phase 1 moisture bounds cannot be mapped onto the terminal state.');
    end
    p1_Ohy = repmat(p1_Ohy, 1, numel(p1_O)/numel(p1_Ohy));
    p1_Osat = repmat(p1_Osat, 1, numel(p1_O)/numel(p1_Osat));
    p1_water_state_physical = all(isfinite(p1_O)) && ...
        all(p1_O >= p1_Ohy - 1e-8) && all(p1_O <= p1_Osat + 1e-8);
    if p1_n ~= n_spinup_years || ~p1_water_state_physical
        error('Run_BG_Spinup:InvalidPhase1State', ...
            'Phase 1 must complete its forcing ring with physically bounded soil moisture.');
    end
    p1_bare = p1_peak_LAI < P1_MIN_PEAK_LAI;
    if any(Ccrown > 0) && any(p1_bare)
        error('Run_BG_Spinup:VegetationNotEstablished', ...
            ['Phase 1 cycle year(s) %s reach a peak leaf area of only %s m2 m-2 ' ...
             '(at least %.2g required): the vegetation did not establish.'], ...
            mat2str(spinup_cycle_years(p1_bare)), mat2str(p1_peak_LAI(p1_bare),3), ...
            P1_MIN_PEAK_LAI);
    end
    clear p1_bare

    p1_avgs.Ta       = p1_Ta_acc   / p1_n;
    p1_avgs.Tdp_bio  = p1_Tdp_acc  / p1_n;
    p1_avgs.Se_bio   = p1_Se_acc   / p1_n;
    p1_avgs.Psi_bio  = p1_Psi_acc  / p1_n;
    p1_avgs.VSUM     = p1_VSUM_acc  / p1_n;
    p1_avgs.VTSUM    = p1_VTSUM_acc / p1_n;
    p1_avgs.WFPS_bio = p1_WFPS_acc  / p1_n;
    %% Mean anaerobic volume fraction of the active layer, the moisture-
    %% dependent input BIOGEO_UNIT's denitrification reads (Biogeo_environment).
    p1_avgs.f_anaer_bio = p1_fanaer_acc / p1_n;
    p1_avgs.Lk_day   = p1_Lk_acc    / p1_n;
    p1_avgs.T_H_day      = p1_T_H_acc      / p1_n;
    p1_avgs.T_L_day      = p1_T_L_acc      / p1_n;
    p1_avgs.LAI_H        = p1_LAI_H_acc    / p1_n;
    p1_avgs.LAI_L        = p1_LAI_L_acc    / p1_n;
    p1_avgs.B_H_root     = p1_B_H_root_acc / p1_n;
    p1_avgs.B_L_root     = p1_B_L_root_acc / p1_n;
    p1_avgs.SupN_H       = p1_SupN_H_acc   / p1_n;
    p1_avgs.SupP_H       = p1_SupP_H_acc   / p1_n;
    p1_avgs.SupK_H       = p1_SupK_H_acc   / p1_n;
    p1_avgs.SupCa_H      = p1_SupCa_H_acc  / p1_n;
    p1_avgs.SupMg_H      = p1_SupMg_H_acc  / p1_n;
    p1_avgs.SupSi_H      = p1_SupSi_H_acc  / p1_n;
    p1_avgs.SupN_L       = p1_SupN_L_acc   / p1_n;
    p1_avgs.SupP_L       = p1_SupP_L_acc   / p1_n;
    p1_avgs.SupK_L       = p1_SupK_L_acc   / p1_n;
    p1_avgs.SupCa_L      = p1_SupCa_L_acc  / p1_n;
    p1_avgs.SupMg_L      = p1_SupMg_L_acc  / p1_n;
    p1_avgs.SupSi_L      = p1_SupSi_L_acc  / p1_n;
    p1_avgs.IS           = p1_IS_acc    / p1_n;
    p1_avgs.REXMY        = p1_REXMY_acc / p1_n;
    p1_avgs.Se_fc        = Se_fc;  % constant soil parameter, preserved from MAIN_FRAME
    p1_avgs.O_out        = O_out;
    p1_avgs.phase1_complete = p1_n == n_spinup_years;
    p1_avgs.phase1_water_state_physical = p1_water_state_physical;

    p1_avgs.spinup_cycle_years = spinup_cycle_years;
    p1_avgs.standing_nutrient_state_vers = STANDING_NUTRIENT_STATE_VERS;
    p1_avgs.spinup_provenance = spinup_provenance;
    p1_avgs.peak_LAI_by_year = p1_peak_LAI;
    save(p1_chk, 'p1_avgs', '-v7.3');
    fprintf('Phase 1 complete. Averages saved to: %s\n', p1_chk);
end

%% Load p1_avgs if Phase 2 still needs them
need_p1_avgs = ~spinup_handoff_only && ~p2_done && ~p3_done;
if need_p1_avgs
    if ~exist('p1_avgs', 'var')
        tmp = load(p1_chk, 'p1_avgs');
        p1_avgs = tmp.p1_avgs;
        clear tmp
    end
    O_out = p1_avgs.O_out;
end

%% =========================================================================
%% Phase 2: Standalone BIOGEO_UNIT; SOM steady state
%% =========================================================================
if ~spinup_handoff_only && ~p2_done && ~p3_done
    fprintf('\n=== Phase 2: Standalone BIOGEO_UNIT (until convergence) ===\n');

    % Establish PARAM_IC variables (Zbio, Bio_Zs, ExEM, B_IO, Ccrown, etc.).
    % PARAM_IC references NN and NNd (MAIN_FRAME-internal) for Lmax_day and
    % per-day pre-allocations.  Phase 2 never calls VEGGIE_UNIT so Lmax_day is
    % irrelevant here; NN=0 makes the daylength loop a no-op (no Datam needed).
    % Phase 3 MAIN_FRAME re-runs PARAM_IC with real NN/Datam, recomputing Lmax_day.
    % Clear per-PFT struct arrays first: indexed assignment in PARAM_IC fails if
    % they already exist with different fields from Phase 1.
    NN  = 0;
    NNd = 1;
    Ta  = p1_avgs.Ta;
    O   = zeros(1, ms);
    P   = zeros(1, 79);
    clear Mpar_L Mpar_H ParEx_L ParEx_H VegH_Param_Dyn VegL_Param_Dyn
    Ca = 400;  % placeholder; Phase 2 does not use Ci/Ca (no VEGGIE_UNIT)
    cd(Directory);
    run(fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']));
    %%% The row's vegetation override, as MAIN_FRAME_BG_ERW_10mins applies it in
    %%% Phases 1 and 3: Phase 2 then equilibrates with the site's own roots,
    %%% deposition and management, not the borrowed row's.  A row without one
    %%% runs as before.
    if istable(MOD_PARAM_input) && ...
            ismember('pft_override',MOD_PARAM_input.Properties.VariableNames) && ...
            strlength(strtrim(string(MOD_PARAM_input.pft_override))) > 0
        run(fullfile(Parameter_Directory, ...
            [char(strtrim(string(MOD_PARAM_input.pft_override))) '.m']));
    end
    cd(current_directory);
    %%% The legacy standalone spin-up runs at unit crown (Biogeo_Spin_Up.m sets
    %%% Ccrown=1) so prescribed litter return and plant uptake are both active.
    %%% The scenario workbook leaves Ccrown at zero because the coupled phases
    %%% set it per management interval in CropHeightType, which Phase 2 never
    %%% calls; left as-is that silences uptake while litter keeps arriving from
    %%% the Phase 1 averages, making every nutrient a one-way source.  Unit
    %%% weight per tile: the Phase 1 per-tile means (root, LAI, Sup*) are raw
    %%% and already carry each tile's rotation duty cycle through its off-year
    %%% zeros, so weighting them again would apply that cycle twice.
    Ccrown = ones(1, cc);
    % Phase 2 never calls VEGGIE_UNIT, so there is no crop uptake to balance
    % a recurring fertiliser input: smearing the annual total across all 366
    % days of every cycled year compounds into an
    % unbounded nutrient pool in a prior site run.
    % Phase 2 represents background/native biogeochemical equilibration;
    % fertiliser management is applied in the actual transient run instead.
    B_IO_p2 = B_IO;
    B_IO_p2.FertN  = zeros(1,366);
    B_IO_p2.FertP  = zeros(1,366);
    B_IO_p2.FertK  = zeros(1,366);
    B_IO_p2.FertCa = zeros(1,366);
    B_IO_p2.FertMg = zeros(1,366);
    B_IO_p2.FertSi = zeros(1,366);
    if isfield(B_IO_p2,'FertAn'); B_IO_p2.FertAn = zeros(1,366); end
    if isfield(B_IO_p2,'f_fert_NH4'); B_IO_p2.f_fert_NH4 = 0.75 * ones(1,366); end
    if isfield(B_IO_p2,'fert_frac_NO3'); B_IO_p2.fert_frac_NO3 = zeros(1,366); end
    if isfield(B_IO_p2,'fert_frac_NH4'); B_IO_p2.fert_frac_NH4 = zeros(1,366); end
    if isfield(B_IO_p2,'fert_frac_urea'); B_IO_p2.fert_frac_urea = zeros(1,366); end
    %%% Event deposition (B_IO.Dep_events, Deposition_Day): Phase 2 runs at the
    %%% ring's constant mean forcing, so it takes the ring's mean daily
    %%% deposition -- each ring year's dry deposition over its days and its wet
    %%% deposition on its precipitation, weighted as Phase 1's averages are.
    if isfield(B_IO_p2,'Dep_events') && ~isempty(B_IO_p2.Dep_events)
        p2_dep_total = zeros(1,6);
        p2_dep_days = 0;
        for p2_dep_year = reshape(spinup_cycle_years,1,[])
            p2_dep_mask = (spinup_all_yr == p2_dep_year);
            p2_dep_n = nnz(p2_dep_mask)*dth/24;
            [~,p2_dep_rate] = Deposition_Day(B_IO_p2,[p2_dep_year 1 1 0 0 0], ...
                sum(spinup_all_met.Pr(p2_dep_mask))*dth,p2_dep_n);
            p2_dep_total = p2_dep_total+p2_dep_n*p2_dep_rate;
            p2_dep_days = p2_dep_days+p2_dep_n;
        end
        p2_dep_rate = p2_dep_total/p2_dep_days;
        B_IO_p2.DepN = p2_dep_rate(1);  B_IO_p2.DepK = p2_dep_rate(3);
        B_IO_p2.DepCa = p2_dep_rate(4); B_IO_p2.DepMg = p2_dep_rate(5);
        B_IO_p2 = rmfield(B_IO_p2,'Dep_events');
        fprintf(['  Phase 2 deposition, the ring''s daily mean [g/m2/d]: N %.4g, ' ...
            'K %.4g, Ca %.4g, Mg %.4g\n'],p2_dep_rate([1 3 4 5]));
    end

    if any(logical(Ccrown))
        OPT_BfixN_p2 = max(ParEx_L(logical(Ccrown)).bfix);
    else
        OPT_BfixN_p2 = 0;
    end

    % Unpack constant forcing from p1_avgs
    p2_Ta        = p1_avgs.Ta;
    p2_Tdp_bio   = p1_avgs.Tdp_bio;
    p2_Se_bio    = p1_avgs.Se_bio;
    p2_Psi_bio   = p1_avgs.Psi_bio;
    p2_VSUM      = p1_avgs.VSUM;
    p2_VTSUM     = p1_avgs.VTSUM;
    p2_WFPS_bio  = p1_avgs.WFPS_bio;
    p2_Lk_day    = p1_avgs.Lk_day;
    p2_T_H_day   = p1_avgs.T_H_day;
    p2_T_L_day   = p1_avgs.T_L_day;
    p2_LAI_H     = p1_avgs.LAI_H;
    p2_LAI_L     = p1_avgs.LAI_L;
    p2_B_H_root  = p1_avgs.B_H_root;
    p2_B_L_root  = p1_avgs.B_L_root;
    p2_SupN_H    = p1_avgs.SupN_H;
    p2_SupP_H    = p1_avgs.SupP_H;
    p2_SupK_H    = p1_avgs.SupK_H;
    p2_SupCa_H   = p1_avgs.SupCa_H;
    p2_SupMg_H   = p1_avgs.SupMg_H;
    p2_SupSi_H   = p1_avgs.SupSi_H;
    p2_SupN_L    = p1_avgs.SupN_L;
    p2_SupP_L    = p1_avgs.SupP_L;
    p2_SupK_L    = p1_avgs.SupK_L;
    p2_SupCa_L   = p1_avgs.SupCa_L;
    p2_SupMg_L   = p1_avgs.SupMg_L;
    p2_SupSi_L   = p1_avgs.SupSi_L;
    p2_IS        = p1_avgs.IS;
    p2_REXMY     = p1_avgs.REXMY;
    p2_Se_fc     = p1_avgs.Se_fc;
    p2_f_anaer_bio = p1_avgs.f_anaer_bio;

    % Initial BG state: resume in-progress or cold start
    p2ip_compatible = false;
    if isfile(p2_inprog_chk)
        tmp = load(p2_inprog_chk, 'p2ip');
        p2ip_compatible = isfield(tmp, 'p2ip') && ...
            isfield(tmp.p2ip, 'spinup_cycle_years') && ...
            isequal(tmp.p2ip.spinup_cycle_years, spinup_cycle_years) && ...
            isfield(tmp.p2ip, 'observation_constraint') && ...
            isequaln(tmp.p2ip.observation_constraint, BG_SPINUP_CONSTRAINT) && ...
            isfield(tmp.p2ip, 'standing_nutrient_state_vers') && ...
            tmp.p2ip.standing_nutrient_state_vers == STANDING_NUTRIENT_STATE_VERS && ...
            isfield(tmp.p2ip, 'P_prev') && ...
            spinup_provenance_matches(tmp.p2ip, spinup_provenance);
        if ~p2ip_compatible
            fprintf('  Ignoring incompatible Phase 2 in-progress checkpoint.\n');
        end
    end
    p2_constraint_cumulative_delta = zeros(1,79);
    p2_constraint_diag = struct('applied', false, 'delta_by_pool', zeros(1,79));
    if p2ip_compatible
        p2_P_now     = tmp.p2ip.P_now;
        p2_RexmyI    = tmp.p2ip.RexmyI;
        p2_NavlI     = tmp.p2ip.NavlI;
        p2_yr_start  = tmp.p2ip.p2_year + 1;
        p2_streak    = tmp.p2ip.streak;
        p2_P_prev    = tmp.p2ip.P_prev;
        if isfield(tmp.p2ip, 'constraint_cumulative_delta')
            p2_constraint_cumulative_delta = tmp.p2ip.constraint_cumulative_delta;
        end
        if isfield(tmp.p2ip, 'constraint_diag')
            p2_constraint_diag = tmp.p2ip.constraint_diag;
        end
        fprintf('  Resumed Phase 2 from year %d\n', p2_yr_start);
    else
        p2_P_now    = P(1,:);
        p2_RexmyI   = RexmyI(1,:);
        p2_NavlI    = NavlI(1,:);
        p2_yr_start = 1;
        p2_streak   = 0;
        p2_P_prev   = p2_P_now;
    end
    clear tmp
    p2_P_now = Project_Small_Negative_BG_Pools(p2_P_now, P2_POOL_NEG_TOL);
    %%% The measured slow stocks enter once, here, as the biogeochemical
    %%% initial condition.  Re-imposing them every year or every ring makes
    %%% them a forcing instead: the drift a convergence gate then measures is
    %%% the size of its own reset, and a pool the constraint hard-sets reads as
    %%% exactly zero drift and can never fail.  A resumed Phase 2 carries a
    %%% state the constraint has already shaped, so it is not applied again.
    if ~p2ip_compatible
        [p2_P_now, p2_constraint_diag] = spinup_apply_observation_constraint( ...
            p2_P_now, BG_SPINUP_CONSTRAINT);
        p2_constraint_cumulative_delta = p2_constraint_diag.delta_by_pool;
        p2_P_prev = p2_P_now;
    end
    clear p2ip_compatible
    p2_P_prev = Project_Small_Negative_BG_Pools(p2_P_prev, P2_POOL_NEG_TOL);

    %%% Phase 2 runs the chemistry-off biogeochemical spin-up at the observed
    %%% bulk pH.  CEC is prescribed or estimated before the Phase 3 handoff, on the
    %%% settled pools rather than this fertilizer-pulse initial state.
    PH = PHs;

    %%% A tracked cell's lateral exchange as one mean daily map over its ring
    %%% (the retained share and the arriving mass of each mobile pool), the
    %%% lateral counterpart of the Phase-1 mean forcing above.
    if lateral_replay_spinup
        [p2_lat_yr,~,~] = datevec(LATERAL_REPLAY.day_dn);
        p2_lat_days = ismember(p2_lat_yr,spinup_cycle_years);
        if ~any(p2_lat_days)
            error('Run_BG_Spinup:LateralRecordNotCovered', ...
                'The lateral record holds no day of the cycle years.');
        end
        p2_lat_idx = LATERAL_REPLAY.pool_idx(:).';
        p2_lat_retain = mean(LATERAL_REPLAY.retain(p2_lat_days,:),1);
        p2_lat_input = mean(LATERAL_REPLAY.input(p2_lat_days,:),1);
        p2_lat_isP = ismember(p2_lat_idx,35:47);
        clear p2_lat_yr p2_lat_days
    end

    p2_converged = false;
    p2_yr_final = p2_yr_start - 1;
    p2_R_het_last_year = NaN;

    p2_yr = p2_yr_start - 1;
    while ~p2_converged
        p2_yr = p2_yr + 1;
        p2_R_het_year_sum = 0;
        p2_R_het_year_n = 0;
        %%% Phase 3 closes a phosphorus budget every loop; Phase 2 never did, so
        %%% a standalone-path source had nothing to trip it.  Same terms and the
        %%% same closure test, accumulated over the mean-forcing year.
        p2_p_terms = struct('fertilizer',0,'deposition',0,'uplift',0, ...
            'planting',0,'litter',0,'harvest',0,'leaching',0, ...
            'soil_start',sum(p2_P_now(35:47)),'soil_end',0, ...
            'plant_start',0,'plant_end',0);
        if lateral_replay_spinup
            p2_p_terms.lateral_in = 0;
            p2_p_terms.lateral_out = 0;
        end
        FertN_avail_state = [];
        for p2_day = 1:365
            dvec = datevec(datenum(2001, 1, p2_day));
            Datam_day = dvec(1:4);
            [B_IO_day, FertN_avail_state] = Fertilizer_N_Availability_Step( ...
                B_IO_p2, Datam_day, FertN_avail_state, 1, p2_Tdp_bio);

            [p2_P_now, ~,~, p2_LEAK_P, ~,~,~, p2_LEAK_DOP, ~,~,~, ...  % 1-11
             ~,~,~, p2_Puptake_H, ~,~,~,~, p2_Puptake_L, ~, ...        % 12-21
             ~,~,~,~,~,~, ...                       % 22-27
             p2_RexmyI, ...                         % 28
             p2_R_litter, p2_R_microbe, p2_R_litter_sur, p2_R_ew, ~,~,~,~, ... % 29-36
             ~,~,~,~,~,~, p2_NavlI, ~, ...         % 37-43, 44
             ~,~,~,~,~,~, p2_P_inputs] = ...        % 45-51 (v1.0: N2Oflx FertN DepN NITflx DENITflx N2Oflx_nit P_inputs)
                BIOGEO_UNIT(p2_P_now, p2_IS, Zbio, sum(Bio_Zs.*rsd), ...
                PH, p2_Tdp_bio, p2_Ta, p2_Psi_bio, p2_Se_bio, p2_Se_fc, p2_VSUM, p2_VTSUM, ...
                Ccrown, Bio_Zs, RfH_Zs, RfL_Zs, ...
                p2_Lk_day, p2_T_H_day, p2_T_L_day, p2_B_H_root, p2_B_L_root, p2_LAI_H, p2_LAI_L, ...
                p2_SupN_H, p2_SupP_H, p2_SupK_H, p2_SupN_L, p2_SupP_L, p2_SupK_L, ...
                p2_SupCa_H, p2_SupMg_H, p2_SupSi_H, p2_SupCa_L, p2_SupMg_L, p2_SupSi_L, ...
                p2_REXMY, p2_RexmyI, ExEM, p2_NavlI, Pcla, Psan, ...
                B_IO_day, Datam_day, 0, 0, ...
                0, OPT_BfixN_p2, ...              % opt_ERW (off in Phase 2), opt_BfixN,
                p2_f_anaer_bio, p2_WFPS_bio, ...  % mean anaerobic fraction and WFPS of the active layer,
                [], [], [], leach_layer);         % exchange scheme (inert with ERW off), default KNmaxd as in Phase 3, no overrides, the leaching water
            if lateral_replay_spinup
                p2_lat_T = p2_P_now(p2_lat_idx);
                p2_P_now(p2_lat_idx) = p2_lat_retain.*p2_lat_T + p2_lat_input;
                p2_p_terms.lateral_in = p2_p_terms.lateral_in + ...
                    sum(p2_lat_input(p2_lat_isP));
                p2_p_terms.lateral_out = p2_p_terms.lateral_out + ...
                    sum((1 - p2_lat_retain(p2_lat_isP)).*p2_lat_T(p2_lat_isP));
            end
            p2_R_het_day = p2_R_microbe + p2_R_litter + p2_R_ew - p2_R_litter_sur;
            if isfinite(p2_R_het_day)
                p2_R_het_year_sum = p2_R_het_year_sum + p2_R_het_day;
                p2_R_het_year_n = p2_R_het_year_n + 1;
            end
            %%% Phase 2 carries no plant pools, so root uptake leaves the soil
            %%% inventory for good and belongs where harvest sits in Phase 3.
            p2_p_terms.fertilizer = p2_p_terms.fertilizer + p2_P_inputs(1);
            p2_p_terms.deposition = p2_p_terms.deposition + p2_P_inputs(2);
            p2_p_terms.uplift     = p2_p_terms.uplift     + p2_P_inputs(3);
            p2_p_terms.leaching   = p2_p_terms.leaching   + p2_LEAK_P + p2_LEAK_DOP;
            %%% ISOIL P components 13:15 (above-ground, woody, below-ground),
            %%% already Ccrown-weighted to a per-m2-ground flux in Phase 1.
            p2_p_terms.litter     = p2_p_terms.litter + sum(p2_IS(13:15));
            %%% BIOGEO_UNIT reports uptake per unit crown (./Ccrown), so the
            %%% areal flux multiplies it back.  Phase 2 never calls
            %%% CropHeightType, so Ccrown can be zero here; a zero-crown tile
            %%% takes up nothing and its non-finite per-crown rate is dropped
            %%% rather than propagated into the budget.
            p2_uptake_day = Ccrown(:).'.*(p2_Puptake_H(:).' + p2_Puptake_L(:).');
            p2_p_terms.harvest = p2_p_terms.harvest + ...
                sum(p2_uptake_day(isfinite(p2_uptake_day)));
        end

        %%% Taken before the negative-pool projection so the residual measures
        %%% the dynamics alone rather than the projection's own rounding.
        p2_p_terms.soil_end = sum(p2_P_now(35:47));
        p2_p_budget = Phase3_Phosphorus_Budget( ...
            p2_p_terms, P3_NUTRIENT_BUDGET_ABS_TOL, P3_NUTRIENT_BUDGET_REL_TOL);
        p2_p_budget.terms = p2_p_terms;

        p2_P_now = Project_Small_Negative_BG_Pools(p2_P_now, P2_POOL_NEG_TOL);

        if p2_R_het_year_n > 0
            p2_R_het_last_year = p2_R_het_year_sum / p2_R_het_year_n;
        end

        if p2_yr > 1
            %%% Constant mean forcing, so there is no ring: each pool's own
            %%% level is its scale here.
            p2_diag = spinup_pool_group_diagnostics(p2_P_now, p2_P_prev, 1, ...
                P2_POOL_SCALE_FLOOR, abs(p2_P_prev));
            p2_pool_pass = p2_diag.pool_ring_rel < P2_ANNUAL_REL_TOL;
            p2_dynamic_pass = p2_diag.dynamic_ring_rel < P2_ANNUAL_REL_TOL;
            p2_slow_pass = p2_diag.slow_ring_rel < P2_ANNUAL_REL_TOL;
            p2_max_rel = max([p2_diag.dynamic_ring_rel p2_diag.slow_ring_rel]);
            %%% Phase 2 measures drift across one mean-forcing year.  Reporting
            %%% it over the ring puts it on Phase 3's scale, so both phases'
            %%% output reads against the same ring limit.  The pass test is
            %%% unchanged: the annual value against the annual limit is the
            %%% same comparison.
            p2_dynamic_ring_rel = p2_diag.dynamic_ring_rel*n_spinup_years;
            p2_slow_ring_rel = p2_diag.slow_ring_rel*n_spinup_years;
            p2_pass = p2_diag.finite_nonnegative && all(p2_pool_pass) && ...
                all(p2_dynamic_pass) && all(p2_slow_pass);
            [p2_streak, p2_converged] = spinup_update_consecutive_passes( ...
                p2_pass, p2_streak, P2_STREAK_REQ);
            if mod(p2_yr, 50) == 0 || p2_converged
                fprintf(['  Phase 2 yr %d: dynamic ring max=%.4f slow ring max=%.4f ' ...
                    'annual max=%.4g yr^-1 streak=%d/%d\n'], p2_yr, ...
                    max(p2_dynamic_ring_rel), max(p2_slow_ring_rel), p2_max_rel, ...
                    p2_streak, P2_STREAK_REQ);
            end
        end
        p2_P_prev  = p2_P_now;
        p2_yr_final = p2_yr;

        % In-progress checkpoint every 50 years
        if mod(p2_yr, 50) == 0 && ~p2_converged
            p2ip = struct('P_now', p2_P_now, 'RexmyI', p2_RexmyI, 'NavlI', p2_NavlI, ...
                'p2_year', p2_yr, 'streak', p2_streak, 'P_prev', p2_P_prev, ...
                'spinup_cycle_years', spinup_cycle_years, ...
                'standing_nutrient_state_vers', STANDING_NUTRIENT_STATE_VERS, ...
                'spinup_provenance', spinup_provenance, ...
                'observation_constraint', BG_SPINUP_CONSTRAINT, ...
                'constraint_diag', p2_constraint_diag, ...
                'constraint_cumulative_delta', p2_constraint_cumulative_delta);
            save(p2_inprog_chk, 'p2ip', '-v7.3');
        end

        if p2_converged
            fprintf('Phase 2 converged at year %d\n', p2_yr);
            break
        end
    end

    % Steady-state annual mean heterotrophic decomposition flux, assembled the
    % same way MAIN_FRAME does (R_microbe+R_litter+R_ew-R_litter_sur, converted
    % gC/m2/d -> umol C/m2/d).  Phase 3 captures its own value from the coupled
    % loop; this is the fallback the transient's soil-air CO2 boundary falls
    % back to.  There is no r_aut counterpart because Phase 2 has no live
    % vegetation (no VEGGIE_UNIT root respiration).
    if ~isfinite(p2_R_het_last_year)
        p2_R_het_last_year = p2_R_microbe + p2_R_litter + p2_R_ew - p2_R_litter_sur;
    end
    p2_r_het_ss = p2_R_het_last_year / 12 * 1e6;

    % Build and save Phase 2 final bg_state
    p2_Prow = [p2_P_now; p2_P_now];
    p2_Orow = [O(1,:); O(1,:)];
    p2_bg_state = build_bg_state(p2_Prow, p2_Orow, ...
        [p2_RexmyI; p2_RexmyI], [p2_NavlI; p2_NavlI], ...
        zeros(2, cc, 6), ...
        zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), ...
        zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), zeros(2, cc), ...
        id_location);
    p2_bg_state.phase2_handoff_r_het_gC_m2_d = p2_R_het_last_year;
    p2_bg_state.phase2_handoff_r_het_umol_m2_d = p2_r_het_ss;
    p2_bg_state.phase2_converged = p2_converged;
    p2_bg_state.phase2_years_completed = p2_yr_final;
    p2_bg_state.phase2_annual_rel_tol = P2_ANNUAL_REL_TOL;
    p2_bg_state.standing_nutrient_state_vers = STANDING_NUTRIENT_STATE_VERS;

    p2_bg_state.phase2_observation_constraint = BG_SPINUP_CONSTRAINT;
    p2_bg_state.phase2_phosphorus_budget = p2_p_budget;
    p2_bg_state.phase2_constraint_diag = p2_constraint_diag;
    p2_bg_state.phase2_constraint_cumulative_delta = p2_constraint_cumulative_delta;
    p2_bg_state.phase2_pool_scale_floor_g_element_m2 = P2_POOL_SCALE_FLOOR;
    p2_bg_state.phase2_pool_negative_tolerance_g_element_m2 = P2_POOL_NEG_TOL;
    p2_bg_state.phase2_streak_required = P2_STREAK_REQ;
    p2_bg_state.spinup_cycle_years = spinup_cycle_years;
    p2_bg_state.spinup_provenance = spinup_provenance;
    save_bg_state(p2_chk, p2_bg_state);
    fprintf('Phase 2 state saved to: %s\n', p2_chk);
    if isfile(p2_inprog_chk); delete(p2_inprog_chk); end
    clear p2ip p2_Prow p2_Orow
end

%% Load the latest state needed for Phase 3 if not in workspace
if ~spinup_handoff_only && (~exist('p2_bg_state', 'var') || ~isfield(p2_bg_state, 'B_out'))
    if p3_done
        p2_bg_state = p3_resume_state;
    else
        error('Run_BG_Spinup:NoPhase2', 'Phase 2 state not found; run from cold start.');
    end
end

%% =========================================================================
%% Phase 3: Coupled no-ERW loop until dynamic equilibrium
%% =========================================================================
if spinup_handoff_only
    fprintf('\n=== Phase 3: handoff only, no convergence loop ===\n');
else
    fprintf('\n=== Phase 3: Coupled no-ERW convergence loop ===\n');
end

if isnan(transient_start) || isnan(transient_end)
    error('Run_BG_Spinup:NoTransientWindow', ...
        'transient_start/transient_end must be set in MOD_PARAM for Phase 3.');
end
[p3_met, ~] = Load_Spinup_Meteo(daily_meteo_path, transient_start, transient_end, ...
    Lat, Lon, DeltaGMT);

% Initialise BG state for Phase 3 from Phase 2 output (a handoff-only run
% already holds its state, loaded from SPINUP_SEED_FILE)
if ~spinup_handoff_only
    B_out      = p2_bg_state.B_out;
    O_out      = p2_bg_state.O_out;
    NavlI_out  = p2_bg_state.NavlI_out;
    RexmyI_out = p2_bg_state.RexmyI_out;
    NupI_L_out = p2_bg_state.NupI_L_out;
    if isfield(p2_bg_state,'Veg_out') && isstruct(p2_bg_state.Veg_out)
        Veg_out = p2_bg_state.Veg_out;
    end
    if isfield(p2_bg_state,'Nreserve_L_out');  Nreserve_L_out  = p2_bg_state.Nreserve_L_out;  end
    if isfield(p2_bg_state,'Preserve_L_out');  Preserve_L_out  = p2_bg_state.Preserve_L_out;  end
    if isfield(p2_bg_state,'Kreserve_L_out');  Kreserve_L_out  = p2_bg_state.Kreserve_L_out;  end
    if isfield(p2_bg_state,'Careserve_L_out'); Careserve_L_out = p2_bg_state.Careserve_L_out; end
    if isfield(p2_bg_state,'Mgreserve_L_out'); Mgreserve_L_out = p2_bg_state.Mgreserve_L_out; end
    if isfield(p2_bg_state,'Sireserve_L_out'); Sireserve_L_out = p2_bg_state.Sireserve_L_out; end
    if isfield(p2_bg_state,'Nreserve_H_out');  Nreserve_H_out  = p2_bg_state.Nreserve_H_out;  end
    if isfield(p2_bg_state,'Preserve_H_out');  Preserve_H_out  = p2_bg_state.Preserve_H_out;  end
    if isfield(p2_bg_state,'Kreserve_H_out');  Kreserve_H_out  = p2_bg_state.Kreserve_H_out;  end
    if isfield(p2_bg_state,'Careserve_H_out'); Careserve_H_out = p2_bg_state.Careserve_H_out; end
    if isfield(p2_bg_state,'Mgreserve_H_out'); Mgreserve_H_out = p2_bg_state.Mgreserve_H_out; end
    if isfield(p2_bg_state,'Sireserve_H_out'); Sireserve_H_out = p2_bg_state.Sireserve_H_out; end
end

% Resume from Phase 3 checkpoint if present
p3_loop_start = 1;
P3_SLOW_ANNUAL_REL_TOL_YR = P3_SLOW_RING_REL_TOL/n_spinup_years; % [yr^-1]
p3_converged  = false;
p3_B_prev     = [];
p3_B_prev_by_ring = cell(1, n_spinup_years);
p3_pool_ring_scale = [];
p3_pass_by_ring = false(1, n_spinup_years);
p3_ring_streak = 0;
p3_loop_last = p3_loop_start - 1;
p3_dynamic_max_by_ring = nan(1, n_spinup_years);
p3_slow_max_by_ring = nan(1, n_spinup_years);
%%% Diagnostic only, never gated: dB46/dt = K3p*B45 and there is no return
%%% flux. If the other P pools are periodic while B45 remains positive, pool 46
%%% still grows. Report full total-P drift per ring, but do not require it.
p3_total_P_rel_by_ring = nan(1, n_spinup_years);
p3_nutrient_budgets = struct();
%%% Each ring year's ERW water loss [L m-2 d-1] (ERW_Mean_Water_Loss); the
%%% ring's mean is handed off as the long-term water loss of the rain Na seed.
p3_water_loss_by_ring = nan(1, n_spinup_years);

if p3_done
    tmp = p3_resume_state;
    p3_loop_start = tmp.p3_loop + 1;
    p3_loop_last  = tmp.p3_loop;
    p3_converged = isfield(tmp,'p3_converged') && logical(tmp.p3_converged);
    p3_pass_by_ring = logical(tmp.p3_pass_by_ring);
    p3_ring_streak = tmp.p3_ring_streak;
    if isfield(tmp,'p3_B_prev_by_ring')
        p3_B_prev_by_ring = tmp.p3_B_prev_by_ring;
    end
    p3_pool_ring_scale = tmp.p3_pool_ring_scale;
    if isfield(tmp,'p3_dynamic_max_by_ring'); p3_dynamic_max_by_ring = tmp.p3_dynamic_max_by_ring; end
    if isfield(tmp,'p3_slow_max_by_ring'); p3_slow_max_by_ring = tmp.p3_slow_max_by_ring; end
    if isfield(tmp,'p3_total_P_rel_by_ring'); p3_total_P_rel_by_ring = tmp.p3_total_P_rel_by_ring; end
    if isfield(tmp,'p3_water_loss_by_ring'); p3_water_loss_by_ring = tmp.p3_water_loss_by_ring; end
    if isfield(tmp,'p3_nutrient_budgets'); p3_nutrient_budgets = tmp.p3_nutrient_budgets; end
    B_out         = tmp.B_out;
    O_out         = tmp.O_out;
    NavlI_out     = tmp.NavlI_out;
    RexmyI_out    = tmp.RexmyI_out;
    NupI_L_out    = tmp.NupI_L_out;
    if isfield(tmp,'Veg_out') && isstruct(tmp.Veg_out)
        Veg_out = tmp.Veg_out;
    end
    if isfield(tmp,'Nreserve_L_out');  Nreserve_L_out  = tmp.Nreserve_L_out;  end
    if isfield(tmp,'Preserve_L_out');  Preserve_L_out  = tmp.Preserve_L_out;  end
    if isfield(tmp,'Kreserve_L_out');  Kreserve_L_out  = tmp.Kreserve_L_out;  end
    if isfield(tmp,'Careserve_L_out'); Careserve_L_out = tmp.Careserve_L_out; end
    if isfield(tmp,'Mgreserve_L_out'); Mgreserve_L_out = tmp.Mgreserve_L_out; end
    if isfield(tmp,'Sireserve_L_out'); Sireserve_L_out = tmp.Sireserve_L_out; end
    if isfield(tmp,'Nreserve_H_out');  Nreserve_H_out  = tmp.Nreserve_H_out;  end
    if isfield(tmp,'Preserve_H_out');  Preserve_H_out  = tmp.Preserve_H_out;  end
    if isfield(tmp,'Kreserve_H_out');  Kreserve_H_out  = tmp.Kreserve_H_out;  end
    if isfield(tmp,'Careserve_H_out'); Careserve_H_out = tmp.Careserve_H_out; end
    if isfield(tmp,'Mgreserve_H_out'); Mgreserve_H_out = tmp.Mgreserve_H_out; end
    if isfield(tmp,'Sireserve_H_out'); Sireserve_H_out = tmp.Sireserve_H_out; end
    fprintf('  Resumed Phase 3 from loop %d\n', p3_loop_start);
    bg_state = tmp;
    clear tmp
end
if spinup_handoff_only
    %%% The reconstructed state carries the convergence record of the
    %%% replayed cell it was built from (Reconstruct_Cell_Pools_RF); its own
    %%% ERW handoff runs below under this cell's inputs.
    p3_handoff_source = bg_state;
    if isfield(p3_handoff_source,'Veg_out') && isstruct(p3_handoff_source.Veg_out)
        Veg_out = p3_handoff_source.Veg_out;
    end
    required_record = {'spinup_converged','spinup_cycles_completed', ...
        'spinup_ring_streak_achieved','spinup_nutrient_budgets'};
    if ~all(isfield(p3_handoff_source,required_record)) || ...
            ~logical(p3_handoff_source.spinup_converged)
        error('Run_BG_Spinup:HandoffOnlyUnconverged', ...
            'A handoff-only state must carry a converged spin-up record.');
    end
    %%% The plant reserves are part of the state handed on (build_bg_state
    %%% below); a reconstructed state takes them from its template.
    reserve_fields = strcat({'N','P','K','Ca','Mg','Si','N','P','K','Ca','Mg','Si'}, ...
        'reserve_',{'L','L','L','L','L','L','H','H','H','H','H','H'},'_out');
    if ~all(isfield(p3_handoff_source,reserve_fields))
        error('Run_BG_Spinup:HandoffOnlyNeedsState', ...
            'A handoff-only state must carry the plant reserves (%s).', ...
            strjoin(reserve_fields(~isfield(p3_handoff_source,reserve_fields)),', '));
    end
    clear required_record reserve_fields
    p3_converged = true;
    p3_loop_last = p3_handoff_source.spinup_cycles_completed;
    p3_ring_streak = p3_handoff_source.spinup_ring_streak_achieved;
    p3_pass_by_ring = true(1, n_spinup_years);
    p3_nutrient_budgets = p3_handoff_source.spinup_nutrient_budgets;
end

B_out = Project_Small_Negative_BG_Pools(B_out, P3_POOL_NEG_TOL);
for ring_index = 1:numel(p3_B_prev_by_ring)
    if ~isempty(p3_B_prev_by_ring{ring_index})
        p3_B_prev_by_ring{ring_index} = Project_Small_Negative_BG_Pools( ...
            p3_B_prev_by_ring{ring_index}, P3_POOL_NEG_TOL);
    end
end
clear ring_index

OPT_SoilBiogeochemistry_saved = OPT_SoilBiogeochemistry;
OPT_SoilBiogeochemistry = 1;
OPT_pH                  = 0;
%%% The v0.9 driver gates its ERW chemistry on OPT_ERW (OPT_pH is a Model-repo
%%% switch it does not read), so the no-ERW loop clears OPT_ERW itself: with it
%%% set, every Phase-3 year would solve the daily ERW chemistry and reseed the
%%% cation pools from it.  The final handoff diagnostic switches it back on.
OPT_ERW_saved_p3        = OPT_ERW;
OPT_ERW                 = 0;
PH                      = PHs;
%%% Every ring year and the handoff extension lie in a tracked cell's record.
if lateral_replay_spinup
    LATERAL_REPLAY(1).required = true;
end

p3_keep = [who; { ...
    'p3_keep'; 'p3_loop'; 'p3_loop_last'; 'p3_loop_start'; 'p3_converged'; 'p3_B_prev'; 'p3_B_prev_by_ring'; ...
    'p3_pool_ring_scale'; 'p3_pass_by_ring'; 'p3_ring_streak'; 'p3_last_convergence_diag'; ...
    'p3_met'; 'p3_chk'; 'p1_chk'; 'p2_chk'; 'p2_inprog_chk'; ...
    'P3_RING_REL_TOL'; 'P3_SLOW_ANNUAL_REL_TOL_YR'; 'P3_POOL_SCALE_FLOOR'; ...
    'P3_POOL_NEG_TOL'; 'P3_RING_STREAK_REQ'; ...
    'OPT_SoilBiogeochemistry_saved'; 'OPT_ERW_saved_p3'; ...
    'B_out'; 'O_out'; 'NavlI_out'; 'RexmyI_out'; 'NupI_L_out'; ...
    'Nreserve_L_out'; 'Preserve_L_out'; 'Kreserve_L_out'; 'Careserve_L_out'; ...
    'Mgreserve_L_out'; 'Sireserve_L_out'; ...
    'Nreserve_H_out'; 'Preserve_H_out'; 'Kreserve_H_out'; 'Careserve_H_out'; ...
    'Mgreserve_H_out'; 'Sireserve_H_out'; ...
    'Bio_Zs'; 'Zbio'; 'rsd'; 'Veg_out'; ...
    'p3_nutrient_budgets'}];
p3_keep = setdiff(p3_keep, {'Mpar_L','Mpar_H','ParEx_L','ParEx_H','VegH_Param_Dyn','VegL_Param_Dyn'}, 'stable');

if spinup_handoff_only
    fprintf('Handoff only: the state handed in carries a converged record (%d coupled years)\n', ...
        p3_loop_last);
elseif p3_converged
    fprintf('Phase 3 checkpoint is already converged at loop %d\n', p3_loop_last);
end

p3_loop = p3_loop_start - 1;
while ~p3_converged
    p3_loop = p3_loop + 1;
    % Ring-cycle over complete pre-transient years (matches Phase 1's own
    % averaging window: the available-forcing-period methodology runs
    % this loop "as a loop until reaching the steady state of a dynamic
    % equilibrium" over full simulation years, not a single repeated partial
    % crop-season window.
    ring_pos = mod(p3_loop - 1, n_spinup_years) + 1;
    cycle_year = spinup_cycle_years(ring_pos);
    if ring_pos == 1
        p3_pass_by_ring(:) = false;
        p3_dynamic_max_by_ring(:) = NaN;
        p3_slow_max_by_ring(:) = NaN;
        p3_total_P_rel_by_ring(:) = NaN;
        if all(~cellfun(@isempty, p3_B_prev_by_ring))
            p3_pool_ring_scale = spinup_pool_ring_scale(p3_B_prev_by_ring, P3_POOL_SCALE_FLOOR);
        else
            p3_pool_ring_scale = [];
        end
    end
    fprintf('  Phase 3 loop %d (meteo %d)\n', p3_loop, cycle_year);
    yr_mask = (spinup_all_yr == cycle_year);
    Date = spinup_all_met.Date(yr_mask);  Pr = spinup_all_met.Pr(yr_mask);  Ta = spinup_all_met.Ta(yr_mask);
    Ws   = spinup_all_met.Ws(yr_mask);    ea = spinup_all_met.ea(yr_mask);  Pre = spinup_all_met.Pre(yr_mask);
    SAD1 = spinup_all_met.SAD1(yr_mask);  SAD2 = spinup_all_met.SAD2(yr_mask);
    SAB1 = spinup_all_met.SAB1(yr_mask);  SAB2 = spinup_all_met.SAB2(yr_mask);
    N    = spinup_all_met.N(yr_mask);     Tdew = spinup_all_met.Tdew(yr_mask);  esat = spinup_all_met.esat(yr_mask);
    PARB = spinup_all_met.PARB(yr_mask);  PARD = spinup_all_met.PARD(yr_mask);
    NN   = numel(Date);
    Ws(Ws <= 0) = 0.01;
    Ds   = max(esat - ea, 0);
    [YE, MO, DA, HO, ~, ~] = datevec(Date);
    Datam = [YE, MO, DA, HO];
    Ca = Load_CO2_Record(Date, Input_Directory, co2_record_file);
    clear YE MO DA HO

    PARAM_IC     = fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']);
    ERW_PARAM_IC = fullfile(Parameter_Directory, 'ERW_MOD_PARAM_excel.m');

    clear Mpar_L Mpar_H ParEx_L ParEx_H VegH_Param_Dyn VegL_Param_Dyn
    %% Series the driver creates only when absent (irrigation, wetland
    %% level): one left from another year has that year's length.
    clear IrD Wlev
    cd(Directory);
    MAIN_FRAME_BG_ERW_10mins;
    cd(current_directory);
    if leach_layer
        %%% v1.0, leaching_water 'layer': the ERW's drainage is the layer's.
        p3_water_loss_by_ring(ring_pos) = ERW_Mean_Water_Loss(T_H, T_L, Lk_bio, Bio_Zs, RfH_Zs, RfL_Zs);
    else
        p3_water_loss_by_ring(ring_pos) = ERW_Mean_Water_Loss(T_H, T_L, Lk, Bio_Zs, RfH_Zs, RfL_Zs);
    end
    p3_rows = 2:j;
    %%% The plants' nutrient stock on day 1, the state the parameter file sets.
    %%% The frame fills its tissue records (TNIT ... TSIL) from day 2, so their
    %%% row 1 is empty and a stand that starts with wood would seem to gain its
    %%% whole stock on day 2; this is that stock, [N P K Ca Mg Si] per m2.
    %%% The frame hands each day's litter to the soil the next day (BIOGEO_UNIT
    %%% takes ISOIL(j-1)), so a loop's last litter is still on its way at the end;
    %%% the soil storage counts it (the first row's is zero: the frame starts
    %%% without any).  idx are an element's three ISOIL litter columns.
    p3_transit = @(r,w,idx) sum(w.*(sum(ISOIL_H(r,:,idx),3) + sum(ISOIL_L(r,:,idx),3)));
    p3_plant1 = zeros(1,6);
    for p3_c = 1:cc
        p3_plant1 = p3_plant1+CcrownFIX(p3_c)*( ...
            Plant_Tissue_Nutrients(B_H(1,p3_c,:),Stoich_H(p3_c),NuLit_H(1,p3_c,:))+ ...
            Plant_Tissue_Nutrients(B_L(1,p3_c,:),Stoich_L(p3_c),NuLit_L(1,p3_c,:)));
    end
    p3_cover = Ccrown_t(p3_rows,:);
    p3_litter_P = sum(ISOIL_H(p3_rows,:,13:15),3) + sum(ISOIL_L(p3_rows,:,13:15),3);
    p3_reserve_harvest_P = reshape(AgrHarNut_H(p3_rows,:,2),numel(p3_rows),cc) + ...
        reshape(AgrHarNut_L(p3_rows,:,2),numel(p3_rows),cc);
    p3_p_terms.fertilizer = sum(P_inputs(p3_rows,1));
    p3_p_terms.deposition = sum(P_inputs(p3_rows,2));
    p3_p_terms.uplift = sum(P_inputs(p3_rows,3));
    p3_p_terms.planting = sum(p3_cover.*(PlantingP_H(p3_rows,:) + PlantingP_L(p3_rows,:)),'all');
    %%% Phase 3 simulates the plants, so litter return is an internal transfer
    %%% and is already carried by the harvest and plant-storage terms.
    p3_p_terms.litter = 0;
    p3_p_terms.harvest = sum(p3_cover.*(TexP_H(p3_rows,:) + TexP_L(p3_rows,:) - ...
        p3_litter_P + p3_reserve_harvest_P),'all');
    p3_p_terms.leaching = sum(LEAK_P(p3_rows) + LEAK_DOP(p3_rows));
    p3_p_terms.soil_start = sum(P(1,35:47)) + p3_transit(1,CcrownFIX,13:15);
    p3_p_terms.soil_end = sum(P(j,35:47)) + p3_transit(j,Ccrown_t(j,:),13:15);
    p3_p_terms.plant_start = p3_plant1(2)+sum(CcrownFIX.*(Preserve_H(1,:) + Preserve_L(1,:)));
    p3_p_terms.plant_end = sum(Ccrown_t(j,:).*(TPHO_H(j,:) + TPHO_L(j,:) + ...
        Preserve_H(j,:) + Preserve_L(j,:)));
    if lateral_replay_spinup
        %%% The tracked cell's lateral exchange (Lateral_Replay_Solute's
        %%% ledger) is an input and an output of each element's budget.
        p3_lat_in = sum(lat_rep.solute_in_g(p3_rows,:),1);
        p3_lat_out = sum(lat_rep.solute_out_g(p3_rows,:),1);
        p3_lat_is = ismember(lat_rep.pool_idx,35:47);
        p3_p_terms.lateral_in = sum(p3_lat_in(p3_lat_is));
        p3_p_terms.lateral_out = sum(p3_lat_out(p3_lat_is));
        fprintf(['  lateral replay: %.4g mm water applied, %.4g mm refused by ' ...
            'the moisture bounds\n'],sum(lat_rep.water_mm),sum(lat_rep.water_clip_mm));
    end
    p3_p_budget = Phase3_Phosphorus_Budget( ...
        p3_p_terms, P3_NUTRIENT_BUDGET_ABS_TOL, P3_NUTRIENT_BUDGET_REL_TOL);
    p3_p_budget.terms = p3_p_terms;
    p3_nutrient_budgets = struct();
    p3_nutrient_budgets.P = p3_p_budget;
    p3_names = {'N','K','Ca','Mg','Si'};
    p3_soil_idx = {23:34,48:55,56:63,64:71,72:79};
    %%% v1.0: BIOGEO_UNIT adds the day's B_IO.FertN rate directly and reports it as FertN_d.
    p3_fert = {FertN_d(p3_rows), ...
        Spinup_Fertilizer_Series(B_IO.FertK,jDay(p3_rows),cycle_year,'Krate'), ...
        Spinup_Fertilizer_Series(B_IO.FertCa,jDay(p3_rows),cycle_year,'Carate'), ...
        Spinup_Fertilizer_Series(B_IO.FertMg,jDay(p3_rows),cycle_year,'Mgrate'), ...
        Spinup_Fertilizer_Series(B_IO.FertSi,jDay(p3_rows),cycle_year,'Sirate')};
    p3_dep = [B_IO.DepN B_IO.DepK B_IO.DepCa B_IO.DepMg B_IO.DepSi];
    %%% With an event record the deposition varies day by day (Deposition_Day);
    %%% the budget takes what the frame applied (Dep_d: N P K Ca Mg Si).
    p3_event_dep = isfield(B_IO,'Dep_events') && ~isempty(B_IO.Dep_events);
    if p3_event_dep
        p3_dep_sum = sum(Dep_d(p3_rows,[1 3 4 5 6]),1);
    end
    p3_uplift = [0 B_IO.Tup_K B_IO.Tup_Ca B_IO.Tup_Mg B_IO.Tup_Si];
    p3_tex_h = {TexN_H,TexK_H,TexCa_H,TexMg_H,TexSi_H};
    p3_tex_l = {TexN_L,TexK_L,TexCa_L,TexMg_L,TexSi_L};
    p3_litter_idx = {10:12,16:18,19:21,22:24,25:27};
    p3_leach = {LEAK_NH4+LEAK_NO3+LEAK_DON,LEAK_K,LEAK_Ca,LEAK_Mg,LEAK_Si};
    p3_plant_h = {TNIT_H,TPOT_H,TCAL_H,TMAG_H,TSIL_H};
    p3_plant_l = {TNIT_L,TPOT_L,TCAL_L,TMAG_L,TSIL_L};
    p3_reserve_h = {Nreserve_H,Kreserve_H,Careserve_H,Mgreserve_H,Sireserve_H};
    p3_reserve_l = {Nreserve_L,Kreserve_L,Careserve_L,Mgreserve_L,Sireserve_L};
    %%% P is budgeted above; remap the remaining five nutrient columns.
    p3_nut_idx = [1 3 4 5 6];
    for p3_x = 1:numel(p3_names)
        p3_terms = struct();
        p3_n_i = p3_nut_idx(p3_x);
        if any(~isfinite(p3_fert{p3_x}))
            error('Run_BG_Spinup:InvalidFertilizerLedger', ...
                'Phase 3 %s fertilizer ledger contains a non-finite value.', ...
                p3_names{p3_x});
        end
        p3_terms.fertilizer = sum(p3_fert{p3_x});
        if p3_event_dep
            p3_terms.deposition = p3_dep_sum(p3_x);
        else
            p3_terms.deposition = numel(p3_rows)*p3_dep(p3_x);
        end
        p3_terms.uplift = numel(p3_rows)*p3_uplift(p3_x);
        p3_terms.planting = sum(p3_cover.*(PlantingNut_H(p3_rows,:,p3_n_i) + PlantingNut_L(p3_rows,:,p3_n_i)),'all');
        p3_terms.fixation = strcmp(p3_names{p3_x},'N')*sum(BfixN(p3_rows));
        p3_litter = sum(ISOIL_H(p3_rows,:,p3_litter_idx{p3_x}),3) + sum(ISOIL_L(p3_rows,:,p3_litter_idx{p3_x}),3);
        p3_reserve_harvest = reshape(AgrHarNut_H(p3_rows,:,p3_n_i),numel(p3_rows),cc) + reshape(AgrHarNut_L(p3_rows,:,p3_n_i),numel(p3_rows),cc);
        p3_terms.harvest = sum(p3_cover.*(p3_tex_h{p3_x}(p3_rows,:) + p3_tex_l{p3_x}(p3_rows,:) - p3_litter + p3_reserve_harvest),'all');
        p3_terms.leaching = sum(p3_leach{p3_x}(p3_rows));
        %%% v1.0 gas partition (BIOGEOCHEMISTRY_DYNAMIC3): N2 from denitrification and
        %%% total N2O (nitrification + denitrification); there is no NO pathway.
        p3_terms.gaseous_loss = strcmp(p3_names{p3_x},'N')*sum(VOL(p3_rows)+N2flx(p3_rows)+N2Oflx(p3_rows));
        p3_terms.soil_start = sum(P(1,p3_soil_idx{p3_x})) + p3_transit(1,CcrownFIX,p3_litter_idx{p3_x});
        p3_terms.soil_end = sum(P(j,p3_soil_idx{p3_x})) + p3_transit(j,Ccrown_t(j,:),p3_litter_idx{p3_x});
        p3_terms.plant_start = p3_plant1(p3_n_i)+sum(CcrownFIX.*(p3_reserve_h{p3_x}(1,:) + p3_reserve_l{p3_x}(1,:)));
        p3_terms.plant_end = sum(Ccrown_t(j,:).*(p3_plant_h{p3_x}(j,:) + p3_plant_l{p3_x}(j,:) + p3_reserve_h{p3_x}(j,:) + p3_reserve_l{p3_x}(j,:)));
        if lateral_replay_spinup
            p3_lat_is = ismember(lat_rep.pool_idx,p3_soil_idx{p3_x});
            p3_terms.lateral_in = sum(p3_lat_in(p3_lat_is));
            p3_terms.lateral_out = sum(p3_lat_out(p3_lat_is));
        end
        p3_nutrient_budgets.(p3_names{p3_x}) = Phase3_Nutrient_Budget( ...
            p3_terms,P3_NUTRIENT_BUDGET_ABS_TOL,P3_NUTRIENT_BUDGET_REL_TOL);
        p3_nutrient_budgets.(p3_names{p3_x}).terms = p3_terms;
        if ~p3_nutrient_budgets.(p3_names{p3_x}).closed
            error('Run_BG_Spinup:NutrientBudgetNotClosed', ...
                ['Phase 3 loop %d %s residual %.6g g m-2 exceeds %.6g (inputs: fertiliser %.6g, ' ...
                 'deposition %.6g, uplift %.6g, planting %.6g, fixation %.6g; outputs: harvest %.6g, ' ...
                 'leaching %.6g, gas %.6g; soil %.6g -> %.6g; plant %.6g -> %.6g).'], ...
                p3_loop,p3_names{p3_x},p3_nutrient_budgets.(p3_names{p3_x}).residual, ...
                p3_nutrient_budgets.(p3_names{p3_x}).tolerance,p3_terms.fertilizer, ...
                p3_terms.deposition,p3_terms.uplift,p3_terms.planting,p3_terms.fixation, ...
                p3_terms.harvest,p3_terms.leaching,p3_terms.gaseous_loss,p3_terms.soil_start, ...
                p3_terms.soil_end,p3_terms.plant_start,p3_terms.plant_end);
        end
    end
    if ~p3_p_budget.closed
        error('Run_BG_Spinup:PhosphorusBudgetNotClosed', ...
            ['Phase 3 loop %d phosphorus residual %.6g g P m-2 exceeds ' ...
            'the %.6g g P m-2 closure tolerance.'], ...
            p3_loop, p3_p_budget.residual, p3_p_budget.tolerance);
    end

    bg_state = build_bg_state(P, O, RexmyI, NavlI, NupI_L, ...
        Nreserve_L, Preserve_L, Kreserve_L, Careserve_L, Mgreserve_L, Sireserve_L, ...
        Nreserve_H, Preserve_H, Kreserve_H, Careserve_H, Mgreserve_H, Sireserve_H, id_location);
    %%% A perennial stand (its vegetation file sets PFT_perennial) goes on
    %%% from this year's end state (Perennial_Vegetation_Capture).
    if exist('PFT_perennial','var') && isequal(PFT_perennial,true)
        Perennial_Vegetation_Capture;
        bg_state.Veg_out = Veg_out;
        fprintf('%s\n',Perennial_Stand_Summary(LAI_H(2:j,:),LAI_L(2:j,:),Veg_out));
    end

    % The next year restarts from this year's END state, so a year boundary
    % moves no mass and the converged ring is a periodic orbit of the model.
    % The last year's mean (Zhang et al. 2025 hand the transient the annual
    % mean) is the ring-position convergence metric only; the handoff
    % extension below hands over its own end state.
    ok = any(isfinite(P) & abs(P) > 0, 2);
    daily_idx = find(ok, 1, 'last');
    window_start = max(1, daily_idx - 364);
    P_avg = mean(P(window_start:daily_idx, :), 1, 'omitnan');
    p3_B_metric = Project_Small_Negative_BG_Pools(P_avg, P3_POOL_NEG_TOL);
    bg_state.B_out = P(daily_idx, :);

    B_out = Project_Small_Negative_BG_Pools(bg_state.B_out, P3_POOL_NEG_TOL);
    bg_state.B_out = B_out;
    O_out       = bg_state.O_out;
    NavlI_out   = bg_state.NavlI_out;
    RexmyI_out  = bg_state.RexmyI_out;
    NupI_L_out  = bg_state.NupI_L_out;
    if isfield(bg_state,'Nreserve_L_out');  Nreserve_L_out  = bg_state.Nreserve_L_out;  end
    if isfield(bg_state,'Preserve_L_out');  Preserve_L_out  = bg_state.Preserve_L_out;  end
    if isfield(bg_state,'Kreserve_L_out');  Kreserve_L_out  = bg_state.Kreserve_L_out;  end
    if isfield(bg_state,'Careserve_L_out'); Careserve_L_out = bg_state.Careserve_L_out; end
    if isfield(bg_state,'Mgreserve_L_out'); Mgreserve_L_out = bg_state.Mgreserve_L_out; end
    if isfield(bg_state,'Sireserve_L_out'); Sireserve_L_out = bg_state.Sireserve_L_out; end
    if isfield(bg_state,'Nreserve_H_out');  Nreserve_H_out  = bg_state.Nreserve_H_out;  end
    if isfield(bg_state,'Preserve_H_out');  Preserve_H_out  = bg_state.Preserve_H_out;  end
    if isfield(bg_state,'Kreserve_H_out');  Kreserve_H_out  = bg_state.Kreserve_H_out;  end
    if isfield(bg_state,'Careserve_H_out'); Careserve_H_out = bg_state.Careserve_H_out; end
    if isfield(bg_state,'Mgreserve_H_out'); Mgreserve_H_out = bg_state.Mgreserve_H_out; end
    if isfield(bg_state,'Sireserve_H_out'); Sireserve_H_out = bg_state.Sireserve_H_out; end

    % Convergence check
    p3_pass = false;
    p3_B_prev = p3_B_prev_by_ring{ring_pos};
    if ~isempty(p3_B_prev)
        p3_diag = spinup_pool_group_diagnostics(p3_B_metric, p3_B_prev, n_spinup_years, ...
            P3_POOL_SCALE_FLOOR, p3_pool_ring_scale, P3_UNGATED_POOLS);
        p3_diag.pool_pass = p3_diag.pool_ring_rel < P3_RING_REL_TOL;
        p3_diag.pool_pass(p3_diag.slow_pool_idx) = ...
            p3_diag.pool_ring_rel(p3_diag.slow_pool_idx) < P3_SLOW_RING_REL_TOL;
        p3_diag.pool_pass(P3_UNGATED_POOLS) = true;
        p3_dynamic_pass = p3_diag.dynamic_ring_rel < P3_RING_REL_TOL;
        p3_slow_pass = p3_diag.slow_ring_rel < P3_SLOW_RING_REL_TOL;
        p3_element_total_pass = p3_diag.element_ring_rel < P3_RING_REL_TOL;
        p3_pass = p3_diag.finite_nonnegative && all(p3_diag.pool_pass) && ...
            all(p3_dynamic_pass) && all(p3_slow_pass) && ...
            all(p3_element_total_pass);
        p3_dynamic_max_by_ring(ring_pos) = max(p3_diag.dynamic_ring_rel);
        p3_slow_max_by_ring(ring_pos) = max(p3_diag.slow_ring_rel);
        %%% Pools 35:47 (total P), against the ring-invariant pool_scale used
        %%% above so a growing pool can't pass by a shrinking relative
        %%% denominator (see the p3_total_P_rel_by_ring init comment for why
        %%% this stays a diagnostic rather than a gate).
        p3_total_P_ring_abs = abs(sum(p3_B_metric(35:47)) - sum(p3_B_prev(35:47)));
        p3_total_P_ring_rel = p3_total_P_ring_abs / sum(p3_diag.pool_scale(35:47));
        p3_total_P_pass = p3_total_P_ring_rel < P3_RING_REL_TOL;
        p3_total_P_rel_by_ring(ring_pos) = p3_total_P_ring_rel;
        p3_last_convergence_diag = struct('loop', p3_loop, 'cycle_year', cycle_year, ...
            'ring_pos', ring_pos, ...
            'ring_rel_tol', P3_RING_REL_TOL, ...
            'slow_annual_rel_tol_yr', P3_SLOW_ANNUAL_REL_TOL_YR, ...
            'pool_ring_rel', p3_diag.pool_ring_rel, ...
            'pool_ring_abs', p3_diag.pool_ring_abs, ...
            'pool_basis_idx', p3_diag.pool_basis_idx, ...
            'pool_scale', p3_diag.pool_scale, ...
            'pool_pass', p3_diag.pool_pass, ...
            'dynamic_group_ring_rel', p3_diag.dynamic_group_ring_rel, ...
            'dynamic_group_ring_abs', p3_diag.dynamic_group_ring_abs, ...
            'dynamic_pool_ring_abs', p3_diag.dynamic_pool_ring_abs, ...
            'dynamic_pool_ring_rel', p3_diag.dynamic_pool_ring_rel, ...
            'slow_group_ring_rel', p3_diag.slow_group_ring_rel, ...
            'slow_group_ring_abs', p3_diag.slow_group_ring_abs, ...
            'slow_pool_ring_abs', p3_diag.slow_pool_ring_abs, ...
            'slow_pool_ring_rel', p3_diag.slow_pool_ring_rel, ...
            'slow_group_annual_rel', p3_diag.slow_group_annual_rel, ...
            'slow_pool_annual_rel', p3_diag.slow_pool_annual_rel, ...
            'dynamic_ring_rel', p3_diag.dynamic_ring_rel, ...
            'dynamic_ring_abs', p3_diag.dynamic_ring_abs, ...
            'slow_ring_abs', p3_diag.slow_ring_abs, ...
            'slow_ring_rel', p3_diag.slow_ring_rel, ...
            'slow_annual_rel', p3_diag.slow_annual_rel, 'pass', p3_pass, ...
            'total_P_ring_rel', p3_total_P_ring_rel, 'total_P_pass', p3_total_P_pass, ...
            'element_names', {p3_diag.element_names}, ...
            'element_total_ring_rel', p3_diag.element_ring_rel, ...
            'element_total_pass', p3_element_total_pass);
        fprintf(['  Phase 3 loop %d: dynamic ring max=%.4f slow ring max=%.4f ' ...
            'slow annual max=%.4f yr^-1 ring pass=%d\n'], p3_loop, ...
            max(p3_diag.dynamic_ring_rel), max(p3_diag.slow_ring_rel), ...
            max(p3_diag.slow_annual_rel), p3_pass);
        if ~all(p3_diag.pool_pass)
            fprintf('    individual pool failures: %s\n', mat2str(find(~p3_diag.pool_pass)));
        end
        if ~all(p3_dynamic_pass)
            fprintf('    dynamic ring max driven by: %s\n', spinup_worst_offender_label( ...
                p3_diag.dynamic_ring_rel, p3_diag.dynamic_group_idx, p3_diag.dynamic_pool_idx));
        end
        if ~all(p3_slow_pass)
            fprintf('    slow ring max driven by: %s\n', spinup_worst_offender_label( ...
                p3_diag.slow_ring_rel, p3_diag.slow_group_idx, p3_diag.slow_pool_idx));
        end
    end
    p3_pass_by_ring(ring_pos) = p3_pass;
    p3_ring_pass = false;
    if ring_pos == n_spinup_years
        p3_ring_pass = all(p3_pass_by_ring);
        [p3_ring_streak, p3_converged] = spinup_update_consecutive_passes( ...
            p3_ring_pass, p3_ring_streak, P3_RING_STREAK_REQ);
        fprintf('  Phase 3 complete-ring streak=%d/%d\n', ...
            p3_ring_streak, P3_RING_STREAK_REQ);

    else
        p3_converged = false;
    end
    p3_B_prev_by_ring{ring_pos} = p3_B_metric;

    p3_loop_last = p3_loop;
    bg_state.p3_loop   = p3_loop;
    bg_state.p3_converged = p3_converged;
    bg_state.spinup_cycle_years = spinup_cycle_years;
    bg_state.p3_handoff_applied = false;
    bg_state.p3_ring_rel_tol = P3_RING_REL_TOL;

    bg_state.p3_slow_ring_rel_tol = P3_SLOW_RING_REL_TOL;
    bg_state.p3_pool_scale_floor_g_element_m2 = P3_POOL_SCALE_FLOOR;
    bg_state.p3_pool_negative_tolerance_g_element_m2 = P3_POOL_NEG_TOL;
    bg_state.p3_ring_streak = p3_ring_streak;
    bg_state.p3_ring_streak_required = P3_RING_STREAK_REQ;
    bg_state.p3_ring_pass = p3_ring_pass;
    bg_state.p3_pass_by_ring = p3_pass_by_ring;
    bg_state.p3_B_prev_by_ring = p3_B_prev_by_ring;
    bg_state.p3_pool_ring_scale = p3_pool_ring_scale;
    bg_state.p3_dynamic_max_by_ring = p3_dynamic_max_by_ring;
    bg_state.p3_slow_max_by_ring = p3_slow_max_by_ring;
    bg_state.p3_total_P_rel_by_ring = p3_total_P_rel_by_ring;
    bg_state.p3_water_loss_by_ring = p3_water_loss_by_ring;
    %%% Gated on by load_compatible_bg_checkpoint, so it must record the
    %%% constraint this loop actually ran under.
    bg_state.p3_observation_constraint = BG_SPINUP_CONSTRAINT;
    bg_state.p3_transient_observation_constraint = BG_TRANSIENT_INITIAL_CONSTRAINT;
    bg_state.p3_phosphorus_budget = p3_p_budget;
    bg_state.p3_nutrient_budget_vers = P3_NUTRIENT_BUDGET_VERS;
    bg_state.p3_nutrient_budgets = p3_nutrient_budgets;
    bg_state.p3_p_budget_abs_tol_gP_m2 = P3_NUTRIENT_BUDGET_ABS_TOL;
    bg_state.p3_p_budget_rel_tol = P3_NUTRIENT_BUDGET_REL_TOL;
    if exist('p3_last_convergence_diag', 'var')
        bg_state.p3_last_convergence_diag = p3_last_convergence_diag;
    end
    %%% The physicality gate below must keep running on a resume
    %%% whose ring loop is skipped as already converged (MAIN_FRAME never
    %%% runs, so Bio_Zs/rsd/Zbio never exist).  Persist the bio-layer soil
    %%% mass here so the resumed state can feed the gate from the checkpoint.
    bg_state.spinup_biogeochemical_soil_mass_kg_m2 = ...
        sum(Bio_Zs.*rsd)*Zbio/1000;
    bg_state.spinup_provenance = spinup_provenance;
    save_bg_state(p3_chk, bg_state);

    clearvars('-except', p3_keep{:});

    if p3_converged
        fprintf('Phase 3 converged at loop %d\n', p3_loop);
        break
    end
end

%%% A completed cycle budget is diagnostic only. The handoff requires a passed
%%% complete ring so no unresolved pool drift can seed the transient.
%%% Physical pool gate: reject an unphysical converged state before handoff.
%%% The bio-layer soil mass lives in MAIN_FRAME's workspace (Bio_Zs/rsd/Zbio),
%%% which a skip-loops-already-converged resume never produces; fall back to
%%% the checkpoint value stored above, then the workbook bulk density (equal,
%%% because rsd is layer-uniform and sum(Bio_Zs) = 1).
if exist('Bio_Zs', 'var') && exist('rsd', 'var') && exist('Zbio', 'var')
    p3_soil_mass = sum(Bio_Zs.*rsd)*Zbio/1000;
elseif isfield(bg_state, 'spinup_biogeochemical_soil_mass_kg_m2') && ...
        isscalar(bg_state.spinup_biogeochemical_soil_mass_kg_m2) && ...
        isfinite(bg_state.spinup_biogeochemical_soil_mass_kg_m2) && ...
        bg_state.spinup_biogeochemical_soil_mass_kg_m2 > 0
    p3_soil_mass = bg_state.spinup_biogeochemical_soil_mass_kg_m2;
elseif ismember('rsd', MOD_PARAM_input.Properties.VariableNames) && ...
        isfinite(MOD_PARAM_input.rsd) && isfinite(MOD_PARAM_input.Zbio) && ...
        MOD_PARAM_input.rsd > 0 && MOD_PARAM_input.Zbio > 0
    p3_soil_mass = MOD_PARAM_input.rsd * MOD_PARAM_input.Zbio / 1000;
else
    error('Run_BG_Spinup:MissingSoilMassForGate', ...
        'Cannot compute the converged Phase 3 soil mass: no MAIN_FRAME layers and no checkpoint value.');
end
Validate_BG_Pool_Physicality(B_out, p3_soil_mass, 'Phase 3 converged state');

% Handoff extension: one non-repeated pass over the real partial-year meteo
% between Jan 1 and the transient start (e.g. 2022-01-01 to 2022-05-04),
% bridging the ring-cycled complete-year convergence loop above to the
% actual observed transient window using real (not cycled) daily forcing.
is_handoff_extension = ~spinup_handoff_only && ...
    exist('handoff_end_dn', 'var') && isfinite(handoff_end_dn);
if is_handoff_extension
    fprintf('  Phase 3 handoff extension: %s - %s (real meteo)\n', ...
        datestr(handoff_start_dn, 'yyyy-mm-dd'), datestr(handoff_end_dn - 1/24, 'yyyy-mm-dd HH:MM'));
    yr_mask = (spinup_all_met.Date >= handoff_start_dn) & (spinup_all_met.Date < handoff_end_dn);
    Date = spinup_all_met.Date(yr_mask);  Pr = spinup_all_met.Pr(yr_mask);  Ta = spinup_all_met.Ta(yr_mask);
    Ws   = spinup_all_met.Ws(yr_mask);    ea = spinup_all_met.ea(yr_mask);  Pre = spinup_all_met.Pre(yr_mask);
    SAD1 = spinup_all_met.SAD1(yr_mask);  SAD2 = spinup_all_met.SAD2(yr_mask);
    SAB1 = spinup_all_met.SAB1(yr_mask);  SAB2 = spinup_all_met.SAB2(yr_mask);
    N    = spinup_all_met.N(yr_mask);     Tdew = spinup_all_met.Tdew(yr_mask);  esat = spinup_all_met.esat(yr_mask);
    PARB = spinup_all_met.PARB(yr_mask);  PARD = spinup_all_met.PARD(yr_mask);
    NN   = numel(Date);
    Ws(Ws <= 0) = 0.01;
    Ds   = max(esat - ea, 0);
    [YE, MO, DA, HO, ~, ~] = datevec(Date);
    Datam = [YE, MO, DA, HO];
    Ca = Load_CO2_Record(Date, Input_Directory, co2_record_file);
    clear YE MO DA HO

    PARAM_IC     = fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']);
    ERW_PARAM_IC = fullfile(Parameter_Directory, 'ERW_MOD_PARAM_excel.m');

    clear Mpar_L Mpar_H ParEx_L ParEx_H VegH_Param_Dyn VegL_Param_Dyn
    %% Series the driver creates only when absent (irrigation, wetland
    %% level): one left from another year has that year's length.
    clear IrD Wlev
    cd(Directory);
    MAIN_FRAME_BG_ERW_10mins;
    cd(current_directory);

    bg_state = build_bg_state(P, O, RexmyI, NavlI, NupI_L, ...
        Nreserve_L, Preserve_L, Kreserve_L, Careserve_L, Mgreserve_L, Sireserve_L, ...
        Nreserve_H, Preserve_H, Kreserve_H, Careserve_H, Mgreserve_H, Sireserve_H, id_location);
    %%% A perennial stand (its vegetation file sets PFT_perennial) goes on
    %%% from this year's end state (Perennial_Vegetation_Capture).
    if exist('PFT_perennial','var') && isequal(PFT_perennial,true)
        Perennial_Vegetation_Capture;
        bg_state.Veg_out = Veg_out;
        fprintf('%s\n',Perennial_Stand_Summary(LAI_H(2:j,:),LAI_L(2:j,:),Veg_out));
    end

    % The real pre-transient extension defines the state at the observed start.
    ok = any(isfinite(P) & abs(P) > 0, 2);
    daily_idx = find(ok, 1, 'last');
    P_handoff = P(daily_idx, :);
    bg_state.B_out = P_handoff;

    B_out = Project_Small_Negative_BG_Pools(bg_state.B_out, P3_POOL_NEG_TOL);
    bg_state.B_out = B_out;
    O_out       = bg_state.O_out;
    NavlI_out   = bg_state.NavlI_out;
    RexmyI_out  = bg_state.RexmyI_out;
    NupI_L_out  = bg_state.NupI_L_out;
    if isfield(bg_state,'Nreserve_L_out');  Nreserve_L_out  = bg_state.Nreserve_L_out;  end
    if isfield(bg_state,'Preserve_L_out');  Preserve_L_out  = bg_state.Preserve_L_out;  end
    if isfield(bg_state,'Kreserve_L_out');  Kreserve_L_out  = bg_state.Kreserve_L_out;  end
    if isfield(bg_state,'Careserve_L_out'); Careserve_L_out = bg_state.Careserve_L_out; end
    if isfield(bg_state,'Mgreserve_L_out'); Mgreserve_L_out = bg_state.Mgreserve_L_out; end
    if isfield(bg_state,'Sireserve_L_out'); Sireserve_L_out = bg_state.Sireserve_L_out; end
    if isfield(bg_state,'Nreserve_H_out');  Nreserve_H_out  = bg_state.Nreserve_H_out;  end
    if isfield(bg_state,'Preserve_H_out');  Preserve_H_out  = bg_state.Preserve_H_out;  end
    if isfield(bg_state,'Kreserve_H_out');  Kreserve_H_out  = bg_state.Kreserve_H_out;  end
    if isfield(bg_state,'Careserve_H_out'); Careserve_H_out = bg_state.Careserve_H_out; end
    if isfield(bg_state,'Mgreserve_H_out'); Mgreserve_H_out = bg_state.Mgreserve_H_out; end
    if isfield(bg_state,'Sireserve_H_out'); Sireserve_H_out = bg_state.Sireserve_H_out; end

end

%%% The recurrent spin-up is fully prognostic. At the observed transient start,
%%% the measured slow stocks are mapped once so the transient begins from the
%%% field state rather than an unobserved legacy-P trajectory.
[B_out, handoff_constraint_diag] = spinup_apply_observation_constraint( ...
    B_out, BG_TRANSIENT_INITIAL_CONSTRAINT);
B_out = Project_Small_Negative_BG_Pools(B_out, P3_POOL_NEG_TOL);

OPT_SoilBiogeochemistry = OPT_SoilBiogeochemistry_saved;

NupI_L_no_erw = NupI_L_out;
if ismatrix(NupI_L_no_erw) && size(NupI_L_no_erw,1) == cc && size(NupI_L_no_erw,2) == 6
    NupI_L_no_erw = reshape(NupI_L_no_erw, 1, cc, 6);
end
bg_state_no_erw = build_bg_state([B_out; B_out], [O_out; O_out], ...
    [RexmyI_out; RexmyI_out], [NavlI_out; NavlI_out], cat(1, NupI_L_no_erw, NupI_L_no_erw), ...
    [Nreserve_L_out; Nreserve_L_out], [Preserve_L_out; Preserve_L_out], ...
    [Kreserve_L_out; Kreserve_L_out], [Careserve_L_out; Careserve_L_out], ...
    [Mgreserve_L_out; Mgreserve_L_out], [Sireserve_L_out; Sireserve_L_out], ...
    [Nreserve_H_out; Nreserve_H_out], [Preserve_H_out; Preserve_H_out], ...
    [Kreserve_H_out; Kreserve_H_out], [Careserve_H_out; Careserve_H_out], ...
    [Mgreserve_H_out; Mgreserve_H_out], [Sireserve_H_out; Sireserve_H_out], id_location);
if exist('Veg_out','var') && isstruct(Veg_out)
    bg_state_no_erw.Veg_out = Veg_out;   % a perennial stand, for the transient
end
%%% The long-term water loss of the ERW layer, which the rain Na seed divides
%%% by (ERW_Seed_Water_Loss): the ring's mean, or a handed-in state's own.
%%% The drivers seed on their second day, when their forcing series hold one
%%% day of transpiration and drainage, so they take it from the state.
if spinup_handoff_only
    ERW_seed_water_loss_l_m2_d = NaN;
    if isfield(p3_handoff_source,'seed_water_loss_l_m2_d')
        ERW_seed_water_loss_l_m2_d = p3_handoff_source.seed_water_loss_l_m2_d;
    end
else
    ERW_seed_water_loss_l_m2_d = mean(p3_water_loss_by_ring(isfinite(p3_water_loss_by_ring)));
end
if isscalar(ERW_seed_water_loss_l_m2_d) && isfinite(ERW_seed_water_loss_l_m2_d) && ...
        ERW_seed_water_loss_l_m2_d > 0
    bg_state_no_erw.seed_water_loss_l_m2_d = ERW_seed_water_loss_l_m2_d;
elseif (exist('ERW_rain_chemistry','var') && ~isempty(strtrim(ERW_rain_chemistry.source_file))) || ...
        (exist('B_IO','var') && isfield(B_IO,'Dep_events') && ...
        isfield(B_IO.Dep_events,'erw_mean_umol_m2_d'))
    error('Run_BG_Spinup:NoSeedWaterLoss', ...
        ['The rain or deposition record seeds Na from the long-term water loss, ' ...
         'which this spin-up has not recorded: rerun Phase 3.']);
else
    clear ERW_seed_water_loss_l_m2_d
end
P_no_erw = bg_state_no_erw.B_out;
O_no_erw = bg_state_no_erw.O_out;

%% =========================================================================
%% Phase 3 final: ERW charge-balance diagnostic
%% =========================================================================
fprintf('\n=== Phase 3 final: ERW charge-balance diagnostic ===\n');
% The repeated phases converge the unamended background with the ERW
% chemistry off, even when the run is configured as erw.  Turn it on here for
% the charge-balance check of the seed the transient will build.
OPT_pH = 1;
OPT_ERW = 1;
% The biogeochemistry module and inversion target both use observed bulk pH;
% pore-water pH remains an ERW equilibrium diagnostic.
PH = PHs;
spinup_erw_readiness = ERW_Handoff_Readiness();
clear erw_readiness_diag
% The driver seeds the ERW state on its second day (the first daily
% boundary) from the first two days of forcing and the handoff pools, exactly
% as the transient will, so those two days are all this pass needs.
ERW_SEED_HOURS = 48;
if numel(p3_met.Date) < ERW_SEED_HOURS
    error('Run_BG_Spinup:ShortSeedWindow', ...
        'The transient window holds %d forcing hours; the ERW seed needs %d.', ...
        numel(p3_met.Date), ERW_SEED_HOURS);
end
seed = 1:ERW_SEED_HOURS;
Date = p3_met.Date(seed);  Pr = p3_met.Pr(seed);    Ta = p3_met.Ta(seed);
Ws   = p3_met.Ws(seed);    ea = p3_met.ea(seed);    Pre = p3_met.Pre(seed);
SAD1 = p3_met.SAD1(seed);  SAD2 = p3_met.SAD2(seed);
SAB1 = p3_met.SAB1(seed);  SAB2 = p3_met.SAB2(seed);
N    = p3_met.N(seed);     Tdew = p3_met.Tdew(seed);  esat = p3_met.esat(seed);
PARB = p3_met.PARB(seed);  PARD = p3_met.PARD(seed);
clear seed
NN   = numel(Date);
Ws(Ws <= 0) = 0.01;
Ds   = max(esat - ea, 0);
[YE, MO, DA, HO, ~, ~] = datevec(Date);
Datam = [YE, MO, DA, HO];
Ca = Load_CO2_Record(Date, Input_Directory, co2_record_file);
clear YE MO DA HO

PARAM_IC     = fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']);
ERW_PARAM_IC = fullfile(Parameter_Directory, 'ERW_MOD_PARAM_excel.m');

clear Mpar_L Mpar_H ParEx_L ParEx_H VegH_Param_Dyn VegL_Param_Dyn
clear ERW_biogeochem ERW_weathering ERW_const ERW_TC ERW_TC_const
%%% The seed check runs the transient's first two days, which no lateral
%%% record covers: a replayed column keeps its terrain radiation there and
%%% receives no lateral exchange.
if lateral_replay_spinup
    LATERAL_REPLAY(1).required = false;
end
try
    %% Series the driver creates only when absent (irrigation, wetland
    %% level): one left from another year has that year's length.
    clear IrD Wlev
    cd(Directory);
    MAIN_FRAME_BG_ERW_10mins;
    cd(current_directory);
    %%% Readiness of the seed the driver built: its initial ERW state, the
    %%% handoff pools it was seeded from, and the dissolved-inventory water of
    %%% the seeding day (the driver's own cov_moll at j = 2).
    if ~(exist('ERW_biogeochem','var') && isstruct(ERW_biogeochem) && ...
            exist('ERW_TC','var') && isstruct(ERW_TC) && isfield(ERW_TC,'s') && ...
            size(ERW_TC.s,1) >= 2)
        error('Run_BG_Spinup:NoERWSeed', ...
            'The two-day pass did not build an ERW seed (is OPT_ERW honoured?).');
    end
    cov_moll_seed = ERW_TC_const.n * Zbio / 1000 * ERW_TC.s(2,:) * 1000;
    erw_readiness_diag = ERW_Handoff_Diagnostic(P(1,:), ERW_biogeochem, ...
        ERW_weathering, ERW_const, cov_moll_seed);
    fprintf(['  ERW handoff diagnostic: pH_pore=%.2f, pH_bulk=%.2f ' ...
        '(forward %.2f, observed %.2f), base_sat=%.3f, Al=%.3g umol/L, ' ...
        'charge=%.3e umolc/L (rel %.3e), f_CEC sum=%.6f\n'], ...
        erw_readiness_diag.pH_pore_initial, erw_readiness_diag.pH_bulk_initial, ...
        erw_readiness_diag.pH_bulk_predicted, erw_readiness_diag.pH_bulk_observed, ...
        erw_readiness_diag.base_saturation, erw_readiness_diag.dissolved_Al_umol_L, ...
        erw_readiness_diag.charge_residual_umolc_L, ...
        erw_readiness_diag.relative_charge_residual, erw_readiness_diag.f_CEC_sum);
    if erw_readiness_diag.status == "warning"
        warning('Run_BG_Spinup:ERWHandoffDiagnostic', ...
            'ERW handoff diagnostic warnings: %s', ...
            char(strjoin(erw_readiness_diag.warnings, '; ')));
    end
catch ME
    cd(current_directory);
    diagnostic_error = MException('Run_BG_Spinup:ERWHandoffDiagnostic', ...
        'Final ERW handoff diagnostic failed.');
    throwAsCaller(addCause(diagnostic_error, ME));
end
OPT_ERW = OPT_ERW_saved_p3;
if exist('erw_readiness_diag','var')
    spinup_erw_readiness = ERW_Handoff_Readiness(erw_readiness_diag);
    fprintf('  ERW handoff readiness: %s - %s\n', ...
        char(spinup_erw_readiness.status), spinup_erw_readiness.message);
    if ~spinup_erw_readiness.ready
        error('Run_BG_Spinup:ERWHandoffNotReady', '%s', spinup_erw_readiness.message);
    end
else
    error('Run_BG_Spinup:MissingERWHandoffDiagnostic', ...
        'Final ERW handoff diagnostic did not produce readiness fields.');
end
P = P_no_erw;
O = O_no_erw;
bg_state = bg_state_no_erw;

%% Save final bg_state
%%% A completed spin-up has one passed meteorological ring and a ready ERW handoff.
bg_state = Attach_Spinup_Metadata(bg_state, p3_converged, ...
    p3_converged, ...
    spinup_cycle_years(1), spinup_cycle_years(end), p3_loop_last, ...
    NaN, NaN, NaN, NaN, ...
    {}, [], [], P3_RING_REL_TOL, P3_RING_REL_TOL, NaN, '', p3_converged, ...
    n_spinup_years, sum(p3_pass_by_ring), spinup_erw_readiness);
bg_state.spinup_ring_rel_tol = P3_RING_REL_TOL;

bg_state.spinup_slow_ring_rel_tol = P3_SLOW_RING_REL_TOL;
bg_state.spinup_slow_annual_rel_tol_yr = P3_SLOW_ANNUAL_REL_TOL_YR;
bg_state.spinup_pool_scale_floor_g_element_m2 = P3_POOL_SCALE_FLOOR;
bg_state.spinup_pool_negative_tolerance_g_element_m2 = P3_POOL_NEG_TOL;
bg_state.spinup_ring_streak_required = P3_RING_STREAK_REQ;
bg_state.spinup_ring_streak_achieved = p3_ring_streak;
bg_state.spinup_ring_ungated_pools = P3_UNGATED_POOLS;
bg_state.spinup_observation_constraint = BG_SPINUP_CONSTRAINT;
bg_state.spinup_transient_observation_constraint = BG_TRANSIENT_INITIAL_CONSTRAINT;
bg_state.spinup_transient_observation_constraint_diag = handoff_constraint_diag;
if exist('Bio_Zs', 'var') && exist('rsd', 'var') && exist('Zbio', 'var')
    spinup_biogeochemical_soil_mass_kg_m2 = sum(Bio_Zs.*rsd)*Zbio/1000;
elseif isfield(bg_state, 'spinup_biogeochemical_soil_mass_kg_m2') && ...
        isscalar(bg_state.spinup_biogeochemical_soil_mass_kg_m2) && ...
        isfinite(bg_state.spinup_biogeochemical_soil_mass_kg_m2) && ...
        bg_state.spinup_biogeochemical_soil_mass_kg_m2 > 0
    spinup_biogeochemical_soil_mass_kg_m2 = ...
        bg_state.spinup_biogeochemical_soil_mass_kg_m2;
elseif ismember('rsd', MOD_PARAM_input.Properties.VariableNames) && ...
        isfinite(MOD_PARAM_input.rsd) && isfinite(MOD_PARAM_input.Zbio) && ...
        MOD_PARAM_input.rsd > 0 && MOD_PARAM_input.Zbio > 0
    spinup_biogeochemical_soil_mass_kg_m2 = ...
        MOD_PARAM_input.rsd * MOD_PARAM_input.Zbio / 1000;
else
    error('Run_BG_Spinup:InvalidSoilMass', ...
        'The biogeochemical layer must have a finite positive soil mass.');
end
bg_state.spinup_biogeochemical_soil_mass_kg_m2 = ...
    spinup_biogeochemical_soil_mass_kg_m2;
% Read the P budget off the nutrient-budget struct, not the loop-local
% p3_p_budget: the ring loop's clearvars keeps p3_nutrient_budgets (it is in
% p3_keep) and drops p3_p_budget, and a converged resume never runs the loop
% body at all.  Both hold the same struct (the loop sets .P = p3_p_budget).
bg_state.spinup_phosphorus_budget = p3_nutrient_budgets.P;
bg_state.spinup_nutrient_budget_vers = P3_NUTRIENT_BUDGET_VERS;
bg_state.spinup_nutrient_budgets = p3_nutrient_budgets;
prev_total_P_warn = warning('off', 'TC_SMEW_N_run:TotalPExceedsScreen');
restore_total_P_warn = onCleanup(@() warning(prev_total_P_warn));
Validate_BG_State_For_Transient(bg_state, bg_file);
clear restore_total_P_warn
bg_state.spinup_provenance = spinup_provenance;
%%% A warm start's cycle count is its own; the earlier run's provenance is
%%% part of this one's (the seed label).  Recorded only on a warm start, so
%%% every other state is as before.
if spinup_warm_start
    bg_state.spinup_warm_start = true;
end
if lateral_replay_spinup
    bg_state.spinup_lateral_record_identity = char(string(LATERAL_REPLAY.identity));
    bg_state.reconstruction_method = "replay";
end
if spinup_handoff_only
    for inherited = {'reconstruction_method','reconstruction_template_cell_id', ...
            'reconstruction_rf_cv_rel_error','reconstruction_tracking_key'}
        if isfield(p3_handoff_source,inherited{1})
            bg_state.(inherited{1}) = p3_handoff_source.(inherited{1});
        end
    end
end
%%% The instant this state hands off at: the first hour of the transient
%%% it seeds.  Load_Field_Spinup_States refuses it for a run that starts
%%% at any other instant.
bg_state.spinup_handoff_dn = transient_start;
%%% A cell's identity is stamped here, before the one atomic write, so no
%%% reader ever sees this cell's state without it.
if exist('RUN_SETTINGS','var') && isstruct(RUN_SETTINGS) && ...
        isfield(RUN_SETTINGS,'CELL_ID')
    bg_state.cell_id = double(RUN_SETTINGS.CELL_ID);
    bg_state.cell_agronomic_identity = ...
        char(string(RUN_SETTINGS.CELL_AGRONOMIC_IDENTITY));
end
save_bg_state(bg_file, bg_state);
fprintf('Spin-up complete. Final state: %s\n', bg_file);

%% =========================================================================
%% Local helper functions
%% =========================================================================
function [B, diag] = spinup_apply_observation_constraint(B, constraint)
B_before = B;
diag = struct('applied', false, 'source', '', 'applied_fields', {{}}, ...
    'delta_by_pool', zeros(1,79));
if isempty(constraint)
    return
end
if ~isstruct(constraint) || ~isscalar(constraint)
    error('Run_BG_Spinup:InvalidObservationConstraint', ...
        'A spin-up observation constraint must be a scalar struct.');
end
if ~isfield(constraint, 'targets')
    return
end
if ~isstruct(constraint.targets) || ~isscalar(constraint.targets)
    error('Run_BG_Spinup:InvalidObservationConstraint', ...
        'A spin-up observation constraint targets field must be a scalar struct.');
end
if isempty(fieldnames(constraint.targets))
    return
end
options = struct();
if isfield(constraint, 'options')
    options = constraint.options;
end
[B, map_diag] = Initialize_BG_State_From_Observations(B, constraint.targets, options);
diag = map_diag;
diag.applied = true;
diag.source = '';
if isfield(constraint, 'source')
    diag.source = char(string(constraint.source));
end
diag.delta_by_pool = B - B_before;
diag.delta_C_gCm2 = sum(diag.delta_by_pool(1:22));
diag.delta_N_gNm2 = sum(diag.delta_by_pool(23:34));
diag.delta_P_gPm2 = sum(diag.delta_by_pool(35:47));
diag.delta_K_gKm2 = sum(diag.delta_by_pool(48:55));
diag.delta_Ca_gCam2 = sum(diag.delta_by_pool(56:63));
diag.delta_Mg_gMgm2 = sum(diag.delta_by_pool(64:71));
diag.delta_Si_gSim2 = sum(diag.delta_by_pool(72:79));
end

function scale = spinup_pool_ring_scale(B_by_ring, pool_scale_floor)
%%% Ring-mean level of each pool, used to non-dimensionalise Phase 3 drift.
%%% The scale must not vary across the ring. Dividing by a pool's level at a
%%% single slot holds a pool that swings by a factor S to REL_TOL/S at its
%%% leanest slot, and the all-slots gate then binds on that slot. Mean rather
%%% than max: tighter, and not swayed by one extreme slot.
populated = B_by_ring(~cellfun(@isempty, B_by_ring));
if isempty(populated)
    error('Run_BG_Spinup:EmptyRingScale', ...
        'A pool ring scale requires at least one populated ring position.');
end
ring = cell2mat(cellfun(@(b) reshape(b(1:79), 1, 79), populated(:), ...
    'UniformOutput', false));
scale = max(mean(abs(ring), 1), pool_scale_floor);
end

function diag = spinup_pool_group_diagnostics(B_now, B_prev, cycle_years, pool_scale_floor, pool_ring_scale, ungated)
%%% Pool 46 has no return flux and is excluded from the ring gate. Total P
%%% drift remains diagnostic; ERW reads mineral P in pool 43 only.  Pools in
%%% ungated (Phase 3: P3_UNGATED_POOLS) leave every group, pool list and
%%% element total, so their drift neither passes nor fails a ring; Phase 2
%%% passes none.
if nargin < 6; ungated = []; end
dynamic_idx = spinup_drop_ungated_pools( ...
    {1:22, 23:34, [35:43 47], 48:53, 56:61, 64:69, 72:77}, ungated);
slow_idx = spinup_drop_ungated_pools({44:45, 54:55, 62:63, 70:71, 78:79}, ungated);
B_now = B_now(:).';
B_prev = B_prev(:).';

if numel(B_now) < 79 || numel(B_prev) < 79
    error('Run_BG_Spinup:BadPoolVector', ...
        'Spin-up convergence requires all 79 biogeochemical pools.');
end
if ~(isscalar(cycle_years) && isfinite(cycle_years) && cycle_years >= 1)
    error('Run_BG_Spinup:BadCycleLength', ...
        'The meteorological cycle length must be at least one year.');
end
if ~(isnumeric(pool_scale_floor) && isreal(pool_scale_floor) && ...
        isscalar(pool_scale_floor) && isfinite(pool_scale_floor) && pool_scale_floor > 0)
    error('Run_BG_Spinup:BadPoolScaleFloor', ...
        'The convergence pool-scale floor must be finite and positive.');
end
if ~(isnumeric(pool_ring_scale) && isreal(pool_ring_scale) && ...
        numel(pool_ring_scale) >= 79 && all(isfinite(pool_ring_scale(1:79))))
    error('Run_BG_Spinup:BadPoolRingScale', ...
        'The convergence pool ring scale must be finite and real for all 79 pools.');
end
%%% Floor each pool on its own mass and each group on the group's mass: summing
%%% floored members would inflate a group denominator by its pool count.
raw_pool_scale = abs(pool_ring_scale(:).');
pool_ring_scale = max(raw_pool_scale, pool_scale_floor);
group_scale = @(idx) max(sum(raw_pool_scale(idx)), pool_scale_floor);

diag.pool_basis_idx = 1:79;
diag.pool_scale = pool_ring_scale(diag.pool_basis_idx);
diag.pool_ring_abs = abs(B_now(diag.pool_basis_idx) - B_prev(diag.pool_basis_idx));
diag.pool_ring_rel = diag.pool_ring_abs ./ diag.pool_scale;
%%% Pools are non-negative, so the summed scale is the group ring-mean.
diag.dynamic_group_ring_abs = zeros(1, numel(dynamic_idx));
diag.dynamic_group_ring_rel = zeros(1, numel(dynamic_idx));
for i = 1:numel(dynamic_idx)
    idx = dynamic_idx{i};
    diag.dynamic_group_ring_abs(i) = abs(sum(B_now(idx)) - sum(B_prev(idx)));
    diag.dynamic_group_ring_rel(i) = diag.dynamic_group_ring_abs(i)/group_scale(idx);
end
dynamic_pool_idx = [dynamic_idx{:}];
diag.dynamic_pool_ring_abs = diag.pool_ring_abs(dynamic_pool_idx);
diag.dynamic_pool_ring_rel = diag.pool_ring_rel(dynamic_pool_idx);
diag.dynamic_ring_abs = [diag.dynamic_group_ring_abs diag.dynamic_pool_ring_abs];
diag.dynamic_ring_rel = [diag.dynamic_group_ring_rel diag.dynamic_pool_ring_rel];
diag.dynamic_group_idx = dynamic_idx;
diag.dynamic_pool_idx = dynamic_pool_idx;

diag.slow_group_ring_abs = zeros(1, numel(slow_idx));
diag.slow_group_ring_rel = zeros(1, numel(slow_idx));
for i = 1:numel(slow_idx)
    idx = slow_idx{i};
    diag.slow_group_ring_abs(i) = abs(sum(B_now(idx)) - sum(B_prev(idx)));
    diag.slow_group_ring_rel(i) = diag.slow_group_ring_abs(i)/group_scale(idx);
end
slow_pool_idx = [slow_idx{:}];
diag.slow_pool_ring_abs = diag.pool_ring_abs(slow_pool_idx);
diag.slow_pool_ring_rel = diag.pool_ring_rel(slow_pool_idx);
diag.slow_ring_abs = [diag.slow_group_ring_abs diag.slow_pool_ring_abs];
diag.slow_ring_rel = [diag.slow_group_ring_rel diag.slow_pool_ring_rel];
diag.slow_group_annual_rel = diag.slow_group_ring_rel / cycle_years;
diag.slow_pool_annual_rel = diag.slow_pool_ring_rel / cycle_years;
diag.slow_annual_rel = diag.slow_ring_rel / cycle_years;
diag.slow_group_idx = slow_idx;
diag.slow_pool_idx = slow_pool_idx;

% Gate total N and non-occluded P as well as individual pools.
element_idx = struct('N', 23:34, 'P', setdiff([35:45 47], ungated, 'stable'));
element_names = fieldnames(element_idx);
diag.element_names = element_names;
diag.element_ring_abs = zeros(1, numel(element_names));
diag.element_ring_rel = zeros(1, numel(element_names));
for e = 1:numel(element_names)
    idx = element_idx.(element_names{e});
    diag.element_ring_abs(e) = abs(sum(B_now(idx)) - sum(B_prev(idx)));
    diag.element_ring_rel(e) = diag.element_ring_abs(e)/group_scale(idx);
end

diag.finite_nonnegative = all(isfinite(B_now(1:79))) && ...
    all(isfinite(B_prev(1:79))) && all(B_now(1:79) >= 0) && ...
    all(B_prev(1:79) >= 0);
end

function groups = spinup_drop_ungated_pools(groups, ungated)
%%% The pool groups without the ungated pools; a group left empty goes.
groups = cellfun(@(g) setdiff(g, ungated, 'stable'), groups, 'UniformOutput', false);
groups = groups(~cellfun(@isempty, groups));
end

function label = spinup_worst_offender_label(rel_vec, group_idx, pool_idx)
%%% Identify which pool group (or individual pool) dominates a Phase 3
%%% convergence metric's max, so a stalled ring-pass gate can be diagnosed
%%% from the log alone.  rel_vec = [group_rel, pool_rel], exactly as built by
%%% spinup_pool_group_diagnostics; group_idx/pool_idx are its matching
%%% *_group_idx/*_pool_idx index maps (dynamic or slow).
n_groups = numel(group_idx);
[worst_val, worst_i] = max(rel_vec);
if isempty(worst_i) || ~isfinite(worst_val)
    label = 'n/a';
    return
end
if worst_i <= n_groups
    label = sprintf('group %d total (pools %s)', worst_i, mat2str(group_idx{worst_i}));
    return
end
pool_num = pool_idx(worst_i - n_groups);
group_num = 0;
for g = 1:n_groups
    if any(group_idx{g} == pool_num)
        group_num = g;
        break
    end
end
label = sprintf('pool #%d%s (group %d)', pool_num, spinup_known_pool_name(pool_num), group_num);
end

function name = spinup_known_pool_name(pool_num)
%%% Chemical labels verified against MAIN_FRAME_BG_ERW_10mins.m's P(.,N)
%%% assignment comments.  An unlabeled pool number returns '' rather than a
%%% guessed name.
persistent known
if isempty(known)
    known = containers.Map('KeyType', 'double', 'ValueType', 'char');
    known(31) = 'NH4'; known(32) = 'NO3'; known(43) = 'PO4';
    known(52) = 'K';   known(60) = 'Ca';  known(68) = 'Mg'; known(76) = 'Si';
end
if isKey(known, pool_num)
    name = sprintf(' (%s)', known(pool_num));
else
    name = '';
end
end

function [streak, sustained_pass] = spinup_update_consecutive_passes(single_pass, streak, required_passes)
if single_pass
    streak = streak + 1;
else
    streak = 0;
end
sustained_pass = streak >= required_passes;
end

function [ok, p1_avgs] = load_compatible_p1_checkpoint(filename, cycle_years, expected_version)
ok = false;
p1_avgs = struct();
if ~isfile(filename)
    return
end
try
    checkpoint = load(filename, 'p1_avgs');
catch ME
    fprintf('  Ignoring unreadable Phase 1 checkpoint: %s\n', ME.message);
    return
end
if isfield(checkpoint, 'p1_avgs')
    p1_avgs = checkpoint.p1_avgs;
    ok = isfield(p1_avgs, 'spinup_cycle_years') && ...
        isequal(p1_avgs.spinup_cycle_years, cycle_years) && ...
        isfield(p1_avgs, 'standing_nutrient_state_vers') && ...
        isequal(p1_avgs.standing_nutrient_state_vers, expected_version) && ...
        isfield(p1_avgs, 'f_anaer_bio') && ...
        isscalar(p1_avgs.f_anaer_bio) && isfinite(p1_avgs.f_anaer_bio) && ...
        p1_avgs.f_anaer_bio >= 0 && p1_avgs.f_anaer_bio <= 1 && ...
        isfield(p1_avgs, 'O_out') && ...
        isnumeric(p1_avgs.O_out) && isreal(p1_avgs.O_out) && ...
        isvector(p1_avgs.O_out) && ~isempty(p1_avgs.O_out) && ...
        all(isfinite(p1_avgs.O_out(:))) && ...
        isfield(p1_avgs, 'phase1_complete') && ...
        isscalar(p1_avgs.phase1_complete) && logical(p1_avgs.phase1_complete) && ...
        isfield(p1_avgs, 'phase1_water_state_physical') && ...
        isscalar(p1_avgs.phase1_water_state_physical) && ...
        logical(p1_avgs.phase1_water_state_physical);
end
if ~ok
    p1_avgs = struct();
end
end

function [ok, bg_state] = load_compatible_p2_checkpoint(filename, cycle_years, constraint, expected_version)
[ok, bg_state] = load_compatible_bg_checkpoint(filename, cycle_years, constraint, ...
    'phase2_observation_constraint');
if ~ok
    return
end
ok = isfield(bg_state, 'standing_nutrient_state_vers') && ...
    isequal(bg_state.standing_nutrient_state_vers, expected_version);
if ~ok
    fprintf('  Ignoring Phase 2 checkpoint: standing nutrient state changed.\n');
    bg_state = struct();
    return
end
phase2_converged = isfield(bg_state, 'phase2_converged') && ...
    (islogical(bg_state.phase2_converged) || isnumeric(bg_state.phase2_converged)) && ...
    isscalar(bg_state.phase2_converged) && isfinite(bg_state.phase2_converged) && ...
    logical(bg_state.phase2_converged);
ok = phase2_converged;
if ~ok
    fprintf('  Ignoring incomplete Phase 2 checkpoint.\n');
    bg_state = struct();
end
end

function [ok, bg_state] = load_compatible_bg_checkpoint(filename, cycle_years, constraint, constraint_field)
ok = false;
bg_state = struct();
if ~isfile(filename)
    return
end
try
    bg_state = load_bg_state(filename);
catch ME
    fprintf('  Ignoring unreadable spin-up checkpoint %s: %s\n', filename, ME.message);
    bg_state = struct();
    return
end
ok = isfield(bg_state, 'spinup_cycle_years') && ...
    isequal(bg_state.spinup_cycle_years, cycle_years);
%%% A biogeochemical checkpoint stores pool state the constraint already shaped,
%%% so a changed constraint makes it a different starting point rather than a
%%% resumable one.  This compares the constraint itself, not a version stamp: a
%%% stamp has to be remembered and bumped by hand, whereas the struct is the
%%% thing that actually determined the stored state.
if ok
    ok = isfield(bg_state, constraint_field) && ...
        isequaln(bg_state.(constraint_field), constraint);
    if ~ok
        fprintf(['  Ignoring spin-up checkpoint %s: written under a different ' ...
            'observation constraint.\n'], filename);
    end
end
if ~ok
    bg_state = struct();
end
end

function ok = spinup_provenance_matches(state, provenance)
%%% Without a requested provenance every checkpoint that passes the other
%%% gates is resumed, as before.
ok = isempty(provenance) || ...
    (isstruct(state) && isfield(state,'spinup_provenance') && ...
    strcmp(char(string(state.spinup_provenance)), provenance));
end

function B = Project_Small_Negative_BG_Pools(B, tolerance)
valid = isnumeric(B) && isreal(B) && ~isempty(B) && all(isfinite(B), 'all') && ...
    isscalar(tolerance) && isfinite(tolerance) && tolerance >= 0;
if ~valid
    error('Run_BG_Spinup:InvalidPoolProjection', ...
        'Pool projection requires finite real pools and a non-negative tolerance.');
end

large_negative = find(B < -tolerance);
if ~isempty(large_negative)
    error('Run_BG_Spinup:NegativeBiogeochemicalPool', ...
        ['Biogeochemical pools %s are below the numerical projection limit ' ...
        '(minimum %.6g g element m-2; limit %.6g).'], ...
        mat2str(large_negative(:).'), min(B(large_negative)), tolerance);
end
B(B < 0) = 0;
end

function ok = spinup_nutrient_budget_resumed(state)
ok = false;
if ~isfield(state, 'p3_nutrient_budgets')
    return
end
budgets = state.p3_nutrient_budgets;
if ~isstruct(budgets) || isempty(fieldnames(budgets))
    return
end
elements = {'P','N','K','Ca','Mg','Si'};
for e = 1:numel(elements)
    if ~isfield(budgets, elements{e})
        return
    end
    b = budgets.(elements{e});
    if ~isstruct(b) || ~all(isfield(b, {'closed','residual','tolerance'}))
        return
    end
    valid_closed = (islogical(b.closed) || isnumeric(b.closed)) && ...
        isreal(b.closed) && isscalar(b.closed) && ...
        isfinite(double(b.closed)) && logical(b.closed);
    valid_residual = isnumeric(b.residual) && isreal(b.residual) && ...
        isscalar(b.residual) && isfinite(b.residual);
    valid_tolerance = isnumeric(b.tolerance) && isreal(b.tolerance) && ...
        isscalar(b.tolerance) && isfinite(b.tolerance) && b.tolerance >= 0;
    if ~valid_closed || ~valid_residual || ~valid_tolerance || ...
            abs(b.residual) > b.tolerance
        return
    end
end
ok = true;
end
