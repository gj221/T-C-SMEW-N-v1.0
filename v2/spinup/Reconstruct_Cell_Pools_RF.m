function [CellBGStates, report] = Reconstruct_Cell_Pools_RF( ...
    topo, resolved, CellInputs, subsample_states, subsample_idx, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct_Cell_Pools_RF
%%
%% Lian et al. (2025) random-forest initialisation of the distributed soil
%% biogeochemistry pools, with ERW charge-balance closure through the
%% model's own validated handoff path.
%%
%%   1. A topographically-stratified subsample of cells is spun up in full
%%      (the caller does this and passes the results here).
%%   2. A random forest is trained on {elevation, slope, aspect, curvature,
%%      TWI, soil class, vegetation} -> the conserved biological / nutrient
%%      pool totals of the subsample cells (Lian et al. Fig. 3e).
%%   3. Every other cell gets its totals predicted, mapped onto the 79-pool
%%      state with Initialize_BG_State_From_Observations (stoichiometry
%%      preserved), seeded from its nearest subsample neighbour's full
%%      state for everything the RF does not predict (exchange complex,
%%      pore-water pH, plant reserves, moisture).
%%   4. That seed is handed to the standard Run_Automatic_Site_Spinup with
%%      the standard convergence criteria: starting near equilibrium it
%%      needs fewer rings than a cold start, and it runs the ERW handoff
%%      gate (ERW_Handoff_Readiness),
%%      so electroneutrality and pH are closed by the same code the
%%      transient run trusts.  Any cell that fails the gate or does not
%%      converge falls back to a full per-cell spin-up.
%%
%% The RF captures the topo-climatic initial condition (elevation ->
%% temperature dominates, Lian et al.); the lateral redistribution signal
%% then develops during the 2-D run itself (total-head Darcy + overland
%% routing).  A flux-tracking spin-up (Lian et al. Fig. 3c-d) would fold
%% that signal into the initial condition too -- future work.
%%
%% INPUTS
%%   topo              topography-attached topology
%%   resolved          Resolve_Cell_Agronomic_Inputs table (per cell)
%%   CellInputs        Configure_Cell_Agronomic_Inputs struct array
%%   subsample_states  struct array of full bg_state for the subsample cells
%%   subsample_idx     indices (into 1..topo.N_cells) of those cells
%%
%% Name/value
%%   'CaseNumber'       workbook scenario (default 5)
%%   'OutputDirectory'  where per-cell seed / state files go
%%   'RFOptions'        cell of name/values passed to Train_Pool_Random_Forest
%%   'CellIndex'        reconstruct ONLY this cell (index into 1..N_cells)
%%                      and return a one-element CellBGStates: the
%%                      cell-parallel mode a SLURM array uses, one cell
%%                      per task, with the forest retrained identically
%%                      (same subsample, same seed) in every task
%%
%% OUTPUTS
%%   CellBGStates  [N_cells x 1] validated bg_state struct array
%%   report        struct: .rf (the trained model), .method_per_cell
%%                 ('full' | 'rf_seed' | 'fallback' | 'failed'),
%%                 .gate_pass_rate, .cv_rel_error, .cells (the indices
%%                 reconstructed), .subsample_idx, .nn_of
%%   Every returned state carries reconstruction_method (string) and
%%   reconstruction_failed (logical), so the struct array stays uniform
%%   and an assembly step can report how each cell was obtained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spinup_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(spinup_dir);
code = fullfile(v2,'T&C_Code');
addpath(spinup_dir, code, fullfile(code,'Parameters'));

opt = struct('CaseNumber',5,'OutputDirectory',tempname,'RFOptions',{{}}, ...
    'CellIndex',[]);
for k = 1:2:numel(varargin)
    if ~isfield(opt,char(varargin{k}))
        error('Reconstruct_Cell_Pools_RF:UnknownOption','Unknown option %s.',char(varargin{k}));
    end
    opt.(char(varargin{k})) = varargin{k+1};
end
if ~exist(opt.OutputDirectory,'dir'); mkdir(opt.OutputDirectory); end

N = topo.N_cells;
subsample_idx = subsample_idx(:);
assert(numel(subsample_states) == numel(subsample_idx), ...
    'One state per subsample cell.');

%% --- features for every cell --------------------------------------
[X_all, feature_names] = cell_features(topo, resolved);
X_sub = X_all(subsample_idx,:);

%% --- conserved-total targets from the subsample states -----------
target_names = {'POC_gCm2','MAOC_gCm2','DOC_gCm2','MBC_gCm2', ...
    'total_N_gNm2','NH4_gNm2','NO3_gNm2','total_P_gPm2'};
Y_sub = zeros(numel(subsample_idx), numel(target_names));
for s = 1:numel(subsample_idx)
    Y_sub(s,:) = pool_totals(subsample_states(s).B_out(end,:), target_names);
end

%% --- train + predict --------------------------------------------
rf = Train_Pool_Random_Forest(X_sub, Y_sub, feature_names, target_names, ...
    opt.RFOptions{:});
Y_all = rf.predict(X_all);

%% --- reconstruct + close every cell -----------------------------
cells_to_do = 1:N;
if ~isempty(opt.CellIndex)
    if ~(isnumeric(opt.CellIndex) && isscalar(opt.CellIndex) && ...
            opt.CellIndex == fix(opt.CellIndex) && opt.CellIndex >= 1 && opt.CellIndex <= N)
        error('Reconstruct_Cell_Pools_RF:BadCellIndex', ...
            'CellIndex must be an integer in 1..%d.', N);
    end
    cells_to_do = opt.CellIndex;
end
CellBGStates = repmat(stamp(subsample_states(1),"full",false), N, 1);
method_per_cell = strings(N,1);
gate_pass = false(N,1);
%% Nearest subsample neighbour in STANDARDISED feature space: the raw
%% features mix metres (elevation), fractions (slope) and log units (TWI),
%% so an unscaled Euclidean distance is just an elevation match.  The
%% scale comes from the WHOLE field, not the subsample, so a covariate
%% the subsample happens not to span still counts at its field-wide
%% scale; a covariate constant across the field contributes nothing
%% (dividing by eps instead would let a 1e-16 spread dominate the metric
%% by thirty orders of magnitude and annihilate every other feature).
mu_f = mean(X_all,1); sd_f = std(X_all,0,1);
flat = sd_f <= 1e-12*max(1,max(abs(X_all),[],1));
sd_f(flat) = 1;
Xs = (X_sub-mu_f)./sd_f; Xa = (X_all-mu_f)./sd_f;
Xs(:,flat) = 0; Xa(:,flat) = 0;
[~,nn_of] = pdist2(Xs, Xa, 'euclidean', 'Smallest', 1);

for c = cells_to_do
    is_sub = find(subsample_idx == c, 1);
    if ~isempty(is_sub)
        CellBGStates(c) = stamp(subsample_states(is_sub),"full",false);
        method_per_cell(c) = "full";
        gate_pass(c) = handoff_ready(subsample_states(is_sub));
        continue
    end
    template = subsample_states(nn_of(c));
    targets = totals_struct(Y_all(c,:), target_names);
    %% Carry the template's OWN phosphorus partition into the mapping:
    %% its SOM C:P sets organic P from the (rescaled) SOM carbon, so
    %% organic P scales with the RF's carbon prediction and the SOM C:P
    %% guard stays live; its inorganic fractions keep the converged
    %% labile/adsorbed/occluded split instead of the cold-start defaults.
    %% preserve_organic_P_pools would instead freeze B(38) at the
    %% template's absolute value while SOC moves, drifting SOM C:P and
    %% switching the guard off in exactly the case that needs it.  With
    %% the template's own totals this mapping is an identity.
    Bt = template.B_out(end,:);
    cell_file = fullfile(opt.OutputDirectory, sprintf('cell_%04d_BG.mat', ...
        resolved.cell_id(c)));
    seed_file = fullfile(opt.OutputDirectory, sprintf('rf_seed_%04d.mat', ...
        resolved.cell_id(c)));
    try
        %% The mapping can refuse a prediction (inconsistent nitrogen
        %% stocks, implausible stoichiometry); that is a reason for the
        %% cold-start fallback below, not for abandoning the field.
        p_opts = struct();
        if Bt(38) > 0 && sum(Bt(9:11)) > 0
            p_opts.organic_P_C_to_P = sum(Bt(9:11))/Bt(38);
        end
        if sum(Bt(43:46)) > 0
            p_opts.inorganic_P_fractions = Bt(43:46)/sum(Bt(43:46));
        end
        B_seed = Initialize_BG_State_From_Observations(Bt, targets, p_opts);
        seed_state = template;
        seed_state.B_out = B_seed;
        save_seed(seed_file, seed_state);
        %% Phase checkpoints left by an earlier attempt for this cell would
        %% be resumed regardless of the seed, so a re-run would carry a cold
        %% start's state under the rf_seed label; start the seeded run clean.
        clear_spinup_checkpoints(cell_file);
        settings = struct('NN_case_study_ERW_MOD',opt.CaseNumber, ...
            'RUN_MODE','spinup','RUN_CHEMISTRY','erw','RUN_DOMAIN','single', ...
            'SAVE_RESULTS',0,'MOD_PARAM_override',CellInputs(c).MOD_PARAM, ...
            'ERW_PARAM_override',CellInputs(c).ERW_PARAM, ...
            'BG_OUTPUT_FILE',cell_file,'CELL_ID',resolved.cell_id(c), ...
            'CELL_AGRONOMIC_IDENTITY',CellInputs(c).agronomic_identity, ...
            'SPINUP_SEED_FILE',seed_file);
        Run_Automatic_Site_Spinup(settings);
        st = load_bg_state(cell_file);
        Validate_BG_State_For_Transient(st, cell_file);
        CellBGStates(c) = stamp(st,"rf_seed",false);
        method_per_cell(c) = "rf_seed";
        gate_pass(c) = handoff_ready(st);
    catch ME
        warning('Reconstruct_Cell_Pools_RF:CellFallback', ...
            'Cell %d RF seed failed (%s); cold-start spin-up fallback.', ...
            resolved.cell_id(c), ME.message);
        %% A real cold start.  Run_BG_Spinup derives its phase
        %% checkpoints from BG_OUTPUT_FILE and resumes them on whatever
        %% state they hold, with no record of whether it came from a
        %% seed -- so simply dropping SPINUP_SEED_FILE would resume the
        %% failed seeded state and fail again identically.
        clear_spinup_checkpoints(cell_file);
        settings = rmfield(settings, 'SPINUP_SEED_FILE');
        try
            Run_Automatic_Site_Spinup(settings);
            st = load_bg_state(cell_file);
            Validate_BG_State_For_Transient(st, cell_file);
            CellBGStates(c) = stamp(st,"fallback",false);
            method_per_cell(c) = "fallback";
            gate_pass(c) = handoff_ready(st);
        catch ME2
            %% One unreconstructable cell must not discard the whole
            %% field: keep its nearest-neighbour template, flag it, and
            %% let the caller decide (report.n_failed).
            warning('Reconstruct_Cell_Pools_RF:CellFailed', ...
                ['Cell %d cold-start spin-up also failed (%s); keeping its ' ...
                 'nearest subsample neighbour state and flagging the cell.'], ...
                resolved.cell_id(c), ME2.message);
            %% Stamp the cell's OWN identity and clear the converged /
            %% handoff flags: a verbatim copy would carry the neighbour's
            %% cell_id (breaking the one-state-per-cell invariant that
            %% Load_Field_Spinup_States enforces) and would look converged
            %% and handoff-ready when it is neither.
            failed_state = template;
            failed_state.cell_id = resolved.cell_id(c);
            failed_state.cell_agronomic_identity = CellInputs(c).agronomic_identity;
            if isfield(failed_state,'spinup_converged')
                failed_state.spinup_converged = false;
            end
            if isfield(failed_state,'spinup_erw_handoff_ready')
                failed_state.spinup_erw_handoff_ready = false;
            end
            CellBGStates(c) = stamp(failed_state,"failed",true);
            method_per_cell(c) = "failed";
            gate_pass(c) = false;
        end
    end
end

if ~isempty(opt.CellIndex)
    CellBGStates = CellBGStates(cells_to_do);
    method_per_cell = method_per_cell(cells_to_do);
    gate_pass = gate_pass(cells_to_do);
end
report = struct();
report.rf = rf;
report.cells = cells_to_do(:);
report.subsample_idx = subsample_idx;
report.nn_of = nn_of(:);
report.method_per_cell = method_per_cell;
report.gate_pass_rate = mean(gate_pass);
report.cv_rel_error = rf.cv_rel_error;
report.n_full = nnz(method_per_cell == "full");
report.n_rf_seed = nnz(method_per_cell == "rf_seed");
report.n_fallback = nnz(method_per_cell == "fallback");
report.n_failed = nnz(method_per_cell == "failed");
fprintf(['Reconstruct_Cell_Pools_RF: %d full, %d RF-seeded, %d fallback, ' ...
    '%d FAILED; handoff gate pass %.1f%%\n'], report.n_full, report.n_rf_seed, ...
    report.n_fallback, report.n_failed, 100*report.gate_pass_rate);
if report.n_failed > 0
    warning('Reconstruct_Cell_Pools_RF:UnreconstructedCells', ...
        ['%d cell(s) could not be spun up and hold a neighbour''s state; ' ...
         'they are listed in report.method_per_cell == "failed".'], ...
        report.n_failed);
end
end

%% =======================================================================
function clear_spinup_checkpoints(cell_file)
%% Remove the phase checkpoints Run_BG_Spinup derives from this output
%% path, so the fallback genuinely cold-starts.
[d,n] = fileparts(cell_file);
pattern = fullfile(d,[n '_checkpoint_*.mat']);
listing = dir(pattern);
for k = 1:numel(listing)
    f = fullfile(listing(k).folder,listing(k).name);
    try
        delete(f);
    catch
        warning('Reconstruct_Cell_Pools_RF:CheckpointNotCleared', ...
            'Could not delete %s; the fallback may resume a seeded state.',f);
    end
end
end

%% =======================================================================
function [X, names] = cell_features(topo, resolved)
%% One feature definition shared with the validation experiment.
[X, names] = Cell_Topographic_Features(topo, resolved);
end

%% =======================================================================
function t = pool_totals(B, names)
%% Conserved element totals from a 79-pool state (g element m-2); the
%% one definition shared with the validation experiment.
t = Conserved_Pool_Totals(B, names);
end

%% =======================================================================
function st = stamp(st, method, failed)
%% Every state carries how it was obtained, so struct arrays stay
%% uniform whichever path produced each element.
st.reconstruction_method = string(method);
st.reconstruction_failed = logical(failed);
end

%% =======================================================================
function s = totals_struct(row, names)
s = struct();
for k = 1:numel(names)
    if isfinite(row(k)) && row(k) > 0
        s.(names{k}) = row(k);
    end
end
end

%% =======================================================================
function tf = handoff_ready(st)
tf = isfield(st,'spinup_erw_handoff_ready') && ...
    logical(st.spinup_erw_handoff_ready);
end

%% =======================================================================
function save_seed(file, state) %#ok<INUSD>
S = state; %#ok<NASGU>
save(file, '-struct', 'S', '-v7.3');
end
