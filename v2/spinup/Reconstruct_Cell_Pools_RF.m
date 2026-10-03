function [CellBGStates, report] = Reconstruct_Cell_Pools_RF( ...
    topo, resolved, CellInputs, replayed_states, tracked_idx, handoff_states, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct_Cell_Pools_RF
%%
%% The last stage of the flux-tracking field spin-up (Lian et al. 2025,
%% Fig. 3e-f), with the ERW handoff closed by the model's own seed path.
%%
%%   1. The tracked cells have reached their dynamic steady state under
%%      their own recorded lateral exchange (the replays; the caller passes
%%      their states).
%%   2. A random forest is trained on {elevation, slope, aspect, curvature,
%%      TWI, soil} -> the conserved biological / nutrient pool totals of the
%%      replayed cells (Conserved_Pool_Totals), with k-fold cross-validation
%%      (Lian et al.: < 10 % relative error).
%%   3. Every other cell's totals are predicted and mapped onto the 79-pool
%%      state of its nearest replayed cell in standardised covariate space
%%      (Initialize_BG_State_From_Observations, the template's own
%%      partition and stoichiometry kept); its soil moisture is its own,
%%      from the 2-D tracking run at the transient start.
%%   4. That state is the cell's initial condition (as in Lian et al., no
%%      further spin-up), but -- the ERW adaptation -- it is handed on only
%%      through the transient's own seed path: Run_Automatic_Site_Spinup
%%      with SPINUP_HANDOFF_ONLY applies the transient-start observation
%%      constraint and builds the ERW exchanger and pore-water seed from the
%%      cell's own bulk pH, CEC and measured base saturation, gated by
%%      ERW_Handoff_Readiness (charge balance, pH, base saturation).  The
%%      forest never predicts the exchanger, pore-water pH, Al, DIC or the
%%      residual anion, so electroneutrality holds by construction.
%%   A cell whose reconstruction fails falls back to its own cold
%%   plot-scale spin-up (no lateral exchange), flagged 'fallback'; one whose
%%   fallback fails too is returned flagged 'failed'.
%%
%% INPUTS
%%   topo             topography-attached topology
%%   resolved         Resolve_Cell_Agronomic_Inputs table (per cell)
%%   CellInputs       Configure_Cell_Agronomic_Inputs struct array
%%   replayed_states  struct array of the tracked cells' replayed states
%%   tracked_idx      indices (into 1..topo.N_cells) of those cells
%%   handoff_states   [N_cells x 1] struct of every cell's 2-D state at the
%%                    transient start (Run_Field_Lateral_Tracking), or []
%%                    to hand on the template's soil moisture
%%
%% Name/value
%%   'CaseNumber'       workbook scenario (default 5)
%%   'OutputDirectory'  where the per-cell seed and state files go
%%   'RFOptions'        cell of name/values for Train_Pool_Random_Forest
%%   'CellIndex'        reconstruct ONLY this cell (index into 1..N_cells)
%%                      and return a one-element CellBGStates: the
%%                      cell-parallel mode of a SLURM array, every task
%%                      retraining the identical forest
%%
%% OUTPUTS
%%   CellBGStates  [N_cells x 1] validated states (or one, with CellIndex)
%%   report        .rf (the trained forest), .method_per_cell ('replay' |
%%                 'rf' | 'fallback' | 'failed'), .gate_pass_rate,
%%                 .cv_rel_error, .cells, .tracked_idx, .nn_of
%%   Every state carries reconstruction_method and reconstruction_failed.
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
tracked_idx = tracked_idx(:);
assert(numel(replayed_states) == numel(tracked_idx), ...
    'One replayed state per tracked cell.');
if ~isempty(handoff_states) && numel(handoff_states) ~= N
    error('Reconstruct_Cell_Pools_RF:BadHandoffStates', ...
        'handoff_states must hold one state per cell (%d).',N);
end

%% --- features for every cell --------------------------------------
[X_all, feature_names] = Cell_Topographic_Features(topo, resolved);
X_sub = X_all(tracked_idx,:);

%% --- conserved-total targets from the replayed states ------------
target_names = {'POC_gCm2','MAOC_gCm2','DOC_gCm2','MBC_gCm2', ...
    'total_N_gNm2','NH4_gNm2','NO3_gNm2','total_P_gPm2'};
Y_sub = zeros(numel(tracked_idx), numel(target_names));
for s = 1:numel(tracked_idx)
    Y_sub(s,:) = Conserved_Pool_Totals(replayed_states(s).B_out(end,:), target_names);
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
%% States from different paths differ in optional fields, so they are
%% gathered as cells and made one struct array at the end.
states = cell(N,1);
method_per_cell = strings(N,1);
gate_pass = false(N,1);
%% Nearest replayed cell in STANDARDISED feature space: the raw features
%% mix metres (elevation), fractions (slope) and log units (TWI), so an
%% unscaled Euclidean distance is just an elevation match.  The scale comes
%% from the WHOLE field, not the subsample, so a covariate the subsample
%% happens not to span still counts at its field-wide scale; a covariate
%% constant across the field contributes nothing (dividing by eps instead
%% would let a 1e-16 spread dominate the metric by thirty orders of
%% magnitude and annihilate every other feature).
mu_f = mean(X_all,1); sd_f = std(X_all,0,1);
flat = sd_f <= 1e-12*max(1,max(abs(X_all),[],1));
sd_f(flat) = 1;
Xs = (X_sub-mu_f)./sd_f; Xa = (X_all-mu_f)./sd_f;
Xs(:,flat) = 0; Xa(:,flat) = 0;
[~,nn_of] = pdist2(Xs, Xa, 'euclidean', 'Smallest', 1);

for c = cells_to_do
    is_tracked = find(tracked_idx == c, 1);
    if ~isempty(is_tracked)
        states{c} = stamp(replayed_states(is_tracked),"replay",false);
        method_per_cell(c) = "replay";
        gate_pass(c) = handoff_ready(replayed_states(is_tracked));
        continue
    end
    template = replayed_states(nn_of(c));
    targets = totals_struct(Y_all(c,:), target_names);
    cell_id = resolved.cell_id(c);
    cell_file = fullfile(opt.OutputDirectory, sprintf('rf_%04d_BG.mat',cell_id));
    seed_file = fullfile(opt.OutputDirectory, sprintf('rf_seed_%04d.mat',cell_id));
    cold_file = fullfile(opt.OutputDirectory, sprintf('cell_%04d_BG.mat',cell_id));
    base = struct('NN_case_study_ERW_MOD',opt.CaseNumber, ...
        'RUN_MODE','spinup','RUN_CHEMISTRY','erw','RUN_DOMAIN','single', ...
        'SAVE_RESULTS',0,'MOD_PARAM_override',CellInputs(c).MOD_PARAM, ...
        'ERW_PARAM_override',CellInputs(c).ERW_PARAM, ...
        'CELL_ID',cell_id,'CELL_AGRONOMIC_IDENTITY',CellInputs(c).agronomic_identity);
    try
        %% Carry the template's OWN phosphorus partition into the mapping:
        %% its SOM C:P sets organic P from the (rescaled) SOM carbon, so
        %% organic P scales with the forest's carbon prediction and the
        %% SOM C:P guard stays live; its inorganic fractions keep the
        %% converged labile/adsorbed/occluded split.  A converged primary-P
        %% fraction of order 1e-14 or a SOM C:P below 103 is the model's
        %% state, not an implausible prior, so the bounds are the
        %% template's own.  The mapping can refuse a prediction
        %% (inconsistent nitrogen stocks); that is a reason for the
        %% fallback below, not for abandoning the field.
        Bt = template.B_out(end,:);
        p_opts = struct();
        if Bt(38) > 0 && sum(Bt(9:11)) > 0
            p_opts.organic_P_C_to_P = sum(Bt(9:11))/Bt(38);
            p_opts.som_C_to_P_bounds = p_opts.organic_P_C_to_P*[1-1e-6 1+1e-6];
        end
        if sum(Bt(43:46)) > 0
            p_opts.inorganic_P_fractions = Bt(43:46)/sum(Bt(43:46));
            p_opts.inorganic_P_fraction_bounds = [0 1; 0 1];
        end
        B_seed = Initialize_BG_State_From_Observations(Bt, targets, p_opts);
        seed_state = template;
        seed_state.B_out = B_seed;
        if ~isempty(handoff_states)
            seed_state.O_out = handoff_states(c).O_out;
        end
        seed_state.cell_id = cell_id;
        seed_state.cell_agronomic_identity = CellInputs(c).agronomic_identity;
        seed_state.reconstruction_method = "rf";
        seed_state.reconstruction_template_cell_id = resolved.cell_id(tracked_idx(nn_of(c)));
        seed_state.reconstruction_rf_cv_rel_error = rf.cv_rel_error;
        save_seed(seed_file, seed_state);
        settings = base;
        settings.BG_OUTPUT_FILE = cell_file;
        settings.SPINUP_SEED_FILE = seed_file;
        settings.SPINUP_HANDOFF_ONLY = true;
        settings.SPINUP_PROVENANCE = Spinup_Provenance(cell_id, ...
            CellInputs(c).agronomic_identity,opt.CaseNumber, ...
            CellInputs(c).MOD_PARAM,CellInputs(c).ERW_PARAM, ...
            ['rf:' mat2str([B_seed(:); seed_state.O_out(:)].',17) ...
             ';template:' sprintf('%d',seed_state.reconstruction_template_cell_id)]);
        st = reuse_or_run(cell_file, settings);
        states{c} = stamp(st,"rf",false);
        method_per_cell(c) = "rf";
        gate_pass(c) = handoff_ready(st);
    catch ME
        warning('Reconstruct_Cell_Pools_RF:CellFallback', ...
            'Cell %d reconstruction failed (%s); cold-start spin-up fallback.', ...
            cell_id, ME.message);
        %% A real cold start: its provenance differs from the
        %% reconstruction's, so neither resumes the other's checkpoints.
        try
            settings = base;
            settings.BG_OUTPUT_FILE = cold_file;
            settings.SPINUP_PROVENANCE = Spinup_Provenance(cell_id, ...
                CellInputs(c).agronomic_identity,opt.CaseNumber, ...
                CellInputs(c).MOD_PARAM,CellInputs(c).ERW_PARAM,'cold');
            st = reuse_or_run(cold_file, settings);
            states{c} = stamp(st,"fallback",false);
            method_per_cell(c) = "fallback";
            gate_pass(c) = handoff_ready(st);
        catch ME2
            %% One unreconstructable cell must not discard the whole
            %% field: keep its nearest replayed state, flag it, and let
            %% the caller decide (report.n_failed).
            warning('Reconstruct_Cell_Pools_RF:CellFailed', ...
                ['Cell %d cold-start spin-up also failed (%s); keeping its ' ...
                 'nearest replayed state and flagging the cell.'], ...
                cell_id, ME2.message);
            %% Stamp the cell's OWN identity and clear the converged /
            %% handoff flags: a verbatim copy would carry the neighbour's
            %% cell_id (breaking the one-state-per-cell invariant that
            %% Load_Field_Spinup_States enforces) and would look converged
            %% and handoff-ready when it is neither.
            failed_state = template;
            failed_state.cell_id = cell_id;
            failed_state.cell_agronomic_identity = CellInputs(c).agronomic_identity;
            if isfield(failed_state,'spinup_converged')
                failed_state.spinup_converged = false;
            end
            if isfield(failed_state,'spinup_erw_handoff_ready')
                failed_state.spinup_erw_handoff_ready = false;
            end
            states{c} = stamp(failed_state,"failed",true);
            method_per_cell(c) = "failed";
            gate_pass(c) = false;
        end
    end
end

CellBGStates = harmonise_states(states(cells_to_do));
if ~isempty(opt.CellIndex)
    method_per_cell = method_per_cell(cells_to_do);
    gate_pass = gate_pass(cells_to_do);
end
report = struct();
report.rf = rf;
report.cells = cells_to_do(:);
report.tracked_idx = tracked_idx;
report.nn_of = nn_of(:);
report.method_per_cell = method_per_cell;
report.gate_pass_rate = mean(gate_pass);
report.cv_rel_error = rf.cv_rel_error;
report.n_replay = nnz(method_per_cell == "replay");
report.n_rf = nnz(method_per_cell == "rf");
report.n_fallback = nnz(method_per_cell == "fallback");
report.n_failed = nnz(method_per_cell == "failed");
fprintf(['Reconstruct_Cell_Pools_RF: %d replayed, %d from the forest, %d fallback, ' ...
    '%d FAILED; handoff gate pass %.1f%%\n'], report.n_replay, report.n_rf, ...
    report.n_fallback, report.n_failed, 100*report.gate_pass_rate);
if report.n_failed > 0
    warning('Reconstruct_Cell_Pools_RF:UnreconstructedCells', ...
        ['%d cell(s) could not be reconstructed and hold a neighbour''s state; ' ...
         'they are listed in report.method_per_cell == "failed".'], ...
        report.n_failed);
end
end

%% =======================================================================
function st = reuse_or_run(cell_file, settings)
%% A finished state of this very run (same provenance key) is reused, so an
%% interrupted array task resumes; anything else is run afresh.
if isfile(cell_file)
    try
        st = load_bg_state(cell_file);
        Validate_BG_State_For_Transient(st, cell_file);
        if isfield(st,'spinup_provenance') && ...
                strcmp(char(string(st.spinup_provenance)), settings.SPINUP_PROVENANCE)
            return
        end
    catch
    end
    delete(cell_file);
end
Run_Automatic_Site_Spinup(settings);
st = load_bg_state(cell_file);
Validate_BG_State_For_Transient(st, cell_file);
end

%% =======================================================================
function states = harmonise_states(cells)
%% The union of the states' fields, so they form one struct array.
names = {};
for c = 1:numel(cells); names = union(names,fieldnames(cells{c}),'stable'); end
for c = 1:numel(cells)
    for k = 1:numel(names)
        if ~isfield(cells{c},names{k}); cells{c}.(names{k}) = []; end
    end
    cells{c} = orderfields(cells{c},names);
end
states = vertcat(cells{:});
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
