%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bg_state = Attach_Spinup_Metadata(bg_state, spinup_converged, spinup_complete, ...
    spinup_start_year, spinup_end_year, spinup_cycles_completed, ...
    spinup_final_C_mean, spinup_final_C_max, spinup_final_N_mean, spinup_final_N_max, ...
    spinup_ion_names, spinup_ion_idx, spinup_final_ion_rel, spinup_mean_rel_tol, ...
    spinup_max_rel_tol, spinup_ion_rel_tol, spinup_seed_file, spinup_final_pass, ...
    spinup_required_consecutive_passes, spinup_final_consecutive_passes, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 18
    error('Attach_Spinup_Metadata:MissingInput', 'Spin-up metadata inputs are incomplete.');
end
if nargin < 19 || isempty(spinup_required_consecutive_passes)
    spinup_required_consecutive_passes = 1;
end
if nargin < 20 || isempty(spinup_final_consecutive_passes)
    spinup_final_consecutive_passes = double(spinup_final_pass);
end
if ~isstruct(bg_state)
    error('Attach_Spinup_Metadata:BadState', 'BG state must be a struct.');
end
if numel(spinup_ion_names) ~= numel(spinup_ion_idx) || numel(spinup_ion_idx) ~= numel(spinup_final_ion_rel)
    error('Attach_Spinup_Metadata:BadIonMetadata', 'Spin-up ion metadata lengths must match.');
end

bg_state.spinup_converged = logical(spinup_converged);
bg_state.spinup_complete = logical(spinup_complete);
bg_state.spinup_start_year = double(spinup_start_year);
bg_state.spinup_end_year = double(spinup_end_year);
bg_state.spinup_cycles_completed = double(spinup_cycles_completed);
bg_state.spinup_final_C_mean = double(spinup_final_C_mean);
bg_state.spinup_final_C_max = double(spinup_final_C_max);
bg_state.spinup_final_N_mean = double(spinup_final_N_mean);
bg_state.spinup_final_N_max = double(spinup_final_N_max);
bg_state.spinup_ion_names = string(spinup_ion_names);
bg_state.spinup_ion_idx = double(spinup_ion_idx);
bg_state.spinup_final_ion_rel = double(spinup_final_ion_rel);
if isempty(spinup_final_ion_rel)
    bg_state.spinup_final_ion_mean = NaN;
    bg_state.spinup_final_ion_max = NaN;
else
    bg_state.spinup_final_ion_mean = double(mean(spinup_final_ion_rel));
    bg_state.spinup_final_ion_max = double(max(spinup_final_ion_rel));
end
bg_state.spinup_mean_rel_tol = double(spinup_mean_rel_tol);
bg_state.spinup_max_rel_tol = double(spinup_max_rel_tol);
bg_state.spinup_ion_rel_tol = double(spinup_ion_rel_tol);
bg_state.spinup_seed_file = char(string(spinup_seed_file));
bg_state.spinup_final_pass = logical(spinup_final_pass);
bg_state.spinup_required_consecutive_passes = double(spinup_required_consecutive_passes);
bg_state.spinup_final_consecutive_passes = double(spinup_final_consecutive_passes);
if ~isempty(varargin) && isstruct(varargin{1})
    R = varargin{1};
    bg_state.spinup_erw_handoff_evaluated = logical(first_field(R, 'evaluated', false));
    bg_state.spinup_erw_handoff_ready = logical(first_field(R, 'ready', false));
    bg_state.spinup_erw_handoff_status = char(string(first_field(R, 'status', "not_evaluated")));
    bg_state.spinup_erw_handoff_message = char(string(first_field(R, 'message', '')));
    bg_state.spinup_erw_handoff_charge_residual_umolc_L = double(first_field(R, 'charge_residual_umolc_L', NaN));
    bg_state.spinup_erw_handoff_relative_charge_residual = double(first_field(R, 'relative_charge_residual', NaN));
    bg_state.spinup_erw_handoff_pH_pore = double(first_field(R, 'pH_pore_initial', NaN));
    bg_state.spinup_erw_handoff_pH_bulk = double(first_field(R, 'pH_bulk_predicted', NaN));
    bg_state.spinup_erw_handoff_diagnostic = first_field(R, 'diagnostic', struct());
end
end

function v = first_field(S, name, default)
v = default;
if isstruct(S) && isfield(S, name)
    v = S.(name);
end
end
