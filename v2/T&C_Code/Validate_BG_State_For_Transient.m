%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Validate_BG_State_For_Transient(bg_state, bg_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2 || isempty(bg_file)
    bg_file = 'BG state';
end
if ~isstruct(bg_state)
    error('TC_SMEW_N_run:InvalidBGState', 'BG state %s is not a struct.', char(string(bg_file)));
end

required_metadata = {'spinup_converged','spinup_complete','spinup_final_pass', ...
    'spinup_erw_handoff_evaluated','spinup_erw_handoff_ready', ...
    'spinup_nutrient_budget_vers','spinup_nutrient_budgets', ...
    'spinup_biogeochemical_soil_mass_kg_m2'};
missing_metadata = required_metadata(~isfield(bg_state, required_metadata));
if ~isempty(missing_metadata)
    error('TC_SMEW_N_run:MissingSpinupMetadata', ...
        'BG state %s is missing required spin-up metadata: %s.', ...
        char(string(bg_file)), strjoin(missing_metadata, ', '));
end

if ~scalar_true(bg_state.spinup_converged)
    error('TC_SMEW_N_run:UnconvergedBGState', ...
        'BG state %s was produced by an unconverged spin-up and cannot seed a transient run.', ...
        char(string(bg_file)));
end
if ~scalar_true(bg_state.spinup_complete)
    error('TC_SMEW_N_run:UnconvergedBGState', ...
        'BG state %s did not complete spin-up and cannot seed a transient run.', ...
        char(string(bg_file)));
end
if ~scalar_true(bg_state.spinup_final_pass)
    error('TC_SMEW_N_run:UnconvergedBGState', ...
        'BG state %s did not pass the final spin-up diagnostic and cannot seed a transient run.', ...
        char(string(bg_file)));
end
if ~scalar_true(bg_state.spinup_erw_handoff_evaluated)
    error('TC_SMEW_N_run:UnconvergedBGState', ...
        'BG state %s did not evaluate the ERW charge-balance/pH handoff readiness gate.', ...
        char(string(bg_file)));
end
if ~scalar_true(bg_state.spinup_erw_handoff_ready)
    error('TC_SMEW_N_run:UnconvergedBGState', ...
        'BG state %s failed the ERW charge-balance/pH handoff readiness gate.', ...
        char(string(bg_file)));
end

budget_version = bg_state.spinup_nutrient_budget_vers;
if ~(isnumeric(budget_version) && isreal(budget_version) && ...
        isscalar(budget_version) && isfinite(budget_version) && budget_version == 12)
    error('TC_SMEW_N_run:InvalidBGState', ...
        'BG state %s has an incompatible nutrient-budget version.', char(string(bg_file)));
end
budgets = bg_state.spinup_nutrient_budgets;
elements = {'N','P','K','Ca','Mg','Si'};
if ~isstruct(budgets) || ~all(isfield(budgets, elements))
    error('TC_SMEW_N_run:InvalidBGState', ...
        'BG state %s does not contain all six nutrient budgets.', char(string(bg_file)));
end
for e = 1:numel(elements)
    budget = budgets.(elements{e});
    required = {'closed','residual','tolerance'};
    if ~isstruct(budget) || ~all(isfield(budget, required))
        error('TC_SMEW_N_run:InvalidBGState', ...
            'BG state %s has an open or invalid %s budget.', ...
            char(string(bg_file)), elements{e});
    end
    valid_closed = (islogical(budget.closed) || isnumeric(budget.closed)) && ...
        isreal(budget.closed) && isscalar(budget.closed) && ...
        isfinite(double(budget.closed)) && logical(budget.closed);
    valid_residual = isnumeric(budget.residual) && isreal(budget.residual) && ...
        isscalar(budget.residual) && isfinite(budget.residual);
    valid_tolerance = isnumeric(budget.tolerance) && isreal(budget.tolerance) && ...
        isscalar(budget.tolerance) && isfinite(budget.tolerance) && budget.tolerance >= 0;
    if ~valid_closed || ~valid_residual || ~valid_tolerance || ...
            abs(budget.residual) > budget.tolerance
        error('TC_SMEW_N_run:InvalidBGState', ...
            'BG state %s has an open or invalid %s budget.', ...
            char(string(bg_file)), elements{e});
    end
end

if ~isfield(bg_state, 'B_out')
    error('TC_SMEW_N_run:InvalidBGState', ...
        'BG state %s is missing the biogeochemical pool vector.', char(string(bg_file)));
end
Validate_BG_Pool_Physicality(bg_state.B_out, ...
    bg_state.spinup_biogeochemical_soil_mass_kg_m2, bg_file);

end
function tf = scalar_true(x)
if islogical(x) || isnumeric(x)
    tf = isscalar(x) && isfinite(double(x)) && logical(x);
elseif isstring(x) || ischar(x)
    txt = strtrim(char(string(x)));
    tf = any(strcmpi(txt, {'true','1','yes','on'}));
else
    tf = false;
end
end
