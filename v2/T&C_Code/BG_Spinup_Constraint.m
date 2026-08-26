function constraint = BG_Spinup_Constraint(MOD_PARAM_input, include_DOC)
%%% DOC remains prognostic during spinup, then uses the observed transient-start value.

if nargin < 2
    include_DOC = false;
end
if ~isscalar(include_DOC) || ~(islogical(include_DOC) || isnumeric(include_DOC)) || ...
        ~isfinite(double(include_DOC)) || ~any(double(include_DOC) == [0 1])
    error('BG_Spinup_Constraint:InvalidDOCOption', ...
        'The DOC option must be a logical scalar.');
end

excluded_pools = {'MBC_gCm2'};
if ~logical(include_DOC)
    excluded_pools{end+1} = 'DOC_gCm2';
end

targets = BG_Observation_Targets(MOD_PARAM_input);
constraint = struct();
constraint.targets = rmfield(targets, intersect( ...
    excluded_pools, fieldnames(targets)));
constraint.options = phosphorus_options(MOD_PARAM_input, constraint.targets);
if logical(include_DOC)
    constraint.source = 'workbook_transient_start_stocks';
else
    constraint.source = 'workbook_spinup_stocks';
end

function options = phosphorus_options(MOD_PARAM_input, targets)
options = struct('preserve_organic_P_pools', true);
if ~isfield(targets, 'total_P_gPm2')
    return
end

columns = {'InitOrganicP_C_to_P','InitP_solution_fraction', ...
    'InitP_primary_fraction','InitP_secondary_fraction', ...
    'InitP_occluded_fraction','InitSOM_C_to_P_min','InitSOM_C_to_P_max'};
present = ismember(columns, MOD_PARAM_input.Properties.VariableNames);
if ~any(present)
    return
end
if ~all(present)
    error('BG_Spinup_Constraint:IncompletePhosphorusOptions', ...
        'Total-P initialization requires all phosphorus partition options.');
end

values = nan(1,numel(columns));
for k = 1:numel(columns)
    value = MOD_PARAM_input.(columns{k});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
        error('BG_Spinup_Constraint:InvalidPhosphorusOptions', ...
            'Phosphorus partition options must be finite real numeric scalars.');
    end
    values(k) = value;
end
if values(1) <= 0 || any(values(2:5) < 0) || ...
        abs(sum(values(2:5)) - 1) > 1e-12 || values(6) <= 0 || values(7) <= values(6)
    error('BG_Spinup_Constraint:InvalidPhosphorusOptions', ...
        'Phosphorus partition options are outside their physical bounds.');
end

options = struct('preserve_organic_P_pools', false, ...
    'organic_P_C_to_P', values(1), ...
    'inorganic_P_fractions', values(2:5), ...
    'som_C_to_P_bounds', values(6:7));
end
end
