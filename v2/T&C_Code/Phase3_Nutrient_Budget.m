function budget = Phase3_Nutrient_Budget(terms, absolute_tolerance, relative_tolerance)
required = {'fertilizer','deposition','uplift','planting','fixation', ...
    'harvest','leaching','gaseous_loss','soil_start','soil_end', ...
    'plant_start','plant_end'};
missing = required(~isfield(terms, required));
if ~isempty(missing)
    error('Phase3_Nutrient_Budget:MissingTerm', ...
        'Missing nutrient budget term: %s.', strjoin(missing, ', '));
end
for i = 1:numel(required)
    value = terms.(required{i});
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
        error('Phase3_Nutrient_Budget:InvalidTerm', ...
            'Nutrient budget term %s must be a finite real scalar.', required{i});
    end
end
if ~(isnumeric(absolute_tolerance) && isreal(absolute_tolerance) && ...
        isscalar(absolute_tolerance) && isfinite(absolute_tolerance) && absolute_tolerance >= 0) || ...
        ~(isnumeric(relative_tolerance) && isreal(relative_tolerance) && ...
        isscalar(relative_tolerance) && isfinite(relative_tolerance) && relative_tolerance >= 0)
    error('Phase3_Nutrient_Budget:InvalidTolerance', ...
        'Nutrient budget tolerances must be finite non-negative scalars.');
end
budget.total_inputs = terms.fertilizer + terms.deposition + terms.uplift + ...
    terms.planting + terms.fixation;
budget.total_outputs = terms.harvest + terms.leaching + terms.gaseous_loss;
%%% A cell replaying its tracked lateral exchange (the flux-tracking field
%%% spin-up) also gains and loses solute laterally: lateral_in and
%%% lateral_out, both optional.
for lateral = {'lateral_in','lateral_out'}
    if isfield(terms,lateral{1})
        value = terms.(lateral{1});
        if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
            error('%s:InvalidTerm','Budget term %s must be a finite real scalar.', ...
                mfilename,lateral{1});
        end
    end
end
if isfield(terms,'lateral_in')
    budget.total_inputs = budget.total_inputs + terms.lateral_in;
end
if isfield(terms,'lateral_out')
    budget.total_outputs = budget.total_outputs + terms.lateral_out;
end
budget.storage_change = terms.soil_end - terms.soil_start + ...
    terms.plant_end - terms.plant_start;
budget.residual = budget.total_inputs - budget.total_outputs - budget.storage_change;
budget.scale = max([abs(budget.total_inputs), ...
    abs(budget.total_outputs) + abs(budget.storage_change), eps]);
budget.tolerance = max(absolute_tolerance, relative_tolerance*budget.scale);
budget.relative_residual = abs(budget.residual)/budget.scale;
budget.closed = abs(budget.residual) <= budget.tolerance;
end
