function budget = Phase3_Phosphorus_Budget(terms, absolute_tolerance, relative_tolerance)
%%% 'litter' is the prescribed litter-nutrient return.  It is an external
%%% input only where the vegetation is prescribed rather than simulated: the
%%% standalone phase takes it from the Phase 1 averages, so it enters the soil
%%% with no plant pool to debit.  Where plants are simulated the same flux is
%%% internal and the term is zero.  Required, not defaulted, so a caller must
%%% state which case it is (cf. legacy Biogeo_Spin_Up.m, which carries
%%% sum(IS(13:15)) in P_exp alongside deposition, uplift and fertiliser).
required = {'fertilizer','deposition','uplift','planting','litter','harvest','leaching', ...
    'soil_start','soil_end','plant_start','plant_end'};
missing = required(~isfield(terms,required));
if ~isempty(missing)
    error('Phase3_Phosphorus_Budget:MissingTerm', ...
        'Missing phosphorus budget term: %s.', strjoin(missing,', '));
end

values = zeros(size(required));
for i = 1:numel(required)
    value = terms.(required{i});
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
        error('Phase3_Phosphorus_Budget:InvalidTerm', ...
            'Phosphorus budget term %s must be a finite real scalar.', required{i});
    end
    values(i) = value;
end
if ~(isnumeric(absolute_tolerance) && isreal(absolute_tolerance) && ...
        isscalar(absolute_tolerance) && isfinite(absolute_tolerance) && absolute_tolerance >= 0) || ...
        ~(isnumeric(relative_tolerance) && isreal(relative_tolerance) && ...
        isscalar(relative_tolerance) && isfinite(relative_tolerance) && relative_tolerance >= 0)
    error('Phase3_Phosphorus_Budget:InvalidTolerance', ...
        'Phosphorus budget tolerances must be finite non-negative scalars.');
end

budget.total_inputs = terms.fertilizer + terms.deposition + terms.uplift + ...
    terms.planting + terms.litter;
budget.total_outputs = terms.harvest + terms.leaching;
budget.storage_change = terms.soil_end - terms.soil_start + terms.plant_end - terms.plant_start;
budget.residual = budget.total_inputs - budget.total_outputs - budget.storage_change;
budget.scale = max([abs(budget.total_inputs), abs(budget.total_outputs) + abs(budget.storage_change), eps]);
budget.tolerance = max(absolute_tolerance, relative_tolerance*budget.scale);
budget.relative_residual = abs(budget.residual)/budget.scale;
budget.closed = abs(budget.residual) <= budget.tolerance;
