function diag = Validate_BG_Pool_Physicality(B_out, soil_mass_kg_m2, label)
% Validate biogeochemical pool shape and a broad total-phosphorus limit.

if nargin < 3 || isempty(label)
    label = 'BG state';
end
label = char(string(label));

valid_B_out = isnumeric(B_out) && isreal(B_out) && ...
    isequal(size(B_out), [1 79]) && all(isfinite(B_out)) && all(B_out >= 0);
if ~valid_B_out
    error('TC_SMEW_N_run:InvalidBGState', ...
        'BG state %s must contain 79 finite, non-negative areal pools.', label);
end
if ~(isnumeric(soil_mass_kg_m2) && isreal(soil_mass_kg_m2) && ...
        isscalar(soil_mass_kg_m2) && isfinite(soil_mass_kg_m2) && soil_mass_kg_m2 > 0)
    error('TC_SMEW_N_run:InvalidBGState', ...
        'BG state %s has an invalid biogeochemical soil mass.', label);
end

% Xu et al. (2020, doi:10.1002/saj2.20089) report about 400 mg total P kg-1
% in Midwestern surface soils. Five times that value is retained as a broad
% screening limit, not a calibration target. Exceeding it is a warning, not
% an error: occluded P (#46) accumulates by design under a persistent
% fertilizer surplus (Parton et al. 1988), and the ERW charge balance reads
% only the mineral-P pool #43, so a high total-P state does not corrupt the
% transient charge solution.
diag = struct();
diag.soil_mass_kg_m2 = soil_mass_kg_m2;
diag.max_total_P_mg_kg = 2000;
diag.total_P_mg_kg = 1000*sum(B_out(35:47))/soil_mass_kg_m2;
if diag.total_P_mg_kg > diag.max_total_P_mg_kg
    warning('TC_SMEW_N_run:TotalPExceedsScreen', ...
        ['BG state %s has %.6g mg total P kg-1 soil; the transient handoff ' ...
        'screening limit is %.6g mg P kg-1. This is the documented known ' ...
        'caveat of by-design occluded-P accumulation; the ERW charge ' ...
        'balance reads mineral P (#43) only and is unaffected.'], ...
        label, diag.total_P_mg_kg, diag.max_total_P_mg_kg);
end
end
