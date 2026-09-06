function targets = BG_Observation_Targets(MOD_PARAM_input)
%%% Workbook observations to g element m-2.

soil_mass_kg_m2 = MOD_PARAM_input.rsd*MOD_PARAM_input.Zbio/1000;
if ~isscalar(soil_mass_kg_m2) || ~isfinite(soil_mass_kg_m2) || soil_mass_kg_m2 <= 0
    error('BG_Observation_Targets:InvalidSoilMass', ...
        'The biogeochemical soil mass must be a positive finite scalar; rsd and Zbio give %s.', ...
        mat2str(soil_mass_kg_m2));
end

%%% Workbook column, target field and unit conversion.
%   g kg-1 * kg m-2                  -> g m-2
%   ug g-1 (= mg kg-1) * kg m-2/1000 -> g m-2
%   mass % * kg m-2 * 10             -> g m-2
%   cmol_c kg-1 * (1/100 mol_c/cmol_c) * equivalent weight [g mol_c-1]
%       * kg m-2                     -> g m-2   (exchangeable bases: the
%       lab's base saturation times its CEC; Ca 40.078/2, Mg 24.305/2,
%       K 39.098/1 g per mole of charge)
spec = {
    'InitPOC_gC_kg',      'POC_gCm2',      soil_mass_kg_m2
    'InitMAOC_gC_kg',     'MAOC_gCm2',     soil_mass_kg_m2
    'InitDOC_ugC_g',      'DOC_gCm2',      soil_mass_kg_m2/1000
    'InitMBC_ugC_g',      'MBC_gCm2',      soil_mass_kg_m2/1000
    'InitTN_pct',         'total_N_gNm2',  soil_mass_kg_m2*10
    'InitNH4_ugN_g',      'NH4_gNm2',      soil_mass_kg_m2/1000
    'InitNO3_ugN_g',      'NO3_gNm2',      soil_mass_kg_m2/1000
    'InitTotalP_mgP_kg',  'total_P_gPm2',  soil_mass_kg_m2/1000
    };

targets = struct();
for i = 1:size(spec, 1)
    column = spec{i,1};
    if ~ismember(column, MOD_PARAM_input.Properties.VariableNames)
        continue
    end
    value = MOD_PARAM_input.(column);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        continue
    end
    if value < 0
        error('BG_Observation_Targets:NegativeObservation', ...
            'Workbook column %s must be nonnegative; found %g.', column, value);
    end
    targets.(spec{i,2}) = value*spec{i,3};
end
end
