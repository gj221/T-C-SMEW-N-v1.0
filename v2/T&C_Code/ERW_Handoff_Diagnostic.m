function D = ERW_Handoff_Diagnostic(B_source, ERW_biogeochem, ERW_weathering, ERW_const, cov_moll_L_m2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transient ERW initialisation check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    cov_moll_L_m2 = NaN;
end

warnings = strings(1,0);
D = struct();
D.status = "pass";
D.warnings = warnings;
D.pH_bulk_initial = first_field(ERW_biogeochem, 'pH_bulk');
D.pH_pore_initial = first_field(ERW_biogeochem, 'pH');
D.pH_bulk_observed = observed_bulk_pH(D.pH_bulk_initial, ERW_const);
D.pH_bulk_predicted = NaN;
D.pH_bulk_residual = NaN;
D.missing_acidity_umolc_L = NaN;
D.missing_acidity_umolc_m2 = NaN;
D.missing_acidity_target_pH = NaN;
D.missing_acidity_pH_check = NaN;
D.base_saturation = NaN;
D.dissolved_Al_umol_L = NaN;
D.charge_residual_umolc_L = NaN;
D.relative_charge_residual = NaN;
D.f_CEC_sum = NaN;
D.f_CEC = NaN(1,6);
D.site_sum = NaN;
D.exchange_residuals = NaN;
D.exchange_mass_residuals_umol_L = NaN;
D.Omega_AlOH3 = first_field(ERW_weathering, 'Omega_AlOH3');
D.solution_conc_umol_L = empty_solution_struct();
D.source_solution_conc_umol_L = empty_solution_struct();
D.source_pools_g_m2 = empty_source_struct();
D.cov_moll_L_m2 = first_value(cov_moll_L_m2);
if isstruct(ERW_const) && isfield(ERW_const, 'bulk_pH_params') && ...
        isstruct(ERW_const.bulk_pH_params) && ...
        isfield(ERW_const.bulk_pH_params, 'field_water_L_m2')
    extract_field_water = first_value(ERW_const.bulk_pH_params.field_water_L_m2);
    volume_scale = max(abs([D.cov_moll_L_m2, extract_field_water]));
    if all(isfinite([D.cov_moll_L_m2, extract_field_water])) && ...
            abs(D.cov_moll_L_m2 - extract_field_water) > max(1e-9, 1e-8*volume_scale)
        warnings = add_warning(warnings, sprintf( ...
            'source and extract field-water volumes differ (%.6g vs %.6g L/m2)', ...
            D.cov_moll_L_m2, extract_field_water));
    end
end

[D.source_pools_g_m2, D.source_solution_conc_umol_L, warnings] = ...
    source_solution_concentrations(B_source, D.cov_moll_L_m2, warnings);
D.solution_conc_umol_L = equilibrated_solution_concentrations(ERW_biogeochem);

f_CEC = [first_field(ERW_biogeochem, 'f_Ca'), first_field(ERW_biogeochem, 'f_Mg'), ...
    first_field(ERW_biogeochem, 'f_K'), first_field(ERW_biogeochem, 'f_Na'), ...
    first_field(ERW_biogeochem, 'f_Al'), first_field(ERW_biogeochem, 'f_H')];
D.f_CEC = f_CEC;
if all(isfinite(f_CEC))
    D.f_CEC_sum = sum(f_CEC);
    D.base_saturation = sum(f_CEC(1:4));
    if abs(D.f_CEC_sum - 1) > 1e-3
        warnings = add_warning(warnings, sprintf('f_CEC sum %.6g differs from 1 by more than 1e-3', D.f_CEC_sum));
    end
else
    warnings = add_warning(warnings, 'f_CEC contains non-finite values');
end

D.dissolved_Al_umol_L = dissolved_al(ERW_biogeochem, ERW_const);
%%% Keep Al_w in model units for the transient exchange seed.
D.initial_Al_w_model_units = first_field(ERW_biogeochem, 'Al_w');
[D.charge_residual_umolc_L, D.relative_charge_residual, warnings] = ...
    charge_diagnostic(ERW_biogeochem, ERW_const, warnings);
[D.pH_bulk_predicted, D.pH_bulk_residual, extract_diag, warnings] = ...
    bulk_pH_diagnostic(ERW_biogeochem, ERW_const, D.pH_bulk_observed, warnings);
D.site_sum = first_field(extract_diag, 'site_sum');
D.exchange_residuals = numeric_field(extract_diag, 'exchange_residuals');
D.exchange_mass_residuals_umol_L = numeric_field(extract_diag, 'exchange_mass_residuals_umol_L');
missing_acidity = ERW_Missing_Acidity_Diagnostic(ERW_biogeochem, ERW_const, ...
    D.pH_bulk_observed, D.pH_bulk_predicted, D.cov_moll_L_m2);
D.missing_acidity_umolc_L = missing_acidity.umolc_L;
D.missing_acidity_umolc_m2 = missing_acidity.umolc_m2;
D.missing_acidity_target_pH = missing_acidity.target_pH;
D.missing_acidity_pH_check = missing_acidity.pH_check;

D.warnings = warnings;
if ~isempty(warnings)
    D.status = "warning";
end
end

function pH_obs = observed_bulk_pH(pH_bulk_initial, C)
pH_obs = pH_bulk_initial;
if isstruct(C) && isfield(C, 'bulk_pH_params') && isstruct(C.bulk_pH_params) && ...
        isfield(C.bulk_pH_params, 'observed_pH_initial') && isfinite(C.bulk_pH_params.observed_pH_initial)
    pH_obs = C.bulk_pH_params.observed_pH_initial;
end
end

function Al = dissolved_al(E, C)
Al = NaN;
Al_w = first_field(E, 'Al_w');
Al_scale = first_field(C, 'Al_CEC_scale');
if isfinite(Al_w) && isfinite(Al_scale) && Al_scale > 0
    Al = Al_w / Al_scale;
end
end

function [pH_bulk, pH_resid, extract_diag, warnings] = bulk_pH_diagnostic(E, C, pH_obs, warnings)
pH_bulk = NaN;
pH_resid = NaN;
extract_diag = struct();
req_e = ["H", "Ca", "Mg", "K", "Na", "Al_w", "Al", "An", "NO3", "PO4", "NH4", ...
    "k1", "k2", "k_w", "k_nh3", "k_H", "f_Ca", "f_Mg", "f_K", "f_Na", ...
    "f_Al", "f_H", "CEC_tot"];
req_c = ["K1", "K2", "K3", "K4", "bulk_pH_params", "rho_bulk", "Zrw", "Al_CEC_scale"];
if ~has_fields(E, req_e) || ~has_fields(C, req_c)
    warnings = add_warning(warnings, 'forward bulk pH diagnostic is missing required ERW fields');
    return
end
CO2_air_extract = co2_extract_boundary(E);
if ~(isfinite(CO2_air_extract) && CO2_air_extract > 0)
    warnings = add_warning(warnings, 'forward bulk pH diagnostic is missing a finite CO2 boundary');
    return
end
try
    [pH_bulk, extract_diag] = ERW_Bulk_pH_From_State(E, C, 0, CO2_air_extract);
    pH_resid = pH_bulk - pH_obs;
catch ME
    warnings = add_warning(warnings, sprintf('forward bulk pH diagnostic failed: %s', ME.message));
end
end

function CO2_air_extract = co2_extract_boundary(E)
CO2_air_extract = first_field(E, 'CO2_atm');
if isfinite(CO2_air_extract) && CO2_air_extract > 0
    return
end
CO2_air_extract = first_field(E, 'CO2_air');
if isfinite(CO2_air_extract) && CO2_air_extract > 0
    return
end
k_H = first_field(E, 'k_H');
CO2_w = first_field(E, 'CO2_w');
if isfinite(k_H) && k_H > 0 && isfinite(CO2_w) && CO2_w > 0
    CO2_air_extract = CO2_w / k_H;
end
end

function [S, C, warnings] = source_solution_concentrations(B_source, Vw, warnings)
names = ["NH4", "NO3", "PO4", "K", "Ca", "Mg", "Si"];
idx = [31, 32, 43, 52, 60, 68, 76];
mw = [14, 14, 31, 39.1, 40, 24.3, 28.1];
S = empty_source_struct();
C = empty_solution_struct();

if ~(isfinite(Vw) && Vw > 0)
    warnings = add_warning(warnings, 'cov_moll_L_m2 must be finite and positive');
    return
end
if ~isnumeric(B_source) || numel(B_source) < max(idx)
    warnings = add_warning(warnings, 'source B pool vector is missing required mobile solution pools');
    return
end

B = B_source(:).';
for ii = 1:numel(names)
    field = char(names(ii));
    pool_g_m2 = B(idx(ii));
    S.(field) = pool_g_m2;
    C.(field) = pool_g_m2 / mw(ii) * 1e6 / Vw;
end

vals = cell2mat(struct2cell(C));   %% struct2array is gone from recent releases
if any(~isfinite(vals)) || any(vals < 0)
    warnings = add_warning(warnings, 'source-derived solution concentrations must be finite and non-negative');
end
end

function [R, rel_R, warnings] = charge_diagnostic(E, C, warnings)
R = NaN;
rel_R = NaN;
%%% The balance is taken over the species the state was solved with.  v0.9
%%% speciates the pore water with activity-conditioned constants and retains
%%% every charged species (ERW_MOD_PARAM_excel, ERW_function); re-speciating
%%% Al with the thermodynamic K1..K4 of ERW_const against the retained Al3+
%%% would report the activity correction itself as an imbalance (0.69
%%% umolc/L on a Kantola-pool Bruce Farm cell whose species close to 2e-12).
%%% A state that does not retain them is re-speciated as before.
held = ["H", "Ca", "Mg", "Na", "K", "NH4", "Al", "AlOH", "AlOH2", "AlOH4", ...
    "k_w", "HCO3", "CO3", "NO3", "PO4", "An"];
if has_fields(E, held) && has_fields(C, "Al_CEC_scale")
    Al_scale = first_field(C, 'Al_CEC_scale');
    H = first_field(E, 'H');
    cations = H + 2*first_field(E, 'Ca') + 2*first_field(E, 'Mg') + ...
        first_field(E, 'Na') + first_field(E, 'K') + first_field(E, 'NH4') + ...
        (3*first_field(E, 'Al') + 2*first_field(E, 'AlOH') + ...
        first_field(E, 'AlOH2'))/Al_scale;
    anions = first_field(E, 'k_w')/H + first_field(E, 'HCO3') + ...
        2*first_field(E, 'CO3') + first_field(E, 'NO3') + ...
        first_field(E, 'PO4') + first_field(E, 'An') + ...
        first_field(E, 'AlOH4')/Al_scale;
    R = cations - anions;
    rel_R = abs(R) / max(abs(cations) + abs(anions), eps);
    if abs(R) > 1e-6 && rel_R > 1e-8
        warnings = add_warning(warnings, sprintf('charge residual %.3g umolc/L exceeds tolerance', R));
    end
    return
end
req_e = ["H", "CO2_w", "Ca", "Mg", "Na", "K", "NH4", "An", "NO3", "PO4", ...
    "Al_w", "Al", "k1", "k2", "k_w", "k_nh3"];
req_c = ["K1", "K2", "K3", "K4", "Al_CEC_scale"];
if ~has_fields(E, req_e) || ~has_fields(C, req_c)
    warnings = add_warning(warnings, 'charge residual diagnostic is missing required ERW fields');
    return
end

Al_scale = first_field(C, 'Al_CEC_scale');
[R, state] = ERW_Aqueous_Charge_Residual( ...
    first_field(E, 'H'), first_field(E, 'CO2_w'), first_field(E, 'Ca'), ...
    first_field(E, 'Mg'), first_field(E, 'Na'), first_field(E, 'K'), ...
    first_field(E, 'NH4'), first_field(E, 'An'), first_field(E, 'NO3'), ...
    first_field(E, 'PO4'), first_field(E, 'Al_w'), first_field(E, 'k1'), ...
    first_field(E, 'k2'), first_field(E, 'k_w'), first_field(E, 'k_nh3'), first_field(C, 'K1'), ...
    first_field(C, 'K2'), first_field(C, 'K3'), first_field(C, 'K4'), ...
    first_field(E, 'Al'), Al_scale);
charge_scale = max(abs(state.cations) + abs(state.anions), eps);
rel_R = abs(R) / charge_scale;
if abs(R) > 1e-6 && rel_R > 1e-8
    warnings = add_warning(warnings, sprintf('charge residual %.3g umolc/L exceeds tolerance', R));
end
end

function S = empty_solution_struct()
S = struct('NH4', NaN, 'NO3', NaN, 'PO4', NaN, 'K', NaN, 'Ca', NaN, 'Mg', NaN, 'Si', NaN);
end

function C = equilibrated_solution_concentrations(E)
C = empty_solution_struct();
names = fieldnames(C);
for ii = 1:numel(names)
    C.(names{ii}) = first_field(E, names{ii});
end
end

function S = empty_source_struct()
S = struct('NH4', NaN, 'NO3', NaN, 'PO4', NaN, 'K', NaN, 'Ca', NaN, 'Mg', NaN, 'Si', NaN);
end

function tf = has_fields(S, names)
tf = isstruct(S) && all(isfield(S, cellstr(names)));
end

function v = first_field(S, name)
v = NaN;
if ~isstruct(S) || ~isfield(S, name)
    return
end
v = first_value(S.(name));
end

function v = numeric_field(S, name)
v = NaN;
if isstruct(S) && isfield(S, name)
    x = S.(name);
    if isnumeric(x) || islogical(x)
        x = x(:).';
        if ~isempty(x)
            v = double(x);
        end
    end
end
end

function v = first_value(x)
v = NaN;
if isnumeric(x) || islogical(x)
    x = x(:);
    if ~isempty(x)
        v = double(x(1));
    end
end
end

function warnings = add_warning(warnings, msg)
warnings(end+1) = string(msg);
end

function [R, state] = ERW_Aqueous_Charge_Residual(H, CO2_w, Ca, Mg, ...
    Na, K_aq, NH4, An, NO3, PO4, Al_w, k1, k2, k_w, ~, K1, K2, K3, K4, ...
    Al3_in, Al_CEC_scale)
if nargin < 21, Al3_in = []; end
if nargin < 21 || isempty(Al_CEC_scale) || ~isfinite(Al_CEC_scale) || Al_CEC_scale <= 0
    Al_CEC_scale = 1;
end
Al_CEC_scale = max(real(Al_CEC_scale), 1);
H = max(real(H), eps); CO2_w = max(real(CO2_w), 0);
Ca = max(real(Ca), 0); Mg = max(real(Mg), 0);
Na = max(real(Na), 0); K_aq = max(real(K_aq), 0);
NH4 = max(real(NH4), 0); An = max(real(An), 0);
NO3 = max(real(NO3), 0); PO4 = max(real(PO4), 0);
Al_w = max(real(Al_w), 0);
Al_w_aq = Al_w / Al_CEC_scale;
OH = k_w./H; HCO3 = k1.*CO2_w./H; CO3 = k2.*k1.*CO2_w./(H.^2);
Al_den = max(real(H.^4 + H.^3.*K1 + H.^2.*K1.*K2 + H.*K1.*K2.*K3 + K1.*K2.*K3.*K4), eps);
Al3_calc = H.^4./Al_den .* Al_w_aq;
AlOH  = H.^3.*K1./Al_den .* Al_w_aq;
AlOH2 = H.^2.*K1.*K2./Al_den .* Al_w_aq;
AlOH3 = H.*K1.*K2.*K3./Al_den .* Al_w_aq;
AlOH4 = max(Al_w_aq - (Al3_calc + AlOH + AlOH2 + AlOH3), 0);
if isempty(Al3_in), Al3 = Al3_calc;
else, Al3 = max(real(Al3_in), 0) / Al_CEC_scale; end
% Retain the legacy TC-SMEW strong-ion convention used by Alk_tot.
NH4_charge = NH4;
cations = H + 2*Ca + 2*Mg + Na + K_aq + NH4_charge + 3*Al3 + 2*AlOH + AlOH2;
anions  = OH + HCO3 + 2*CO3 + An + NO3 + PO4 + AlOH4;
R = cations - anions;
if nargout > 1
    state = struct('OH',OH,'HCO3',HCO3,'CO3',CO3,'Al3',Al3,'AlOH',AlOH, ...
        'AlOH2',AlOH2,'AlOH3',AlOH3,'AlOH4',AlOH4,'Al_w_aq',Al_w_aq, ...
        'Al_CEC_scale',Al_CEC_scale,'cations',cations,'anions',anions);
end
end

function [pH_bulk, extract_diag] = ERW_Bulk_pH_From_State(E, C, An_extra, CO2_air_extract)
if nargin < 3 || isempty(An_extra), An_extra = 0; end
if nargin < 4 || isempty(CO2_air_extract)
    CO2_air_extract = co2_extract_boundary(E);
end
conv_mol = first_field(C, 'conv_mol');
if ~(isfinite(conv_mol) && conv_mol > 0), conv_mol = 1e6; end
[pH_bulk, extract_diag] = Porewater_to_Bulk_pH( ...
    first_field(E,'H'), first_field(E,'Ca'), first_field(E,'Mg'), ...
    first_field(E,'K'), first_field(E,'Na'), first_field(E,'Al_w'), ...
    first_field(E,'Al'), first_field(E,'An') + An_extra, first_field(E,'NO3'), ...
    first_field(E,'PO4'), first_field(E,'NH4'), CO2_air_extract, ...
    first_field(E,'k1'), first_field(E,'k2'), first_field(E,'k_w'), ...
    first_field(E,'k_nh3'), first_field(E,'k_H'), first_field(C,'K1'), ...
    first_field(C,'K2'), first_field(C,'K3'), first_field(C,'K4'), ...
    first_field(E,'f_Ca'), first_field(E,'f_Mg'), first_field(E,'f_K'), ...
    first_field(E,'f_Na'), first_field(E,'f_Al'), first_field(E,'f_H'), ...
    first_field(E,'CEC_tot'), first_field(C,'rho_bulk'), first_field(C,'Zrw'), ...
    conv_mol, C.bulk_pH_params, first_field(C,'Al_CEC_scale'));
end
