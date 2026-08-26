function R = ERW_Handoff_Readiness(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ERW_HANDOFF_READINESS Hard gate for pH/charge-balance handoff safety.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = default_readiness();
if nargin < 1 || isempty(D)
    return
end
if ~isstruct(D)
    error('ERW_Handoff_Readiness:BadInput', 'ERW handoff diagnostic must be a struct.');
end
R.diagnostic = D;

R.evaluated = true;
R.status = "pass";
warnings = strings(1,0);
charge_abs = abs(first_field(D, 'charge_residual_umolc_L'));
charge_rel = abs(first_field(D, 'relative_charge_residual'));
R.charge_residual_umolc_L = first_field(D, 'charge_residual_umolc_L');
R.relative_charge_residual = charge_rel;
if ~(isfinite(charge_abs) && isfinite(charge_rel) && ...
        (charge_abs <= R.charge_abs_tol_umolc_L || charge_rel <= R.charge_rel_tol))
    warnings = add_warning(warnings, sprintf( ...
        'charge residual exceeds readiness tolerance (abs %.3g umolc/L, rel %.3g)', ...
        charge_abs, charge_rel));
end

pH_pore = first_field(D, 'pH_pore_initial');
pH_bulk = first_field(D, 'pH_bulk_predicted');
R.pH_pore_initial = pH_pore;
R.pH_bulk_predicted = pH_bulk;
if ~(isfinite(pH_pore) && pH_pore > R.pH_min && pH_pore < R.pH_max)
    warnings = add_warning(warnings, sprintf('pore pH %.3g is outside readiness range', pH_pore));
end
if ~(isfinite(pH_bulk) && pH_bulk > R.pH_min && pH_bulk < R.pH_max)
    warnings = add_warning(warnings, sprintf('bulk pH %.3g is outside readiness range', pH_bulk));
end

f_sum = first_field(D, 'f_CEC_sum');
R.f_CEC_sum = f_sum;
if ~(isfinite(f_sum) && abs(f_sum - 1) <= R.f_CEC_sum_tol)
    warnings = add_warning(warnings, sprintf('f_CEC sum %.6g fails readiness tolerance', f_sum));
end

site_sum = first_field(D, 'site_sum');
R.site_sum = site_sum;
if ~(isfinite(site_sum) && abs(site_sum - 1) <= R.site_sum_tol)
    warnings = add_warning(warnings, ...
        sprintf('forward extract site closure %.6g fails readiness tolerance', site_sum));
end

exchange_residuals = numeric_field(D, 'exchange_residuals');
R.exchange_residuals = exchange_residuals;
exchange_residual_max = max_abs_residual(exchange_residuals);
if ~(isfinite(exchange_residual_max) && exchange_residual_max <= R.exchange_residual_tol)
    warnings = add_warning(warnings, sprintf( ...
        'forward extract exchange equilibrium max residual %.3g fails readiness tolerance', ...
        exchange_residual_max));
end

exchange_mass_residuals = numeric_field(D, 'exchange_mass_residuals_umol_L');
R.exchange_mass_residuals_umol_L = exchange_mass_residuals;
exchange_mass_residual_max = max_abs_residual(exchange_mass_residuals);
if ~(isfinite(exchange_mass_residual_max) && exchange_mass_residual_max <= R.exchange_mass_residual_tol_umol_L)
    warnings = add_warning(warnings, sprintf( ...
        'forward extract exchange mass max residual %.3g umol/L fails readiness tolerance', ...
        exchange_mass_residual_max));
end

% Bulk pH residual: predicted vs observed
pH_bulk_obs  = first_field(D, 'pH_bulk_observed');
R.pH_bulk_observed = pH_bulk_obs;
if isfinite(pH_bulk_obs) && isfinite(pH_bulk)
    pH_resid = pH_bulk - pH_bulk_obs;
    R.pH_bulk_residual = pH_resid;
    if abs(pH_resid) > R.pH_bulk_residual_tol
        warnings = add_warning(warnings, ...
            sprintf('bulk pH residual %.4g exceeds tolerance %.4g', pH_resid, R.pH_bulk_residual_tol));
    end
end

% Dissolved Ca sanity check (umol/L from diagnostic struct)
Ca_conc = NaN;
if isstruct(D) && isfield(D, 'solution_conc_umol_L')
    Ca_conc = first_field(D.solution_conc_umol_L, 'Ca');
end
Ca_mM = Ca_conc / 1000;
R.dissolved_Ca_mM = Ca_mM;
if isfinite(Ca_mM)
    if Ca_mM < R.Ca_sanity_min_mM
        warnings = add_warning(warnings, ...
            sprintf('dissolved Ca %.4g mM below sanity minimum %.4g mM', Ca_mM, R.Ca_sanity_min_mM));
    elseif Ca_mM > R.Ca_sanity_max_mM
        warnings = add_warning(warnings, ...
            sprintf('dissolved Ca %.4g mM exceeds sanity maximum %.4g mM', Ca_mM, R.Ca_sanity_max_mM));
    end
end

%%% Pore-water and extract pH are separate diagnostics.
if isfinite(pH_pore) && isfinite(pH_bulk_obs)
    R.pH_pore_bulk_separation = abs(pH_pore - pH_bulk_obs);
end

Al_umol_L = first_field(D, 'dissolved_Al_umol_L');
R.dissolved_Al_umol_L = Al_umol_L;
if ~(isfinite(Al_umol_L) && Al_umol_L >= 0)
    warnings = add_warning(warnings, ...
        sprintf('dissolved Al %.4g umol/L must be finite and non-negative', Al_umol_L));
end

BS = first_field(D, 'base_saturation');
R.base_saturation = BS;
if ~(isfinite(BS) && BS >= 0 && BS <= 1)
    warnings = add_warning(warnings, ...
        sprintf('base saturation %.4g must be within the physical range [0, 1]', BS));
end

R.warnings = warnings;
R.ready = isempty(warnings);
if ~R.ready
    R.status = "fail";
    R.message = char(strjoin(warnings, '; '));
else
    R.message = 'ERW handoff charge balance and pH diagnostics passed.';
end
end

function R = default_readiness()
R = struct();
R.evaluated = false;
R.ready = false;
R.status = "not_evaluated";
R.message = 'ERW handoff readiness was not evaluated.';
R.diagnostic = struct();
R.warnings = strings(1,0);
R.charge_abs_tol_umolc_L = 1e-6;
R.charge_rel_tol = 1e-8;
R.pH_min = 2;
R.pH_max = 12;
R.pH_bulk_residual_tol = 0.15;  % |predicted - observed| bulk pH
R.Ca_sanity_min_mM = 0.01;
R.Ca_sanity_max_mM = 20;
R.f_CEC_sum_tol = 1e-3;
R.site_sum_tol = 1e-8;
R.exchange_residual_tol = 1e-8;
R.exchange_mass_residual_tol_umol_L = 1e-8;
R.charge_residual_umolc_L = NaN;
R.relative_charge_residual = NaN;
R.pH_pore_initial = NaN;
R.pH_bulk_predicted = NaN;
R.pH_bulk_observed = NaN;
R.pH_bulk_residual = NaN;
R.pH_pore_bulk_separation = NaN;
R.dissolved_Ca_mM = NaN;
R.dissolved_Al_umol_L = NaN;
R.base_saturation = NaN;
R.f_CEC_sum = NaN;
R.site_sum = NaN;
R.exchange_residuals = NaN;
R.exchange_mass_residuals_umol_L = NaN;
end

function v = first_field(S, name)
v = NaN;
if isstruct(S) && isfield(S, name)
    x = S.(name);
    if isnumeric(x) || islogical(x)
        x = x(:);
        if ~isempty(x)
            v = double(x(1));
        end
    end
end
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

function max_abs = max_abs_residual(v)
v = v(:);
if isempty(v) || any(~isfinite(v))
    max_abs = NaN;
else
    max_abs = max(abs(v));
end
end

function warnings = add_warning(warnings, msg)
warnings(end+1) = string(msg);
end
