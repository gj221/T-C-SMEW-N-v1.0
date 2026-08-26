function D = ERW_Missing_Acidity_Diagnostic(E, C, pH_obs, pH_pred, Vw, CO2_air_extract)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ERW_MISSING_ACIDITY_DIAGNOSTIC Strong-acid anion needed to match bulk pH.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    Vw = NaN;
end
if nargin < 6
    CO2_air_extract = [];
end

D = struct('umolc_L', NaN, 'umolc_m2', NaN, 'target_pH', pH_obs, ...
    'pH_check', NaN, 'status', "not_run", 'warning', "");

if ~(isfinite(pH_obs))
    D.warning = "observed bulk pH is not finite";
    return
end
if ~(isfinite(pH_pred))
    try
        pH_pred = ERW_Bulk_pH_From_State(E, C, 0, CO2_air_extract);
    catch ME
        D.warning = string(ME.message);
        return
    end
end
if pH_pred <= pH_obs
    D.umolc_L = 0;
    D.pH_check = pH_pred;
    D.status = "none_needed";
    if isfinite(Vw) && Vw > 0
        D.umolc_m2 = 0;
    end
    return
end

try
    residual = @(dAn) ERW_Bulk_pH_From_State(E, C, dAn, CO2_air_extract) - pH_obs;
    r0 = residual(0);
    if ~(isfinite(r0) && r0 > 0)
        D.warning = "initial forward pH does not exceed observed pH";
        return
    end
    hi = 1;
    r_hi = residual(hi);
    while isfinite(r_hi) && r_hi > 0 && hi < 1e8
        hi = hi * 2;
        r_hi = residual(hi);
    end
    if ~(isfinite(r_hi) && r_hi <= 0)
        D.warning = "could not bracket observed bulk pH";
        return
    end
    lo = 0;
    for iter = 1:80
        mid = 0.5*(lo + hi);
        r_mid = residual(mid);
        if ~(isfinite(r_mid))
            hi = mid;
            continue
        end
        if abs(r_mid) <= 1e-5
            lo = mid;
            hi = mid;
            break
        end
        if r_mid > 0
            lo = mid;
        else
            hi = mid;
        end
    end
    D.umolc_L = 0.5*(lo + hi);
    D.pH_check = ERW_Bulk_pH_From_State(E, C, D.umolc_L, CO2_air_extract);
    D.status = "solved";
    if isfinite(Vw) && Vw > 0
        D.umolc_m2 = D.umolc_L * Vw;
    end
catch ME
    D.warning = string(ME.message);
end
end

function pH_bulk = ERW_Bulk_pH_From_State(E, C, An_extra, CO2_air_extract)
if nargin < 3 || isempty(An_extra), An_extra = 0; end
if nargin < 4 || isempty(CO2_air_extract)
    CO2_air_extract = co2_extract_boundary(E);
end
conv_mol = first_field(C, 'conv_mol');
if ~(isfinite(conv_mol) && conv_mol > 0), conv_mol = 1e6; end
pH_bulk = Porewater_to_Bulk_pH( ...
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

function CO2_air_extract = co2_extract_boundary(E)
CO2_air_extract = first_field(E, 'CO2_atm');
if isfinite(CO2_air_extract) && CO2_air_extract > 0, return; end
CO2_air_extract = first_field(E, 'CO2_air');
if isfinite(CO2_air_extract) && CO2_air_extract > 0, return; end
k_H = first_field(E, 'k_H'); CO2_w = first_field(E, 'CO2_w');
if isfinite(k_H) && k_H > 0 && isfinite(CO2_w) && CO2_w > 0
    CO2_air_extract = CO2_w / k_H;
end
end

function v = first_field(S, name)
v = NaN;
if ~isstruct(S) || ~isfield(S, name), return; end
x = S.(name);
if (isnumeric(x) || islogical(x)) && ~isempty(x)
    v = double(x(1));
end
end
