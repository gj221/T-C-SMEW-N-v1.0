function [H_seed, Na_seed, Al_w_seed, Al_seed, f_Ca_seed, f_Mg_seed, ...
    f_K_seed, f_Na_seed, f_Al_seed, f_H_seed] = ...
    ERW_Initial_H_From_Anion_Budget(An_target, H_guess, Ca, Mg, Na_in, K_aq, ...
    NO3, PO4, NH4, CO2_w, k1, k2, k_w, k_nh3, ...
    K1, K2, K3, K4, K_CEC, K_AlOH3, Al_CEC_scale, conv_mol)
%%% Solve the pore-water H state from the independently sourced strong-anion
%%% budget.  Observed bulk pH is not an input.
%%%
%%% Ported verbatim from Model/TC-SMEW-N branch lian-spinup-convergence
%%% (commit b8ba95b, "Predict calibration pH from independent anion
%%% budget").  Its dependencies -- ERW_Residual_An_From_Charge,
%%% Initialize_ERW_Equilibrium, Exchange_Fractions_GT -- already exist in
%%% this tree with matching signatures.  Used by Close_ERW_Charge_Balance
%%% to re-derive a charge-balanced exchange/pH state for RF-reconstructed
%%% cells (Reconstruct_Cell_Pools_RF).

if ~(isscalar(An_target) && isfinite(An_target) && An_target >= 0)
    error('ERW:invalidState', ...
        'Initial sourced An concentration must be finite and non-negative.');
end
if ~(isscalar(Na_in) && isfinite(Na_in) && Na_in >= 0)
    Na_in = 0;
end
concentrations = [Ca, Mg, Na_in, K_aq, NO3, PO4, NH4, CO2_w];
if numel(concentrations) ~= 8 || any(~isfinite(concentrations)) || ...
        any(concentrations < 0) || Ca <= 0
    error('ERW:invalidState', ...
        'Initial anion-budget concentrations must be finite and non-negative, with positive Ca.');
end
equilibrium_constants = [k1, k2, k_w, k_nh3, K1, K2, K3, K4];
if numel(equilibrium_constants) ~= 8 || any(~isfinite(equilibrium_constants)) || ...
        any(equilibrium_constants <= 0) || ...
        numel(K_CEC) ~= 5 || any(~isfinite(K_CEC)) || any(K_CEC <= 0) || ...
        ~(isscalar(K_AlOH3) && isfinite(K_AlOH3)) || ...
        ~(isscalar(Al_CEC_scale) && isfinite(Al_CEC_scale) && Al_CEC_scale > 0) || ...
        ~(isscalar(conv_mol) && isfinite(conv_mol) && conv_mol > 0)
    error('ERW:invalidParam', ...
        'Initial anion-budget equilibrium parameters must be positive and finite.');
end

residual = @(H) initial_anion_residual(H, An_target, Ca, Mg, Na_in, K_aq, ...
    NO3, PO4, NH4, CO2_w, k1, k2, k_w, k_nh3, ...
    K1, K2, K3, K4, K_AlOH3, Al_CEC_scale, conv_mol);

pH_scan = linspace(2, 12, 201);
H_scan = 10.^(-pH_scan) * conv_mol;
R_scan = nan(size(H_scan));
for ii = 1:numel(H_scan)
    try
        R_scan(ii) = residual(H_scan(ii));
    catch ME
        if startsWith(ME.identifier, 'ERW:')
            R_scan(ii) = NaN;
        else
            rethrow(ME)
        end
    end
end
valid = isfinite(R_scan);
bracket_idx = find(valid(1:end-1) & valid(2:end) & ...
    (R_scan(1:end-1) == 0 | R_scan(2:end) == 0 | ...
     sign(R_scan(1:end-1)) ~= sign(R_scan(2:end))));
if isempty(bracket_idx)
    error('ERW:noInitialPHRoot', ...
        'Could not bracket initial porewater H from the sourced anion budget.');
end
if isscalar(H_guess) && isfinite(H_guess) && H_guess > 0
    pH_guess = -log10(H_guess / conv_mol);
else
    pH_guess = 7;
end
[~, i_best] = min(abs(0.5*(pH_scan(bracket_idx) + ...
    pH_scan(bracket_idx+1)) - pH_guess));
lo = H_scan(bracket_idx(i_best));
hi = H_scan(bracket_idx(i_best)+1);
R_lo = residual(lo);
R_hi = residual(hi);
if abs(R_lo) <= 1e-8
    H_seed = lo;
elseif abs(R_hi) <= 1e-8
    H_seed = hi;
else
    for iter = 1:100
        mid = 0.5*(lo + hi);
        R_mid = residual(mid);
        if abs(R_mid) <= max(1e-8, 1e-10*max(An_target, 1))
            break
        end
        if sign(R_lo) == sign(R_mid)
            lo = mid;
            R_lo = R_mid;
        else
            hi = mid;
        end
    end
    H_seed = 0.5*(lo + hi);
end

Al_w_seed = initial_al_from_aloh3( ...
    H_seed, K1, K2, K3, K4, K_AlOH3, Al_CEC_scale, conv_mol);
Al_den = H_seed^4 + H_seed^3*K1 + H_seed^2*K1*K2 + ...
    H_seed*K1*K2*K3 + K1*K2*K3*K4;
Al_seed = (H_seed^4 / Al_den) * Al_w_seed;
[f_Ca_seed, f_Mg_seed, f_K_seed, f_Na_seed, f_Al_seed, f_H_seed] = ...
    Exchange_Fractions_GT(H_seed, Ca, Mg, K_aq, Na_in, Al_w_seed, ...
    K_CEC, K1, K2, K3, K4, Al_CEC_scale);
Na_seed = Na_in;
end

function R = initial_anion_residual(H, An_target, Ca, Mg, Na, K_aq, ...
    NO3, PO4, NH4, CO2_w, k1, k2, k_w, k_nh3, ...
    K1, K2, K3, K4, K_AlOH3, Al_CEC_scale, conv_mol)
Al_w = initial_al_from_aloh3( ...
    H, K1, K2, K3, K4, K_AlOH3, Al_CEC_scale, conv_mol);
An_resid = ERW_Residual_An_From_Charge(H, CO2_w, Ca, Mg, Na, K_aq, ...
    NH4, NO3, PO4, Al_w, k1, k2, k_w, k_nh3, ...
    K1, K2, K3, K4, Al_CEC_scale);
R = An_resid - An_target;
end

function Al_w = initial_al_from_aloh3( ...
    H, K1, K2, K3, K4, K_AlOH3, Al_CEC_scale, conv_mol)
Al_w = Initialize_ERW_Equilibrium(H, K1, K2, K3, K4, ...
    conv_mol, K_AlOH3, Al_CEC_scale);
if ~(isfinite(Al_w) && Al_w > 0)
    error('ERW:invalidState', 'Cannot initialize Al_w from Al(OH)3 equilibrium.');
end
end
