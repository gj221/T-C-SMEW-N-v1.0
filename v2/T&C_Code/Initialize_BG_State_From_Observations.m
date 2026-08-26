function [B, diag] = Initialize_BG_State_From_Observations(B, targets, options)
% Map available observed areal stocks onto the legacy 79-pool state.

if ~isnumeric(B) || ~isreal(B) || ~isequal(size(B), [1 79]) || any(~isfinite(B))
    error('Initialize_BG_State_From_Observations:InvalidPoolVector', ...
        'The initial biogeochemistry state must be a finite real 1-by-79 vector.');
end
if any(B < -1e-6)
    error('Initialize_BG_State_From_Observations:NegativePool', ...
        'The legacy initial state contains a materially negative pool.');
end
B = max(B, 0);
initial_total_N = sum(B(23:34));
initial_total_P = sum(B(35:47));

if ~isstruct(targets) || ~isscalar(targets)
    error('Initialize_BG_State_From_Observations:InvalidTargets', ...
        'Observed stock targets must be supplied as a scalar struct.');
end
if nargin < 3 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('Initialize_BG_State_From_Observations:InvalidOptions', ...
        'Initialisation options must be supplied as a scalar struct.');
end

recognized_targets = {'POC_gCm2','MAOC_gCm2','DOC_gCm2','MBC_gCm2', ...
    'total_N_gNm2','NH4_gNm2','NO3_gNm2','total_P_gPm2'};
target_present = isfield(targets, recognized_targets);
applied_fields = recognized_targets(target_present);
for ii = find(target_present)
    value = targets.(recognized_targets{ii});
    if ~isscalar(value) || ~isfinite(value) || value < 0
        error('Initialize_BG_State_From_Observations:InvalidTargets', ...
            'Observed stock target %s must be a finite nonnegative scalar.', ...
            recognized_targets{ii});
    end
end
if isfield(targets, 'total_N_gNm2') && targets.total_N_gNm2 <= 0
    error('Initialize_BG_State_From_Observations:InvalidTargets', ...
        'A supplied total N target must be positive.');
end
if isfield(targets, 'total_P_gPm2') && targets.total_P_gPm2 <= 0
    error('Initialize_BG_State_From_Observations:InvalidTargets', ...
        'A supplied total P target must be positive.');
end

if isfield(targets, 'POC_gCm2')
    B(9:10) = scaled_pair(B(9:10), targets.POC_gCm2);
end
if isfield(targets, 'MAOC_gCm2')
    B(11) = targets.MAOC_gCm2;
end
if isfield(targets, 'DOC_gCm2')
    B(12:13) = scaled_pair(B(12:13), targets.DOC_gCm2);
end

par = Biogeochemistry_Parameter(1);
microbial_target = isfield(targets, 'MBC_gCm2');
if microbial_target
    B(18:19) = scaled_pair(B(18:19), targets.MBC_gCm2);
    B(27) = B(18)/par.R_CN_bac;
    B(28) = B(19)/par.R_CN_fun;
    B(29) = B(20)/par.R_CN_myc;
    B(30) = B(21)/par.R_CN_myc;
    B(34) = B(22)/par.R_CN_ew;
    B(39) = B(18)/par.R_CP_bac;
    B(40) = B(19)/par.R_CP_fun;
    B(41) = B(20)/par.R_CP_myc;
    B(42) = B(21)/par.R_CP_myc;
end
if isfield(targets, 'NH4_gNm2')
    B(31) = targets.NH4_gNm2;
end
if isfield(targets, 'NO3_gNm2')
    B(32) = targets.NO3_gNm2;
end

N_changed = microbial_target || isfield(targets, 'NH4_gNm2') || ...
    isfield(targets, 'NO3_gNm2') || isfield(targets, 'total_N_gNm2');
if isfield(targets, 'total_N_gNm2')
    total_N_reference = targets.total_N_gNm2;
else
    total_N_reference = initial_total_N;
end
if N_changed
    B(26) = total_N_reference - sum(B([23:25 27:34]));
    if B(26) <= 0
        error('Initialize_BG_State_From_Observations:InconsistentNitrogenStocks', ...
            'Total N is smaller than the explicitly initialised N pools.');
    end
end

organic_fraction = NaN;
inorganic_fractions = nan(1,4);
preserve_organic_P_pools = false;
if isfield(options, 'preserve_organic_P_pools')
    value = options.preserve_organic_P_pools;
    valid_value = (islogical(value) || isnumeric(value)) && ...
        isscalar(value) && isfinite(double(value)) && any(double(value) == [0 1]);
    if ~valid_value
        error('Initialize_BG_State_From_Observations:InvalidOptions', ...
            'The preserve-organic-P-pools option must be a logical scalar.');
    end
    preserve_organic_P_pools = logical(value);
end
%%% Bounds on the primary and secondary inorganic P fractions.  The occluded
%%% pool (B46) accumulates monotonically via K3p*B(45) with no return flux, so
%%% its fraction is not constrained. It is excluded from the convergence
%%% metric in Run_BG_Spinup.m for the same reason.
inorganic_P_fraction_bounds = [0.005 0.90; 0.01 0.95];
if isfield(options, 'inorganic_P_fraction_bounds')
    ipfb = options.inorganic_P_fraction_bounds;
    if ~isnumeric(ipfb) || ~isreal(ipfb) || ~isequal(size(ipfb), [2 2]) || ...
            any(~isfinite(ipfb(:))) || any(ipfb(:,1) < 0) || ...
            any(ipfb(:,2) <= ipfb(:,1))
        error('Initialize_BG_State_From_Observations:InvalidOptions', ...
            'Inorganic P fraction bounds must be a 2-by-2 numeric matrix with positive finite ranges.');
    end
    inorganic_P_fraction_bounds = ipfb;
end
if isfield(targets, 'total_P_gPm2')
    if preserve_organic_P_pools
        organic_P = sum(B([35:42 47]));
        inorganic_P = targets.total_P_gPm2 - organic_P;
        if inorganic_P < 0
            error('Initialize_BG_State_From_Observations:InconsistentPhosphorusStocks', ...
                'Total P is smaller than the preserved organic P pools.');
        end
        inorganic_total = sum(B(43:46));
        if inorganic_P == 0
            B(43:46) = 0;
            inorganic_fractions = zeros(1,4);
        elseif inorganic_total > 0
            inorganic_fractions = B(43:46)/inorganic_total;
            B(43:46) = inorganic_P*inorganic_fractions;
        else
            error('Initialize_BG_State_From_Observations:InconsistentPhosphorusStocks', ...
                'A positive inorganic P target requires a nonzero inorganic P partition.');
        end
        organic_fraction = organic_P/targets.total_P_gPm2;
    else
        [organic_P, organic_fraction, inorganic_fractions] = phosphorus_options( ...
            options, B, targets.total_P_gPm2);
        inorganic_P = targets.total_P_gPm2 - organic_P;
        B(38) = organic_P - sum(B([35:37 39:42 47]));
        B(43:46) = inorganic_P*inorganic_fractions;
    end
    if all(inorganic_fractions > 0)
        pool_names = {'solution','primary','secondary'};
        for k = 2:3
            if inorganic_fractions(k) < inorganic_P_fraction_bounds(k-1,1) || ...
                    inorganic_fractions(k) > inorganic_P_fraction_bounds(k-1,2)
                error('Initialize_BG_State_From_Observations:ImplausibleInorganicPFractions', ...
                    'Inorganic P fraction %s = %.4g is outside the plausible range [%.3g, %.3g].', ...
                    pool_names{k}, inorganic_fractions(k), ...
                    inorganic_P_fraction_bounds(k-1,1), inorganic_P_fraction_bounds(k-1,2));
            end
        end
    end
    total_P_reference = targets.total_P_gPm2;
elseif microbial_target
    B(38) = initial_total_P - sum(B([35:37 39:47]));
    total_P_reference = initial_total_P;
else
    total_P_reference = initial_total_P;
end
P_changed = isfield(targets, 'total_P_gPm2') || microbial_target;
if P_changed && B(38) <= 0
    error('Initialize_BG_State_From_Observations:InconsistentPhosphorusStocks', ...
        'Organic P is smaller than the explicitly initialised organic P pools.');
end

if B(38) > 0
    som_C_to_P = sum(B(9:11))/B(38);
else
    som_C_to_P = NaN;
end
%%% Cropland surface-soil SOM C:P spans roughly 100-250 by mass, so a mapped
%%% state outside that window signals an inconsistent stock combination.
som_C_to_P_bounds = [103 244];
if isfield(options, 'som_C_to_P_bounds')
    som_C_to_P_bounds = options.som_C_to_P_bounds;
    if ~isnumeric(som_C_to_P_bounds) || ~isreal(som_C_to_P_bounds) || ...
            ~isequal(size(som_C_to_P_bounds), [1 2]) || ...
            any(~isfinite(som_C_to_P_bounds)) || som_C_to_P_bounds(1) <= 0 || ...
            som_C_to_P_bounds(2) <= som_C_to_P_bounds(1)
        error('Initialize_BG_State_From_Observations:InvalidOptions', ...
            'SOM C:P bounds must be positive finite lower and upper limits.');
    end
end
%%% The check belongs to the branch that chose the organic/inorganic split from
%%% a total P prior, because only there is the resulting ratio this function's
%%% own construction.  Preserving measured organic P, or conserving an
%%% unmeasured inventory, inherits whatever stoichiometry the state arrived
%%% with and is not this function's to reject.
chose_P_partition = isfield(targets, 'total_P_gPm2') && ~preserve_organic_P_pools;
if chose_P_partition && isfinite(som_C_to_P) && ...
        (som_C_to_P < som_C_to_P_bounds(1) || som_C_to_P > som_C_to_P_bounds(2))
    error('Initialize_BG_State_From_Observations:ImplausibleStoichiometry', ...
        'The initial SOM C:P ratio %.3g is outside the plausible bounds.', som_C_to_P);
end
if any(~isfinite(B)) || any(B < 0)
    error('Initialize_BG_State_From_Observations:InvalidResult', ...
        'The initialised state contains a nonfinite or negative pool.');
end

diag = struct('applied_fields', {applied_fields}, ...
    'total_N_reference_gNm2', total_N_reference, ...
    'total_P_reference_gPm2', total_P_reference, ...
    'organic_P_fraction', organic_fraction, ...
    'inorganic_P_fractions', inorganic_fractions, ...
    'som_C_to_P_bounds', som_C_to_P_bounds, ...
    'som_C_to_P', som_C_to_P);
end

function [organic_P, organic_fraction, inorganic_fractions] = phosphorus_options(options, B, total_P)
%%% Generic defaults preserve the legacy P initialization when no scenario
%%% specific observed fractionation is supplied.
organic_fraction = 0.35;
inorganic_fractions = [0.01 0.55 0.24 0.20];
if isfield(options, 'organic_P_C_to_P')
    organic_P_C_to_P = options.organic_P_C_to_P;
    if ~isnumeric(organic_P_C_to_P) || ~isreal(organic_P_C_to_P) || ...
            ~isscalar(organic_P_C_to_P) || ~isfinite(organic_P_C_to_P) || ...
            organic_P_C_to_P <= 0
        error('Initialize_BG_State_From_Observations:InvalidOptions', ...
            'The organic P C:P ratio must be a positive finite scalar.');
    end
    if isfield(options, 'organic_P_fraction')
        error('Initialize_BG_State_From_Observations:InvalidOptions', ...
            'Specify either organic P C:P or organic P fraction, not both.');
    end
    organic_P = sum(B(9:11))/organic_P_C_to_P + sum(B([35:37 39:42 47]));
    organic_fraction = organic_P/total_P;
elseif isfield(options, 'organic_P_fraction')
    organic_fraction = options.organic_P_fraction;
    organic_P = organic_fraction*total_P;
else
    organic_P = organic_fraction*total_P;
end
if isfield(options, 'inorganic_P_fractions')
    inorganic_fractions = options.inorganic_P_fractions;
end
if ~isscalar(organic_fraction) || ~isfinite(organic_fraction) || ...
        organic_fraction <= 0 || organic_fraction >= 1
    error('Initialize_BG_State_From_Observations:InvalidOptions', ...
        'The organic P fraction must be finite and strictly between zero and one.');
end
if ~isnumeric(inorganic_fractions) || ~isreal(inorganic_fractions) || ...
        ~isequal(size(inorganic_fractions), [1 4]) || ...
        any(~isfinite(inorganic_fractions)) || any(inorganic_fractions < 0) || ...
        abs(sum(inorganic_fractions) - 1) > 1e-12
    error('Initialize_BG_State_From_Observations:InvalidOptions', ...
        'The four inorganic P fractions must be nonnegative and sum to one.');
end
end

function values = scaled_pair(values, target)
if sum(values) > 0
    values = values*(target/sum(values));
else
    values = [0.5 0.5]*target;
end
end
