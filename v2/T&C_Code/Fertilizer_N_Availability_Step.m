%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B_IO_day,state,diag] = Fertilizer_N_Availability_Step(B_IO,Date,state,dt_day,Ts_C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Converts applied fertiliser N to bulk-available N before BIOGEO_UNIT.
%%% Released N enters B31/B32.
%%% Optional B_IO fields:
%%%   fert_frac_NO3, fert_frac_NH4, fert_frac_urea                 [-]
%%%   fert_k_NO3_release, fert_k_NH4_release, fert_k_urea_hydrolysis [1/d]
%%%   NaN fert_k_urea_hydrolysis uses BiogeoPar urea temperature response.

if nargin < 3 || isempty(state)
    state = fert_init_state();
else
    state = fert_upgrade_state(state);
end
if nargin < 4 || isempty(dt_day)
    dt_day = 1;
end
if nargin < 5
    Ts_C = [];
end
if ~isscalar(dt_day) || ~isfinite(dt_day) || dt_day <= 0
    error('Fertilizer_N_Availability:TimeStep', ...
        'dt_day must be a positive finite scalar.');
end

if isdatetime(Date)
    Date_vec = datevec(Date);
elseif isscalar(Date)
    Date_vec = datevec(Date);
else
    Date_vec = Date;
end
Date_vec = reshape(Date_vec,1,[]);
jDay = JULIAN_DAY(Date_vec);

raw_N = fert_daily_value(B_IO,'FertN',Date_vec,jDay,'Nrate',0,1);
if ~isscalar(raw_N) || ~isfinite(raw_N) || raw_N < 0
    error('Fertilizer_N_Availability:FertN', ...
        'B_IO.FertN must resolve to a non-negative finite daily scalar.');
end

has_frac = isfield(B_IO,'fert_frac_NO3') || isfield(B_IO,'fert_frac_NH4') || isfield(B_IO,'fert_frac_urea');
legacy_f_NH4 = fert_daily_value(B_IO,'f_fert_NH4',Date_vec,jDay,'',NaN,0);

if has_frac
    f_NO3 = fert_daily_value(B_IO,'fert_frac_NO3',Date_vec,jDay,'',0,0);
    f_NH4 = fert_daily_value(B_IO,'fert_frac_NH4',Date_vec,jDay,'',0,0);
    f_urea = fert_daily_value(B_IO,'fert_frac_urea',Date_vec,jDay,'',0,0);
    f_sum = f_NO3 + f_NH4 + f_urea;

    if f_sum > 0
        if any(~isfinite([f_NO3 f_NH4 f_urea])) || any([f_NO3 f_NH4 f_urea] < 0) || abs(f_sum - 1) > 1e-8
            error('Fertilizer_N_Availability:Fractions', ...
                'fert_frac_NO3 + fert_frac_NH4 + fert_frac_urea must be non-negative and sum to 1.');
        end
        k_NO3 = fert_daily_value(B_IO,'fert_k_NO3_release',Date_vec,jDay,'',Inf,0);
        k_NH4 = fert_daily_value(B_IO,'fert_k_NH4_release',Date_vec,jDay,'',Inf,0);
        k_urea = fert_daily_value(B_IO,'fert_k_urea_hydrolysis',Date_vec,jDay,'',NaN,0);
        if isnumeric(k_urea) && isscalar(k_urea) && isnan(k_urea)
            k_urea = fert_urea_hydrolysis_rate(B_IO,Ts_C);
        end
        fert_check_rate(k_NO3);
        fert_check_rate(k_NH4);
        fert_check_rate(k_urea);

        state = fert_add_cohort(state,'NO3',raw_N*f_NO3,k_NO3);
        state = fert_add_cohort(state,'NH4',raw_N*f_NH4,k_NH4);
        state = fert_add_cohort(state,'urea',raw_N*f_urea,k_urea);
    elseif raw_N > 0
        if ~isfinite(legacy_f_NH4)
            legacy_f_NH4 = 0.75;
        end
        if isfinite(legacy_f_NH4)
            if legacy_f_NH4 < 0 || legacy_f_NH4 > 1
                error('Fertilizer_N_Availability:f_fert_NH4', ...
                    'B_IO.f_fert_NH4 must be in [0,1] when used.');
            end
            state = fert_add_cohort(state,'NO3',raw_N*(1-legacy_f_NH4),Inf);
            state = fert_add_cohort(state,'NH4',raw_N*legacy_f_NH4,Inf);
        end
    end
elseif raw_N > 0
    B_IO_day = fert_empty_day(B_IO,jDay);
    B_IO_day.FertN(jDay) = raw_N;
    if isfinite(legacy_f_NH4)
        if legacy_f_NH4 < 0 || legacy_f_NH4 > 1
            error('Fertilizer_N_Availability:f_fert_NH4', ...
                'B_IO.f_fert_NH4 must be in [0,1] when used.');
        end
        B_IO_day.f_fert_NH4(jDay) = legacy_f_NH4;
    end
    diag = fert_diag(jDay,raw_N,raw_N,raw_N*max(0,1-legacy_f_NH4),raw_N*min(max(legacy_f_NH4,0),1),0,state,B_IO_day.f_fert_NH4(jDay));
    return
end

[state,release_NO3] = fert_release_cohorts(state,'NO3',dt_day);
[state,release_NH4] = fert_release_cohorts(state,'NH4',dt_day);
[state,release_urea] = fert_release_cohorts(state,'urea',dt_day);
%%% Urea release is collapsed to bulk NH4 availability.
%%% Transient urea alkalinity and urea-C to CO2 are not represented here.
release_NH4 = release_NH4 + release_urea;
release_N = release_NO3 + release_NH4;

B_IO_day = fert_empty_day(B_IO,jDay);
B_IO_day.FertN(jDay) = release_N;
if release_N > 0
    B_IO_day.f_fert_NH4(jDay) = release_NH4/release_N;
end

diag = fert_diag(jDay,raw_N,release_N,release_NO3,release_NH4,release_urea,state,B_IO_day.f_fert_NH4(jDay));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = fert_init_state()
state.NO3_mass = [];
state.NO3_k = [];
state.NH4_mass = [];
state.NH4_k = [];
state.urea_mass = [];
state.urea_k = [];
state.pending_NO3 = 0;
state.pending_NH4 = 0;
state.pending_urea = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = fert_upgrade_state(state)
base = fert_init_state();
names = fieldnames(base);
for ii = 1:numel(names)
    if ~isfield(state,names{ii})
        state.(names{ii}) = base.(names{ii});
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B_IO_day = fert_empty_day(B_IO,jDay)
nDay = max(366,jDay);
B_IO_day = B_IO;
B_IO_day.FertN = zeros(1,nDay);
B_IO_day.f_fert_NH4 = NaN(1,nDay);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fert_daily_value(B_IO,field_name,Date_vec,jDay,var_name,default_value,sum_values)
val = default_value;
if ~isfield(B_IO,field_name)
    return
end
x = B_IO.(field_name);
if istimetable(x)
    ix = find(abs(datenum(x.Time) - datenum(Date_vec(1:3))) < 0.01);
    if isempty(ix)
        return
    end
    if ~isempty(var_name) && any(strcmp(x.Properties.VariableNames,var_name))
        vals = x.(var_name)(ix);
    else
        vals = x{ix,1};
    end
    if sum_values
        val = sum(vals(:),'omitnan');
    else
        val = vals(1);
    end
elseif isscalar(x)
    val = x;
elseif jDay <= numel(x)
    val = x(jDay);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fert_check_rate(k)
if ~(isscalar(k) && ((isfinite(k) && k > 0) || isinf(k)))
    error('Fertilizer_N_Availability:Rate', ...
        'fertiliser transfer rates must be positive [1/d] or Inf.');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = fert_urea_hydrolysis_rate(B_IO,Ts_C)
BiogeoPar = fert_biogeo_params(B_IO);
k_ref = fert_scalar_param(BiogeoPar,'urea_hydrolysis_k_ref',1.0);
T_ref = fert_scalar_param(BiogeoPar,'urea_hydrolysis_T_ref',20.0);
Q10 = fert_scalar_param(BiogeoPar,'urea_hydrolysis_Q10',2.0);
if k_ref <= 0 || Q10 <= 0
    error('Fertilizer_N_Availability:Parameter', ...
        'urea hydrolysis k_ref and Q10 must be positive finite scalars.');
end
if isempty(Ts_C)
    Ts_use = T_ref;
elseif isnumeric(Ts_C) && all(isfinite(Ts_C(:)))
    Ts_use = mean(Ts_C(:));
else
    error('Fertilizer_N_Availability:Temperature', ...
        'Ts_C must be finite when deriving the urea hydrolysis rate.');
end
k = k_ref * Q10^((Ts_use - T_ref)/10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BiogeoPar = fert_biogeo_params(B_IO)
BiogeoPar = Biogeochemistry_Parameter(1);
if isfield(B_IO,'BiogeoPar')
    fn = fieldnames(B_IO.BiogeoPar);
    for ii = 1:numel(fn)
        BiogeoPar.(fn{ii}) = B_IO.BiogeoPar.(fn{ii});
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fert_scalar_param(BiogeoPar,name,default_value)
val = default_value;
if isfield(BiogeoPar,name)
    val = BiogeoPar.(name);
end
if ~(isscalar(val) && isfinite(val))
    error('Fertilizer_N_Availability:Parameter', ...
        'BiogeoPar.%s must be a finite scalar.', name);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = fert_add_cohort(state,name,mass,k)
if mass <= 0
    return
end
mass_name = [name '_mass'];
k_name = [name '_k'];
mass_vec = state.(mass_name);
k_vec = state.(k_name);
mass_vec(end+1) = mass;
k_vec(end+1) = k;
state.(mass_name) = mass_vec;
state.(k_name) = k_vec;
state = fert_refresh_pending(state);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state,release] = fert_release_cohorts(state,name,dt_day)
mass_name = [name '_mass'];
k_name = [name '_k'];
mass = state.(mass_name);
k = state.(k_name);
if isempty(mass)
    release = 0;
    state = fert_refresh_pending(state);
    return
end

rel_frac = zeros(size(mass));
is_fast = isinf(k);
rel_frac(is_fast) = 1;
rel_frac(~is_fast) = 1 - exp(-k(~is_fast)*dt_day);
release_vec = mass .* rel_frac;
mass = max(0,mass - release_vec);
keep = mass > 1e-12;
state.(mass_name) = mass(keep);
state.(k_name) = k(keep);
release = sum(release_vec);
state = fert_refresh_pending(state);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = fert_refresh_pending(state)
state.pending_NO3 = sum(state.NO3_mass);
state.pending_NH4 = sum(state.NH4_mass);
state.pending_urea = sum(state.urea_mass);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diag = fert_diag(jDay,raw_N,release_N,release_NO3,release_NH4,release_urea,state,f_NH4)
diag.jDay = jDay;
diag.raw_N = raw_N;
diag.released_N = release_N;
diag.released_NO3 = release_NO3;
diag.released_NH4 = release_NH4;
diag.released_urea_to_NH4 = release_urea;
diag.pending_NO3 = state.pending_NO3;
diag.pending_NH4 = state.pending_NH4;
diag.pending_urea = state.pending_urea;
diag.pending_N = state.pending_NO3 + state.pending_NH4 + state.pending_urea;
diag.f_fert_NH4 = f_NH4;
end
