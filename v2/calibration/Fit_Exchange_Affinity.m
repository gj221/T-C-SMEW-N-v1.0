function result = Fit_Exchange_Affinity(seed_inputs_file, base_saturation, bulk_pH, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's Appelo exchange constants, fitted the way the Energy Farm's
%% were (runtime/Base_1D/Campaign_Parameters): the reference logKH
%% (SCEPTER's defaults) shifted by z*delta, one shift per unit charge so the
%% relative selectivities stay the reference's, at a chosen alpha, with
%% delta set so that the site's measured base saturation holds at the pH it
%% was measured at.
%%
%% The base saturation fitted is the one that emerges at the handoff: the
%% seed keeps the site's spun-up soil solution and puts the exchanger in
%% equilibrium with it (spinup_equilibrium) at the pore-water pH that
%% reproduces the bulk pH.  The pinned Al(OH)3 phase sets the dissolved Al,
%% which competes for the sites the bases need, so the shift is fitted for
%% each phase in AlPinLogKso; a phase too soluble for the site cannot reach
%% the target.  Given the measured exchangeable Al, the phase's log *Kso is
%% fitted as well, between the listed phases, so that the exchanger holds
%% both the measured bases and the measured Al.
%%
%% The fractions are of the exchanger's fixed capacity, which is a buffered
%% CEC (at pH 7, or by summation with buffer acidity, as the Energy Farm's
%% 13 cmolc/kg): the model carries the pH-dependent charge as exchangeable H.
%%
%% seed_inputs_file  the ERW seed's inputs at a handoff, as the plot-scale
%%                   driver saves them when asked (ERW_SEED_INPUTS_FILE in
%%                   Run_Automatic_Site_Spinup's settings)
%% base_saturation   the measured (Ca + Mg + K + Na)/CEC
%% bulk_pH           the bulk pH that exchanger was measured at (the row's
%%                   extraction ratio holds)
%% Options:
%%   'Alpha'           3.4, Appelo (1994) and SCEPTER's default, as the
%%                     Energy Farm campaign adopted
%%   'ReferenceLogKH'  SCEPTER's [5.9 4.8 10.47 10.786 16.47] (Na K Ca Mg Al)
%%   'AlPinLogKso'     [8.04 9.66 10.8]: crystalline gibbsite,
%%                     microcrystalline gibbsite, amorphous Al(OH)3
%%   'ExchangeableAl'  [] or the measured KCl-extractable Al / CEC, to fit
%%                     the phase too
%%   'CEC'             [] keeps the row's CEC_tot_in [cmolc/kg]
%%   'DeltaGrid'       -8:0.5:2, the shifts scanned for a bracket
%% result  a table, one row per listed phase and, given ExchangeableAl, a
%%   last row for the fitted phase (phase_fitted true): Al_pin_logKso,
%%   fitted (false when the target lies beyond the scanned shifts), delta,
%%   logKH, alpha, the emergent exchanger (base_saturation, f_Ca, f_Mg, f_K,
%%   f_Na, f_Al, f_H), pore_pH, dissolved Ca_umol_L and Al_umol_L, the
%%   residual strong anion An_umolc_L and the handoff readiness.  The fitted
%%   logKH and alpha go in the site's bulk_pH_params.appelo, the phase in
%%   bulk_pH_params.Al_pin_logKso.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('Alpha',3.4,'ReferenceLogKH',[5.9 4.8 10.47 10.786 16.47], ...
    'AlPinLogKso',[8.04 9.66 10.8],'ExchangeableAl',[],'CEC',[],'DeltaGrid',-8:0.5:2);
if mod(numel(varargin),2) ~= 0
    error('Fit_Exchange_Affinity:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k})
        error('Fit_Exchange_Affinity:BadOptions','Unknown option %s.',varargin{k});
    end
    opt.(varargin{k}) = varargin{k+1};
end
if ~(isscalar(base_saturation) && base_saturation > 0 && base_saturation < 1)
    error('Fit_Exchange_Affinity:BadTarget','The base saturation must lie in (0, 1).');
end
if ~(isscalar(bulk_pH) && isfinite(bulk_pH) && bulk_pH > 0 && bulk_pH < 14)
    error('Fit_Exchange_Affinity:BadTarget','The bulk pH must be finite and within (0, 14).');
end
if ~isempty(opt.ExchangeableAl) && ~(isscalar(opt.ExchangeableAl) && ...
        opt.ExchangeableAl > 0 && base_saturation + opt.ExchangeableAl < 1)
    error('Fit_Exchange_Affinity:BadTarget', ...
        'The exchangeable Al must be positive and leave sites for H beside the bases.');
end
code = fullfile(fileparts(fileparts(mfilename('fullpath'))),'T&C_Code');
addpath(code,fullfile(code,'Parameters'));
S = load(seed_inputs_file);
needed = {'t','TT_m','ERW_TC_const','ERW_biogeochem','ERW_PARAM_input', ...
    'opt_ERW_preprocessing','B_out','cov_moll'};
if ~all(isfield(S,needed))
    error('Fit_Exchange_Affinity:BadSeedInputs','%s lacks %s.',seed_inputs_file, ...
        strjoin(needed(~isfield(S,needed)),', '));
end
E = S.ERW_PARAM_input;
E.fX_conc_convention = "spinup_equilibrium";
if ~isempty(opt.CEC)
    E.CEC_tot_in = opt.CEC;
end
phases = sort(opt.AlPinLogKso(:))';
n = numel(phases);
rows = cell(n,1);
for p = 1:n
    rows{p} = fit_shift(S,E,bulk_pH,base_saturation,phases(p),opt);
end
result = vertcat(rows{:});
result.phase_fitted = false(n,1);
if isempty(opt.ExchangeableAl)
    return
end
%% With the base saturation held, the exchangeable Al rises with the phase's
%% solubility (about tenfold per unit of log *Kso, as dissolved Al3+ does)
%% until the bases can no longer be held: bracket the measured fraction
%% between a listed phase that holds the bases with less Al and the next,
%% which holds more or cannot hold the bases, then solve.
excess = @(row) excess_Al(row,opt.ExchangeableAl);
e = arrayfun(@(p) excess(result(p,:)),1:n);
k = find(e(1:end-1) < 0 & e(2:end) >= 0,1);
if isempty(k)
    error('Fit_Exchange_Affinity:NoPhase', ...
        ['No two adjacent phases in %s bracket the exchangeable Al %.4g with the base ' ...
        'saturation held (f_Al %s, fitted %s).'],mat2str(phases),opt.ExchangeableAl, ...
        mat2str(result.f_Al',4),mat2str(result.fitted'));
end
logKso = fzero(@(x) excess(fit_shift(S,E,bulk_pH,base_saturation,x,opt)), ...
    phases([k k+1]),optimset('TolX',1e-5));
row = fit_shift(S,E,bulk_pH,base_saturation,logKso,opt);
if ~row.fitted || abs(row.f_Al/opt.ExchangeableAl - 1) > 1e-3
    error('Fit_Exchange_Affinity:NoPhase', ...
        'The phases that hold the base saturation give at most %.4g exchangeable Al, below %.4g.', ...
        row.f_Al,opt.ExchangeableAl);
end
row.phase_fitted = true;
result = [result; row];
end

function e = excess_Al(row,target)
%% log10 of the exchangeable Al over the measured; +1 where the phase is too
%% soluble to hold the base saturation.
e = 1;
if row.fitted
    e = log10(row.f_Al/target);
end
end

function row = fit_shift(S,E,bulk_pH,target,logKso,opt)
%% The uniform shift that holds the target base saturation with one phase.
z = [1 1 2 2 3];
C = S.ERW_TC_const;
C.PH_in = bulk_pH;
C.bulk_pH_params.Al_solubility_pin = true;
C.bulk_pH_params.Al_pin_logKso = logKso;
state = @(d) seed_state(S,C,E,opt.ReferenceLogKH+z*d,opt.Alpha);
bs = arrayfun(@(d) state(d).base_saturation,opt.DeltaGrid);
k = find((bs(1:end-1)-target).*(bs(2:end)-target) <= 0,1);
if isempty(k)
    [~,best] = min(abs(bs-target));
    d = opt.DeltaGrid(best);
    fitted = false;
else
    d = fzero(@(x) state(x).base_saturation-target,opt.DeltaGrid([k k+1]), ...
        optimset('TolX',1e-6));
    fitted = true;
end
st = state(d);
row = table(logKso,fitted,d,opt.ReferenceLogKH+z*d,opt.Alpha,st.base_saturation, ...
    st.f(1),st.f(2),st.f(3),st.f(4),st.f(5),st.f(6),st.pore_pH,st.Ca_umol_L, ...
    st.Al_umol_L,st.An_umolc_L,st.readiness,'VariableNames',{'Al_pin_logKso', ...
    'fitted','delta','logKH','alpha','base_saturation','f_Ca','f_Mg','f_K','f_Na', ...
    'f_Al','f_H','pore_pH','Ca_umol_L','Al_umol_L','An_umolc_L','readiness'});
end

function st = seed_state(S,C,E,logKH,alpha)
%% The emergent seed and its handoff readiness under trial constants.
C.bulk_pH_params.appelo = struct('logKH',logKH,'alpha',alpha);
st = struct('base_saturation',NaN,'f',NaN(1,6),'pore_pH',NaN,'Ca_umol_L',NaN, ...
    'Al_umol_L',NaN,'An_umolc_L',NaN,'readiness',"error");
try
    [B,~,W,EC] = ERW_MOD_PARAM_excel(S.t,S.TT_m,C,S.ERW_biogeochem,E,S.opt_ERW_preprocessing);
    P1 = ERW_Initialize_TC_Cation_Pools(S.B_out,B,EC,S.TT_m.s(1,:));
    D = ERW_Handoff_Diagnostic(P1,B,W,EC,S.cov_moll);
    R = ERW_Handoff_Readiness(D);
    f = [B.f_Ca(1) B.f_Mg(1) B.f_K(1) B.f_Na(1) B.f_Al(1) B.f_H(1)];
    st = struct('base_saturation',sum(f(1:4)),'f',f,'pore_pH',B.pH(1), ...
        'Ca_umol_L',B.Ca(1),'Al_umol_L',D.dissolved_Al_umol_L,'An_umolc_L',B.An(1), ...
        'readiness',string(R.status));
catch ME
    %% a seed the trial constants cannot build is no fit; a missing
    %% function is a fault, not a trial
    if any(strcmp(ME.identifier,{'MATLAB:UndefinedFunction','MATLAB:undefinedVarOrFunction'}))
        rethrow(ME);
    end
    st.readiness = "error: " + string(ME.message);
end
end
