function test_fit_exchange_affinity()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's exchange calibration (v2/calibration/Fit_Exchange_Affinity) and
%% the pinned Al(OH)3 phase it sets, on Hamden's ERW seed inputs
%% (Inputs/Hamden/Calibration/ERW_seed_inputs_classseed_01.mat):
%%   1  bulk_pH_params.Al_pin_logKso: absent, the pin is the amorphous
%%      phase (10.8) as in v0.9; given, the seed holds dissolved Al at it
%%   2  the fit refuses bad options, targets and seed inputs
%%   3  one phase: the uniform shift found holds the target base saturation
%%      on the reference selectivities (logKH = SCEPTER + z*delta)
%%   4  a phase list that cannot bracket the measured exchangeable Al is
%%      refused
%% (test_hamden_field_configuration checks Hamden's fitted constants.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_fit_exchange_affinity ===\n');
test_directory = fileparts(mfilename('fullpath'));
code = fileparts(test_directory);
v2 = fileparts(code);
addpath(code,fullfile(code,'Parameters'),fullfile(v2,'calibration'));
seed_file = fullfile(code,'Inputs','Hamden','Calibration','ERW_seed_inputs_classseed_01.mat');
S = load(seed_file);
chem = readtable(fullfile(code,'Inputs','Hamden','Soils','Representative_soil_chemistry.csv'));
bs = chem.bs_Ca_fraction + chem.bs_Mg_fraction + chem.bs_K_fraction + chem.bs_Na_fraction;

%% --- 1
C = S.ERW_TC_const;
C.bulk_pH_params.Al_solubility_pin = true;
C.bulk_pH_params = rmfield_if(C.bulk_pH_params,'Al_pin_logKso');
[B0,D0] = seed(S,C);
C.bulk_pH_params.Al_pin_logKso = 10.8;
[B1,D1] = seed(S,C);
C.bulk_pH_params.Al_pin_logKso = 8.04;
[B2,D2] = seed(S,C);
assert(isequaln(B0,B1) && isequaln(D0,D1), ...
    'Without Al_pin_logKso the pin is the amorphous phase, as in v0.9.');
assert(D2.dissolved_Al_umol_L < D0.dissolved_Al_umol_L/100 && B2.f_Al(1) < B0.f_Al(1), ...
    'Crystalline gibbsite holds less Al in solution and on the exchanger.');
fprintf('  1 Al pin default amorphous; gibbsite: Al %.3g -> %.3g umol/L\n', ...
    D0.dissolved_Al_umol_L,D2.dissolved_Al_umol_L);

%% --- 2
assert_error(@() Fit_Exchange_Affinity(seed_file,1.2,5.2),'Fit_Exchange_Affinity:BadTarget');
assert_error(@() Fit_Exchange_Affinity(seed_file,0.5,NaN),'Fit_Exchange_Affinity:BadTarget');
assert_error(@() Fit_Exchange_Affinity(seed_file,0.7,5.2,'ExchangeableAl',0.3), ...
    'Fit_Exchange_Affinity:BadTarget');
assert_error(@() Fit_Exchange_Affinity(seed_file,0.5,5.2,'Bogus',1),'Fit_Exchange_Affinity:BadOptions');
assert_error(@() Fit_Exchange_Affinity(seed_file,0.5,5.2,'Alpha'),'Fit_Exchange_Affinity:BadOptions');
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d,'s'));
t = S.t; save(fullfile(d,'partial.mat'),'t');
assert_error(@() Fit_Exchange_Affinity(fullfile(d,'partial.mat'),0.5,5.2), ...
    'Fit_Exchange_Affinity:BadSeedInputs');
fprintf('  2 bad targets, options and seed inputs refused\n');

%% --- 3
ref = [5.9 4.8 10.47 10.786 16.47];
R = Fit_Exchange_Affinity(seed_file,bs,chem.layers_pH,'AlPinLogKso',8.04,'DeltaGrid',-3:1:0);
assert(height(R) == 1 && R.fitted && ~R.phase_fitted && R.Al_pin_logKso == 8.04 && ...
    abs(R.base_saturation - bs) < 1e-6 && ...
    max(abs(R.logKH - (ref + [1 1 2 2 3]*R.delta))) < 1e-12 && R.alpha == 3.4 && ...
    abs(R.f_Ca + R.f_Mg + R.f_K + R.f_Na + R.f_Al + R.f_H - 1) < 1e-9 && ...
    R.readiness == "pass", ...
    'A one-phase fit holds the target base saturation by one uniform shift.');
fprintf('  3 gibbsite: delta %.4f holds base saturation %.4f (Al %.4f)\n', ...
    R.delta,R.base_saturation,R.f_Al);

%% --- 4
assert_error(@() Fit_Exchange_Affinity(seed_file,bs,chem.layers_pH,'AlPinLogKso',8.04, ...
    'DeltaGrid',0,'ExchangeableAl',chem.Al_fraction),'Fit_Exchange_Affinity:NoPhase');
fprintf('  4 an unbracketed exchangeable Al is refused\n');

fprintf('test_fit_exchange_affinity: all checks passed\n');
end

function [B,D] = seed(S,C)
[B,~,W,EC] = ERW_MOD_PARAM_excel(S.t,S.TT_m,C,S.ERW_biogeochem,S.ERW_PARAM_input, ...
    S.opt_ERW_preprocessing);
D = ERW_Handoff_Diagnostic(ERW_Initialize_TC_Cation_Pools(S.B_out,B,EC,S.TT_m.s(1,:)), ...
    B,W,EC,S.cov_moll);
end

function s = rmfield_if(s,name)
if isfield(s,name)
    s = rmfield(s,name);
end
end

function assert_error(f,id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,id),'Expected %s, got %s: %s',id,ME.identifier,ME.message);
    return
end
error('test_fit_exchange_affinity:NoError','Expected %s.',id);
end
