function test_automatic_spinup_port
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The source spin-up gates and stock conversions remain in the v1.0 port.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
addpath(code_directory,fullfile(v2_directory,'spinup'));

row = table(1200,300,2.5,0.12,40, ...
    'VariableNames',{'rsd','Zbio','InitPOC_gC_kg','InitTN_pct','InitNO3_ugN_g'});
targets = BG_Observation_Targets(row);
soil_mass = 1200*0.3;
assert(abs(targets.POC_gCm2-2.5*soil_mass) < 1e-12);
assert(abs(targets.total_N_gNm2-0.12*soil_mass*10) < 1e-12);
assert(abs(targets.NO3_gNm2-40*soil_mass/1000) < 1e-12);

spinup_source = fileread(fullfile(code_directory,'Run_BG_Spinup.m'));
validator_source = fileread(fullfile(code_directory, ...
    'Validate_BG_State_For_Transient.m'));
assert(contains(spinup_source,'while ~p2_converged'));
assert(contains(spinup_source,'while ~p3_converged'));
assert(~contains(spinup_source,'P3_MAX_LOOPS'));
written = regexp(spinup_source,'P3_NUTRIENT_BUDGET_VERS = (\d+);', ...
    'tokens','once');
accepted = regexp(validator_source,'budget_version == (\d+)\)', ...
    'tokens','once');
assert(isequal(written,accepted));

assert_error(@() Run_Automatic_Site_Spinup(struct( ...
    'NN_case_study_ERW_MOD',5,'RUN_MODE','spinup', ...
    'RUN_CHEMISTRY','erw','SAVE_RESULTS',0)), ...
    'Run_Automatic_Site_Spinup:MissingSpinupMeteo');

workbook = fullfile(code_directory,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
MOD = readtable(workbook,'Sheet','MOD_PARAM','VariableNamingRule','preserve');
ERW = readtable(workbook,'Sheet','ERW_MOD','VariableNamingRule','preserve');
MOD = MOD(MOD.Number == 5,:);
ERW = ERW(ERW.Number == 5,:);
settings = struct('NN_case_study_ERW_MOD',5,'RUN_MODE','spinup', ...
    'RUN_CHEMISTRY','erw','SAVE_RESULTS',0, ...
    'MOD_PARAM_override',MOD,'ERW_PARAM_override',ERW);
assert_error(@() Run_Automatic_Site_Spinup(settings), ...
    'Run_Automatic_Site_Spinup:MissingSpinupMeteo');

bad_MOD = removevars(MOD,MOD.Properties.VariableNames{end});
settings.MOD_PARAM_override = bad_MOD;
assert_error(@() Run_Automatic_Site_Spinup(settings), ...
    'Run_Automatic_Site_Spinup:BadParameterOverride');

fprintf('test_automatic_spinup_port: all checks passed\n');
end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end
