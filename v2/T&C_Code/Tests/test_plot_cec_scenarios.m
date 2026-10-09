function test_plot_cec_scenarios()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's CEC scenarios (v2/run/Run_Plot_CEC_Scenarios), configured on
%% Hamden's class-seed cell without running:
%%   1  the scenario bounds: the representative layers' CEC at pH 7 spans
%%      6.4-13.7 cmolc/kg around their mean, 10.8
%%   2  each scenario runs the cell's own rows from the one spun-up state,
%%      with its own CEC, its own scenario name (its own results folder) and
%%      the scenario in its identity, over the row's transient window
%%   3  a shorter window, and one scenario per process
%%   4  a non-positive CEC and a missing state are refused
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_plot_cec_scenarios ===\n');
code = fileparts(fileparts(mfilename('fullpath')));
v2 = fileparts(code);
addpath(code,fullfile(code,'Parameters'),fullfile(v2,'spinup'),fullfile(v2,'run'), ...
    fullfile(v2,'interpolation'));
grid_dir = fullfile(code,'Inputs','Field_Spinup','Hamden');
shp = fullfile(v2,'Spatial inputs','Hamden Blue Trail','Hamden_Blue_Trail_Parcel.shp');
agro = fullfile(grid_dir,'Cell_Agronomic_Hamden_Blue_Trail_Parcel_polygon_01.csv');
tmp = tempname; mkdir(tmp); cleanup = onCleanup(@() rmdir(tmp,'s'));

%% --- 1
chem = readtable(fullfile(code,'Inputs','Hamden','Soils','Representative_soil_chemistry.csv'));
assert(chem.CEC_min_cmolc_kg == 6.4 && chem.CEC_max_cmolc_kg == 13.7 && ...
    chem.CEC_min_cmolc_kg < chem.CEC_cmolc_kg && chem.CEC_cmolc_kg < chem.CEC_max_cmolc_kg, ...
    'The scenario bounds are the representative layers'' CEC range.');
cec = [chem.CEC_min_cmolc_kg chem.CEC_cmolc_kg chem.CEC_max_cmolc_kg];
fprintf('  1 scenarios %s cmolc/kg\n',mat2str(cec,4));

%% --- 2
M = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_MOD_columns.csv'));
E = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_ERW_columns.csv'));
[~,CI] = Run_Field_Cell_Spinup(shp,agro,'GridCellSizeM',5,'MinCoverage','majority', ...
    'Method','per_cell','MODParamColumns',M,'ERWParamColumns',E, ...
    'OutputFile',fullfile(tmp,'CellBGStates.mat'),'DryRun',true);
seed_cell = CI([CI.cell_id] == 1);
state = fullfile(tmp,'classseed_01.mat');
[files,S] = Run_Plot_CEC_Scenarios(seed_cell,state,cec,'DryRun',true);
names = cellfun(@(s) string(s.MOD_PARAM_override.scenario_transient),S);
assert(isequal(names,["Hamden_transient_CEC_6p4" "Hamden_transient_CEC_10p8" "Hamden_transient_CEC_13p7"]) && ...
    all(cellfun(@(s,c) s.ERW_PARAM_override.CEC_tot_in == c,S,num2cell(cec))) && ...
    all(cellfun(@(s) strcmp(s.RUN_MODE,'transient') && strcmp(s.RUN_CHEMISTRY,'erw') && ...
    strcmp(s.BG_OUTPUT_FILE,state) && s.MOD_PARAM_override.NumberOfDays == 2464,S)) && ...
    all(cellfun(@(s,c) endsWith(s.CELL_AGRONOMIC_IDENTITY,"CEC_" + strrep(sprintf('%g',c),'.','p')), ...
    S,num2cell(cec))) && all(cellfun(@isempty,files)), ...
    'Each scenario runs the cell from the one state with its own CEC, name and identity.');
same = @(field) all(cellfun(@(s) isequaln(rmvars(s.(field)),rmvars(S{1}.(field))),S));
assert(same('MOD_PARAM_override') && same('ERW_PARAM_override'), ...
    'The scenarios differ in nothing but the CEC and the name.');
fprintf('  2 three runs from one state, differing in CEC alone\n');

%% --- 3
[~,S3] = Run_Plot_CEC_Scenarios(seed_cell,state,cec,'DryRun',true,'NDays',3,'Only',2);
assert(isempty(S3{1}) && isempty(S3{3}) && S3{2}.MOD_PARAM_override.NumberOfDays == 3 && ...
    S3{2}.MOD_PARAM_override.transient_end == "2020-01-03", ...
    'A shorter window, and one scenario per process.');
fprintf('  3 shorter window; one scenario per process\n');

%% --- 4
assert_error(@() Run_Plot_CEC_Scenarios(seed_cell,state,[6.4 0],'DryRun',true), ...
    'Run_Plot_CEC_Scenarios:BadCEC');
assert_error(@() Run_Plot_CEC_Scenarios(seed_cell,state,6.4),'Run_Plot_CEC_Scenarios:MissingState');
fprintf('  4 bad CEC and missing state refused\n');
fprintf('test_plot_cec_scenarios: all checks passed\n');
end

function T = rmvars(T)
%% A row without what the scenarios set.
T = removevars(T,intersect({'CEC_tot_in','scenario_transient'},T.Properties.VariableNames));
end

function assert_error(f,id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,id),'Expected %s, got %s: %s',id,ME.identifier,ME.message);
    return
end
error('test_plot_cec_scenarios:NoError','Expected %s.',id);
end
