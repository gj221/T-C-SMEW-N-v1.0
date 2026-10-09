function [files,settings] = Run_Plot_CEC_Scenarios(cell_inputs,state_file,cec_values,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One plot-scale column's transient for each CEC scenario, every run from
%% the same spun-up state, so that the runs differ in the exchanger's
%% capacity alone.
%%
%% The spin-up does not depend on the CEC: the ERW chemistry is off in every
%% spin-up stage, and the biogeochemistry's own exchange partitions by mass.
%% Each scenario keeps the state's handoff exchanger fractions at its own
%% capacity (CEC_tot_in), and the seed re-partitions them at the pore-water
%% pH that reproduces the bulk pH.  On Hamden's converged class seed that
%% comes within 1.4e-4 of a fresh spinup_equilibrium seed at the scenario's
%% CEC.  The scenarios therefore start at one base saturation, with
%% exchangeable bases in proportion to their CEC.
%%
%% cell_inputs  one cell's configured inputs, as Run_Field_Cell_Spinup
%%              ('DryRun',true) returns them: MOD_PARAM, ERW_PARAM,
%%              agronomic_identity, cell_id
%% state_file   its spun-up state (the cell's own, or its class seed's),
%%              handed off with the site's current exchange constants
%% cec_values   [cmolc/kg] one run per value
%% Name/value
%%   'Case'    the workbook case the site borrows (default 5)
%%   'NDays'   a shorter window from the row's transient start (default the
%%             row's NumberOfDays)
%%   'Only'    indices of cec_values to run here (default all), so that the
%%             scenarios can run in separate processes
%%   'DryRun'  true returns the settings without running
%% files     each scenario's result file, under the code tree's
%%           SIM/<scenario_transient>_CEC_<value>_<descriptor>/ ('' until run)
%% settings  each scenario's Run_Automatic_Site_Spinup settings
%% Summarize_CEC_Scenarios reads the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('Case',5,'NDays',[],'Only',[],'DryRun',false);
if mod(numel(varargin),2) ~= 0
    error('Run_Plot_CEC_Scenarios:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k})
        error('Run_Plot_CEC_Scenarios:BadOptions','Unknown option %s.',varargin{k});
    end
    opt.(varargin{k}) = varargin{k+1};
end
if isempty(cec_values) || any(~isfinite(cec_values) | cec_values <= 0)
    error('Run_Plot_CEC_Scenarios:BadCEC','Every scenario CEC must be finite and positive.');
end
if ~opt.DryRun && ~isfile(state_file)
    error('Run_Plot_CEC_Scenarios:MissingState','Spun-up state not found: %s',state_file);
end
only = opt.Only;
if isempty(only)
    only = 1:numel(cec_values);
end
v2 = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(v2,'spinup'),fullfile(v2,'T&C_Code'),fullfile(v2,'T&C_Code','Parameters'), ...
    fullfile(v2,'T&C_Code','Inputs'));
M0 = cell_inputs.MOD_PARAM;
base = string(M0.Location) + "_transient";
if ismember('scenario_transient',M0.Properties.VariableNames) && ...
        strlength(strtrim(string(M0.scenario_transient))) > 0
    base = string(M0.scenario_transient);
end
files = repmat({''},1,numel(cec_values));
settings = cell(1,numel(cec_values));
for k = only
    tag = "CEC_" + strrep(sprintf('%g',cec_values(k)),'.','p');
    M = M0;
    M.scenario_transient = base + "_" + tag;
    if ~isempty(opt.NDays)
        M.NumberOfDays = opt.NDays;
        d0 = datetime(char(string(M0.transient_start)),'InputFormat','yyyy-MM-dd');
        M.transient_end = string(datestr(d0 + days(opt.NDays - 1),'yyyy-mm-dd'));
    end
    E = cell_inputs.ERW_PARAM;
    E.CEC_tot_in = cec_values(k);
    settings{k} = struct('NN_case_study_ERW_MOD',opt.Case,'RUN_MODE','transient', ...
        'RUN_CHEMISTRY','erw','RUN_DOMAIN','single','SAVE_RESULTS',1, ...
        'MOD_PARAM_override',M,'ERW_PARAM_override',E,'BG_OUTPUT_FILE',state_file, ...
        'CELL_ID',cell_inputs.cell_id, ...
        'CELL_AGRONOMIC_IDENTITY',string(cell_inputs.agronomic_identity) + ";" + tag);
    if opt.DryRun
        continue
    end
    Run_Automatic_Site_Spinup(settings{k});
    found = dir(fullfile(v2,'SIM',char(M.scenario_transient + "_*"),'TC_ERW_results_*.mat'));
    if isempty(found)
        error('Run_Plot_CEC_Scenarios:NoResult','No result saved for %s.',M.scenario_transient);
    end
    [~,newest] = max([found.datenum]);
    files{k} = fullfile(found(newest).folder,found(newest).name);
end
end
