%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function settings = Run_Settings_From_Workspace(defaults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(defaults)
    defaults = struct();
end
if ~isstruct(defaults)
    error('Run_Settings_From_Workspace:BadDefaults', 'Default run settings must be a struct.');
end
overrides = defaults;

if evalin('caller', 'exist(''RUN_SETTINGS'',''var'') && ~isempty(RUN_SETTINGS)')
    user_settings = evalin('caller', 'RUN_SETTINGS');
    if ~isstruct(user_settings)
        error('Run_Settings_From_Workspace:BadRunSettings', 'RUN_SETTINGS must be a struct.');
    end
    overrides = merge_settings(overrides, user_settings);
end

names = {'RUN_MODE','RUN_CHEMISTRY','NN_case_study_ERW_MOD', ...
    'RUN_DOMAIN','BUFFER_PLANTING','SPINUP_SEED_FILE', ...
    'SAVE_RESULTS','OPT_PARALLEL_GRID', ...
    'lateral_bc','ERW_region_coverage'};
for ii = 1:numel(names)
    has_var = evalin('caller', sprintf('exist(''%s'',''var'') && ~isempty(%s)', names{ii}, names{ii}));
    if has_var
        overrides.(names{ii}) = evalin('caller', names{ii});
    end
end

settings = Resolve_Run_Settings(overrides);
end

function out = merge_settings(out, in)
fields = fieldnames(in);
for ii = 1:numel(fields)
    out.(fields{ii}) = in.(fields{ii});
end
end
