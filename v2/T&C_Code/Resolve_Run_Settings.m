%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function settings = Resolve_Run_Settings(overrides)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(overrides)
    overrides = struct();
end
if ~isstruct(overrides)
    error('Resolve_Run_Settings:BadInput', 'Run settings overrides must be a struct.');
end

settings = struct();
settings.run_mode = normalize_run_mode(string_setting(overrides, {'RUN_MODE','run_mode'}, 'transient'));
settings.is_spinup = settings.run_mode == "spinup";

raw_chemistry = string_setting(overrides, {'RUN_CHEMISTRY','run_chemistry'}, '');
settings.run_chemistry = normalize_chemistry(raw_chemistry);
if ~settings.is_spinup && settings.run_chemistry == "off"
    error('Resolve_Run_Settings:ChemistryOffTransient', ...
        'RUN_CHEMISTRY=''off'' is not permitted for transient runs.');
end
settings.case_number = integer_setting(overrides, {'NN_case_study_ERW_MOD','case_number'}, 1, ...
    'Resolve_Run_Settings:BadCaseNumber', 'NN_case_study_ERW_MOD');
settings.run_domain = normalize_label(string_setting(overrides, {'RUN_DOMAIN','run_domain'}, 'single'));
settings.buffer_planting = normalize_label(string_setting(overrides, {'BUFFER_PLANTING','buffer_planting'}, 'crop'));
settings.spinup_seed_file = char(string_setting(overrides, {'SPINUP_SEED_FILE','spinup_seed_file'}, ''));
settings.save_results = flag_setting(overrides, {'SAVE_RESULTS','save_results'}, 1, 'SAVE_RESULTS');
settings.opt_parallel_grid = flag_setting(overrides, {'OPT_PARALLEL_GRID','opt_parallel_grid'}, 0, 'OPT_PARALLEL_GRID');
settings.lateral_bc = char(string_setting(overrides, {'lateral_bc','LATERAL_BC'}, 'closed_perimeter'));
settings.ERW_region_coverage = struct_setting(overrides, {'ERW_region_coverage','ERW_REGION_COVERAGE'}, ...
    struct('apply_main', true, 'apply_buffer', true));
end

function mode = normalize_run_mode(raw)
mode = lower(strrep(strtrim(string(raw)), "-", "_"));
switch mode
    case "transient"
        mode = "transient";
    case {"spinup","bg_spinup"}
        mode = "spinup";
    otherwise
        error('Resolve_Run_Settings:BadRunMode', 'RUN_MODE must be transient or spinup.');
end
end

function chemistry = normalize_chemistry(raw)
chemistry = lower(strrep(strtrim(string(raw)), "-", "_"));
switch chemistry
    case "off"
        chemistry = "off";
    case "ph_only"
        chemistry = "pH_only";
    case "erw"
        chemistry = "erw";
    otherwise
        error('Resolve_Run_Settings:BadChemistry', ...
            'RUN_CHEMISTRY must be one of: off, pH_only, erw (got ''%s'').', char(raw));
end
end

function label = normalize_label(raw)
label = lower(strrep(strtrim(string(raw)), "-", "_"));
end

function val = integer_setting(S, names, default, err_id, label)
val = numeric_setting(S, names, default, label);
if ~isfinite(val) || val < 1 || abs(val - round(val)) > eps(max(abs(val),1))
    error(err_id, '%s must be a positive integer.', label);
end
val = round(val);
end

function val = numeric_setting(S, names, default, label)
[found, raw] = get_setting(S, names);
if ~found || isempty(raw)
    val = default;
    return
end
if iscell(raw)
    raw = raw{1};
end
if isstring(raw) || ischar(raw)
    raw = str2double(strtrim(char(raw)));
end
if ~(isnumeric(raw) || islogical(raw)) || ~isscalar(raw) || ~isfinite(double(raw))
    error('Resolve_Run_Settings:BadNumericSetting', '%s must be a finite scalar.', label);
end
val = double(raw);
end

function val = flag_setting(S, names, default, label)
[found, raw] = get_setting(S, names);
if ~found || isempty(raw)
    val = default;
    return
end
if iscell(raw)
    raw = raw{1};
end
if isstring(raw) || ischar(raw)
    txt = lower(strtrim(char(raw)));
    if any(strcmp(txt, {'1','true','yes','on'}))
        val = 1;
        return
    elseif any(strcmp(txt, {'0','false','no','off'}))
        val = 0;
        return
    end
end
if ~(isnumeric(raw) || islogical(raw)) || ~isscalar(raw) || ~isfinite(double(raw))
    error('Resolve_Run_Settings:BadFlagSetting', '%s must be 0 or 1.', label);
end
val = double(raw ~= 0);
end

function val = string_setting(S, names, default)
[found, raw] = get_setting(S, names);
if ~found || isempty(raw)
    val = string(default);
    return
end
if iscell(raw)
    raw = raw{1};
end
if ~(isstring(raw) || ischar(raw) || isnumeric(raw) || islogical(raw))
    error('Resolve_Run_Settings:BadStringSetting', 'Run setting must be string-compatible.');
end
val = string(raw);
val = val(1);
end

function val = struct_setting(S, names, default)
[found, raw] = get_setting(S, names);
if ~found || isempty(raw)
    val = default;
    return
end
if ~isstruct(raw)
    error('Resolve_Run_Settings:BadStructSetting', 'Run setting must be a struct.');
end
val = raw;
end

function [found, raw] = get_setting(S, names)
found = false;
raw = [];
fields = fieldnames(S);
for ii = 1:numel(names)
    jj = find(strcmpi(fields, names{ii}), 1);
    if ~isempty(jj)
        found = true;
        raw = S.(fields{jj});
        return
    end
end
end
