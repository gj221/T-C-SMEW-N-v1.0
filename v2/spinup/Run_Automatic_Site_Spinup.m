function bg_file = Run_Automatic_Site_Spinup(RUN_SETTINGS)
%% TC-SMEW-N  |  automatic single-cell site spin-up
%  Runs a single-site biogeochemistry simulation (transient or spin-up) for
%  a scenario defined in ERW_MOD_PARAM_parameters.xlsx.
%
%  TYPICAL WORKFLOW
%    Step 1 - spin-up:   set RUN_MODE = 'spinup', run this script.
%                        Saves an equilibrated pool state to TC_Code/Inputs/.
%    Step 2 - transient: set RUN_MODE = 'transient', run again.
%                        Writes results to SIM/<scenario>/.


if nargin < 1 || isempty(RUN_SETTINGS)
    RUN_SETTINGS = struct();
end
if ~isstruct(RUN_SETTINGS) || ~isscalar(RUN_SETTINGS)
    error('Run_Automatic_Site_Spinup:BadSettings', ...
        'RUN_SETTINGS must be a scalar struct.');
end

spinup_directory = fileparts(mfilename('fullpath'));
current_directory = fileparts(spinup_directory);
Directory = fullfile(current_directory,'T&C_Code');
Input_Directory = fullfile(Directory,'Inputs');
Data_Directory = fullfile(current_directory,'Data');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory, Input_Directory, Parameter_Directory);

%% ========================================================================
%% CONFIGURATION  - edit only this section
%% ============-============================================================
RunDefaults = struct( ...
    'NN_case_study_ERW_MOD', 1,   ...  workbook scenario row
    'RUN_MODE', 'spinup',            ...  'transient' | 'spinup'
    'RUN_CHEMISTRY', 'erw',       ...
    'RUN_DOMAIN', 'single',       ...  'single' | 'field_buffer'
    'BUFFER_PLANTING', 'crop',    ...  'crop' | 'native_perennial'  (field_buffer only)
    'SAVE_RESULTS', 1,            ...  1 = save workspace to SIM/, 0 = skip
    'OPT_PARALLEL_GRID', 0        ...  1 = parallel workers
);

%% ========================================================================
%% DO NOT EDIT BELOW THIS LINE  -  implementation only
%% ========================================================================

RunSettings = Run_Settings_From_Workspace(RunDefaults);
RUN_MODE = RunSettings.run_mode;
RUN_CHEMISTRY = RunSettings.run_chemistry;
SPINUP_SEED_FILE = RunSettings.spinup_seed_file;
NN_case_study_ERW_MOD = RunSettings.case_number;
RUN_DOMAIN = RunSettings.run_domain;
BUFFER_PLANTING = RunSettings.buffer_planting;
SAVE_RESULTS = RunSettings.save_results;
OPT_PARALLEL_GRID = RunSettings.opt_parallel_grid;
lateral_bc = RunSettings.lateral_bc;  % options: 'closed_perimeter', 'open_outflow'
ERW_region_coverage = RunSettings.ERW_region_coverage;

OPT_SoilBiogeochemistry = 1;
IS_BG_SPINUP = RunSettings.is_spinup;

RunConfig = Run_Configuration(RUN_CHEMISTRY, RUN_DOMAIN, BUFFER_PLANTING);
OPT_GRID = RunConfig.OPT_GRID;
OPT_pH = RunConfig.OPT_pH;
OPT_ERW = RunConfig.OPT_ERW;
OPT_BUFFER_PERENNIAL = RunConfig.OPT_BUFFER_PERENNIAL;

%% ========================================================================
%% READ SCENARIO PARAMETERS  (from ERW_MOD_PARAM_parameters.xlsx)
%% ========================================================================
MOD_PARAM_file = fullfile(Input_Directory,'ERW_MOD_PARAM_parameters.xlsx');
%% One reader for the workbook (Read_Parameter_Workbook): the field
%% launcher builds its per-cell override rows from the same typing, so
%% the override validation below compares like with like.
[MOD_PARAM_input, ERW_PARAM_input] = Read_Parameter_Workbook(MOD_PARAM_file);

ERW_PARAM_input = ERW_PARAM_input(ERW_PARAM_input.Number == NN_case_study_ERW_MOD, :);
MOD_PARAM_input = MOD_PARAM_input(MOD_PARAM_input.Number == NN_case_study_ERW_MOD, :);
if height(MOD_PARAM_input) ~= 1
    error('MOD_PARAM case %d must resolve to one row; found %d.', NN_case_study_ERW_MOD, height(MOD_PARAM_input));
end
if height(ERW_PARAM_input) ~= 1
    error('ERW_MOD case %d must resolve to one row; found %d.', NN_case_study_ERW_MOD, height(ERW_PARAM_input));
end
if isfield(RUN_SETTINGS,'MOD_PARAM_override') && ~isempty(RUN_SETTINGS.MOD_PARAM_override)
    MOD_PARAM_input = validate_parameter_override(RUN_SETTINGS.MOD_PARAM_override, ...
        MOD_PARAM_input,'MOD_PARAM_override');
end
if isfield(RUN_SETTINGS,'ERW_PARAM_override') && ~isempty(RUN_SETTINGS.ERW_PARAM_override)
    ERW_PARAM_input = validate_parameter_override(RUN_SETTINGS.ERW_PARAM_override, ...
        ERW_PARAM_input,'ERW_PARAM_override');
end
OPT_ROCK = RunConfig.OPT_ROCK && isfinite(ERW_PARAM_input.M_rock_in) && ERW_PARAM_input.M_rock_in > 0;
RunConfig.OPT_ROCK = OPT_ROCK;

domain_area_m2 = table_value(MOD_PARAM_input,'DomainArea_m2',NaN);
if OPT_GRID == 0 && (isnan(domain_area_m2) || domain_area_m2 <= 0)
    domain_area_m2 = 1;
elseif isnan(domain_area_m2) || domain_area_m2 <= 0
    error('MOD_PARAM case %d requires a positive DomainArea_m2 value.', NN_case_study_ERW_MOD);
end
field_side_m = sqrt(domain_area_m2);  % square domain side, main + buffer
buffer_width_m = table_value(MOD_PARAM_input,'BufferWidth_m',NaN);
if OPT_GRID == 1
    if isnan(buffer_width_m) || buffer_width_m <= 0
        error('MOD_PARAM case %d requires a positive BufferWidth_m value for field_buffer runs.', NN_case_study_ERW_MOD);
    end
    if field_side_m <= 2*buffer_width_m
        error('MOD_PARAM case %d requires DomainArea_m2 side length > 2*BufferWidth_m.', NN_case_study_ERW_MOD);
    end
else
    buffer_width_m = NaN;
end

%% The workbook Location names the forcing and state files, as in
%% Run_One_Subfield (Data_<Location>.mat, Data_<Location>_BG.mat).
id_location = char(table_value(MOD_PARAM_input,'id_location',MOD_PARAM_input.Location));
forcing_file = char(table_value(MOD_PARAM_input,'forcing_file',"Data_" + string(id_location) + ".mat"));
daily_meteo_file = char(table_value(MOD_PARAM_input,'meteo_file',""));
pH_dilution_factor = table_value(MOD_PARAM_input,'pH_DilutionFactor',NaN);
if ~isnan(pH_dilution_factor)
    if pH_dilution_factor <= 0
        error('MOD_PARAM case %d has invalid pH_DilutionFactor.', NN_case_study_ERW_MOD);
    end
    bulk_pH_params.n_extract = pH_dilution_factor;
end
transient_start = table_datenum(MOD_PARAM_input,'transient_start',NaN);
transient_end = table_datenum(MOD_PARAM_input,'transient_end',NaN);
spinup_start = table_datenum(MOD_PARAM_input,'meteo_start',NaN);
spinup_end = table_datenum(MOD_PARAM_input,'meteo_end',NaN);
if IS_BG_SPINUP && (isempty(strtrim(daily_meteo_file)) || ...
        ~isfinite(spinup_start) || ~isfinite(spinup_end) || ...
        spinup_end < spinup_start)
    error('Run_Automatic_Site_Spinup:MissingSpinupMeteo', ...
        ['Automatic spin-up requires meteo_file and a valid ' ...
         'meteo_start to meteo_end window in MOD_PARAM.']);
end
scenario_transient = char(table_value(MOD_PARAM_input,'scenario_transient',string(id_location) + "_transient"));

base_scenario = scenario_transient;
if IS_BG_SPINUP
    base_scenario = [base_scenario '_spinup'];
end
OutputNames = Output_Names(base_scenario, id_location, RunConfig, ERW_region_coverage, ERW_PARAM_input);
RUN_DESCRIPTOR = char(OutputNames.descriptor);
SCENARIO_TEST = scenario_path(OutputNames.scenario);
RESULT_FILE_NAME = char(OutputNames.result_file);
Lat = MOD_PARAM_input.Lat;
Lon = MOD_PARAM_input.Lon;
DeltaGMT = round(Lon/15);

fprintf('=== TC-SMEW-N Configuration ===\n');
fprintf('  Case: %d, run = %s, chemistry = %s, domain = %s\n', ...
    NN_case_study_ERW_MOD, char(RUN_MODE), RunConfig.run_chemistry, RunConfig.domain_label);
fprintf('  Flags: pH = %d, ERW = %d, rock = %d\n', OPT_pH, OPT_ERW, OPT_ROCK);
fprintf('  Site: %s, Scenario: %s\n', id_location, strtrim(SCENARIO_TEST));
fprintf('  Parameter file: %s\n', char(MOD_PARAM_input.MOD_PARAM));
fprintf('  Forcing file: %s\n', forcing_file);
if ~isempty(strtrim(daily_meteo_file))
    fprintf('  Daily meteo file: %s\n', daily_meteo_file);
end
fprintf('  Window: %s - %s\n', datestr(transient_start,'yyyy-mm-dd'), datestr(transient_end,'yyyy-mm-dd'));
if OPT_GRID == 1
    fprintf('  Domain area: %.0f m2, BufferWidth_m: %.2f, Parallel grid cells: %d\n', ...
        domain_area_m2, buffer_width_m, OPT_PARALLEL_GRID);
end
if IS_BG_SPINUP
    fprintf('  Spin-up termination: convergence\n');
    if isfinite(spinup_start) && isfinite(spinup_end)
        fprintf('  Spin-up meteo window: %s - %s\n', datestr(spinup_start,'yyyy-mm-dd'), datestr(spinup_end,'yyyy-mm-dd'));
    end
end
fprintf('==============================\n');

%% ========================================================================
%% LOAD INITIAL POOL STATE  (BG state for transient; seed or cold-start for spinup)
%% ========================================================================
bg_file = BG_State_File(Input_Directory, id_location);
if isfield(RUN_SETTINGS,'BG_OUTPUT_FILE') && ...
        ~isempty(strtrim(char(string(RUN_SETTINGS.BG_OUTPUT_FILE))))
    bg_file = char(string(RUN_SETTINGS.BG_OUTPUT_FILE));
end
%% A spin-up without an explicit output writes the location's default state
%% file, the seed the transient campaign reads; it never replaces one unasked.
if IS_BG_SPINUP && ~(isfield(RUN_SETTINGS,'BG_OUTPUT_FILE') && ...
        ~isempty(strtrim(char(string(RUN_SETTINGS.BG_OUTPUT_FILE))))) && isfile(bg_file)
    error('Run_Automatic_Site_Spinup:WouldOverwriteState', ...
        ['A spin-up without BG_OUTPUT_FILE writes %s, which already exists. ' ...
         'Pass BG_OUTPUT_FILE, or move the existing state away first.'], bg_file);
end
if IS_BG_SPINUP
    bg_state = [];
    if exist('SPINUP_SEED_FILE','var') && ~isempty(strtrim(char(SPINUP_SEED_FILE)))
        if isfile(char(SPINUP_SEED_FILE))
            bg_state = load_bg_state(char(SPINUP_SEED_FILE));
            fprintf('BG spin-up seed: %s\n', SPINUP_SEED_FILE);
        else
            error('TC_SMEW_N_run:MissingSeedFile', 'SPINUP_SEED_FILE not found: %s', SPINUP_SEED_FILE);
        end
    else
        fprintf('BG state: spin-up mode starts from the parameter initial state.\n');
    end
    if ~isempty(bg_state)
        B_out      = bg_state.B_out;
        O_out      = bg_state.O_out;
        NavlI_out  = bg_state.NavlI_out;
        RexmyI_out = bg_state.RexmyI_out;
        if ismatrix(bg_state.NupI_L_out)
            NupI_L_out = reshape(bg_state.NupI_L_out, 1, size(bg_state.NupI_L_out,1), size(bg_state.NupI_L_out,2));
        else
            NupI_L_out = bg_state.NupI_L_out;
        end
        if isfield(bg_state,'Nreserve_L_out');  Nreserve_L_out  = bg_state.Nreserve_L_out;  end
        if isfield(bg_state,'Preserve_L_out');   Preserve_L_out   = bg_state.Preserve_L_out;   end
        if isfield(bg_state,'Kreserve_L_out');   Kreserve_L_out   = bg_state.Kreserve_L_out;   end
        if isfield(bg_state,'Careserve_L_out');  Careserve_L_out  = bg_state.Careserve_L_out;  end
        if isfield(bg_state,'Mgreserve_L_out');  Mgreserve_L_out  = bg_state.Mgreserve_L_out;  end
        if isfield(bg_state,'Sireserve_L_out');  Sireserve_L_out  = bg_state.Sireserve_L_out;  end
        if isfield(bg_state,'Nreserve_H_out');   Nreserve_H_out   = bg_state.Nreserve_H_out;   end
        if isfield(bg_state,'Preserve_H_out');   Preserve_H_out   = bg_state.Preserve_H_out;   end
        if isfield(bg_state,'Kreserve_H_out');   Kreserve_H_out   = bg_state.Kreserve_H_out;   end
        if isfield(bg_state,'Careserve_H_out');  Careserve_H_out  = bg_state.Careserve_H_out;  end
        if isfield(bg_state,'Mgreserve_H_out');  Mgreserve_H_out  = bg_state.Mgreserve_H_out;  end
        if isfield(bg_state,'Sireserve_H_out');  Sireserve_H_out  = bg_state.Sireserve_H_out;  end
    end
else
    bg_state = load_bg_state(bg_file);
    Validate_BG_State_For_Transient(bg_state, bg_file);
    B_out = bg_state.B_out;
    NavlI_out = bg_state.NavlI_out;
    NupI_L_out = bg_state.NupI_L_out;
    O_out = bg_state.O_out;
    RexmyI_out = bg_state.RexmyI_out;
    %% Respiration, the Na-seed water loss and the plant reserves start the
    %% way the v0.9 driver starts them, from the run and the parameter file.
    % Keep the converged exchange state for the transient pH inversion.
    if isfield(bg_state, 'spinup_erw_handoff_diagnostic') && ...
            isstruct(bg_state.spinup_erw_handoff_diagnostic)
        D_erw = bg_state.spinup_erw_handoff_diagnostic;
        f_CEC_seed = [];
        if isfield(D_erw, 'f_CEC') && isnumeric(D_erw.f_CEC) && ...
                numel(D_erw.f_CEC) == 6 && ...
                all(isfinite(D_erw.f_CEC)) && all(D_erw.f_CEC >= 0) && ...
                abs(sum(D_erw.f_CEC) - 1) <= 1e-3
            f_CEC_seed = double(D_erw.f_CEC(:).');
        end
        if ~isempty(f_CEC_seed)
            if ~exist('bulk_pH_params', 'var') || ~isstruct(bulk_pH_params)
                bulk_pH_params = struct();
            end
            bulk_pH_params.initial_f_CEC = f_CEC_seed;
            if isfield(D_erw, 'initial_Al_w_model_units') && ...
                    isnumeric(D_erw.initial_Al_w_model_units) && ...
                    isfinite(D_erw.initial_Al_w_model_units) && ...
                    D_erw.initial_Al_w_model_units >= 0
                bulk_pH_params.initial_Al_w = double(D_erw.initial_Al_w_model_units);
            end
        end
    end
    fprintf('BG state: %s\n', bg_file);

end

data_input_file = resolve_input_file(Input_Directory, forcing_file);
use_daily_meteo = ~isempty(strtrim(daily_meteo_file));

%% ========================================================================
%% LOAD METEOROLOGICAL FORCING
%% ========================================================================
NN = MOD_PARAM_input.NumberOfDays * 24;
dt = 3600;
dth = 1;
eval(append('Zs', '=', eval('MOD_PARAM_input.Zs'), ';'));
ms = length(Zs) - 1;
cc = MOD_PARAM_input.cc;

site_Lat = Lat; site_Lon = Lon; site_DeltaGMT = DeltaGMT;
load(data_input_file);
%% The forcing file carries its own site; the case row decides where the
%% run is.  A file prepared for another site is refused rather than
%% silently moving the run (and its solar geometry) there.
if exist('Lat','var') && exist('Lon','var') && ...
        (abs(double(Lat) - site_Lat) > 0.01 || abs(double(Lon) - site_Lon) > 0.01)
    error('Run_Automatic_Site_Spinup:SiteMismatch', ...
        ['%s holds forcing for %.4f N %.4f E; the case row is at ' ...
         '%.4f N %.4f E.'], data_input_file, Lat, Lon, site_Lat, site_Lon);
end
Lat = site_Lat; Lon = site_Lon; DeltaGMT = site_DeltaGMT;
clear site_Lat site_Lon site_DeltaGMT
% Find the observation struct generically (the only workspace struct with Date and pH)
obs_vars = whos;
field_obs = struct();
for iv = 1:numel(obs_vars)
    if strcmp(obs_vars(iv).class, 'struct')
        candidate = eval(obs_vars(iv).name);
        if isstruct(candidate) && isfield(candidate, 'Date') && isfield(candidate, 'pH')
            field_obs = candidate;
            break;
        end
    end
end
if isempty(fieldnames(field_obs)) && IS_BG_SPINUP
    PHs = table_value(MOD_PARAM_input,'PH',NaN);
    if ~(isnumeric(PHs) && isscalar(PHs) && isfinite(PHs) && PHs > 0 && PHs <= 14)
        error('TC_SMEW_N_run:InvalidObservedPH', ...
            'Spin-up requires a finite cell bulk pH in MOD_PARAM.PH.');
    end
    pH_sem = NaN;
elseif isempty(fieldnames(field_obs))
    error('TC_SMEW_N_run:MissingObservedPH', ...
        'The forcing file must contain a struct with Date and pH fields.');
else
    obs_dates = field_obs.Date;
    if ~isdatetime(obs_dates)
        obs_dates = datetime(obs_dates, 'ConvertFrom', 'datenum');
    end
    target_date = dateshift(datetime(transient_start, 'ConvertFrom', 'datenum'), ...
        'start', 'day');
    pH_row = find(dateshift(obs_dates, 'start', 'day') == target_date);
    if numel(pH_row) ~= 1
        error('TC_SMEW_N_run:MissingObservedPH', ...
            'Expected one field-mean pH observation on %s; found %d.', ...
            datestr(target_date, 'yyyy-mm-dd'), numel(pH_row));
    end
    PHs = field_obs.pH(pH_row);
    if ~(isscalar(PHs) && isfinite(PHs) && PHs > 0)
        error('TC_SMEW_N_run:InvalidObservedPH', ...
            'The field-mean pH on the transient start date must be finite and positive.');
    end
    pH_sem = NaN;
    if isfield(field_obs, 'pH_sem')
        pH_sem = field_obs.pH_sem(pH_row);
    end
end
fprintf('  pH target: %.4f (field mean; SEM %.4f)\n', ...
    PHs, pH_sem);
if use_daily_meteo
    daily_meteo_path = resolve_meteo_file(daily_meteo_file,Data_Directory,Input_Directory);
    if IS_BG_SPINUP
        if isnan(spinup_start) || isnan(spinup_end)
            error('MOD_PARAM case %d requires meteo_start/meteo_end for daily-meteo spin-up.', NN_case_study_ERW_MOD);
        end
        if spinup_end < spinup_start
            error('TC_SMEW_N_run:SpinupWindow', 'meteo_end must be greater than or equal to meteo_start.');
        end
        [spinup_all_met, spinup_met_info] = Load_Spinup_Meteo(daily_meteo_path, ...
            spinup_start, spinup_end, Lat, Lon, DeltaGMT);
        [spinup_all_yr, ~, ~] = datevec(spinup_all_met.Date);
        %% The ring and the extension, defined once with the field
        %% spin-up's tracking run (Run_Field_Lateral_Tracking).
        [spinup_cycle_years, handoff_start_dn, handoff_end_dn, spinup_avail_years] = ...
            Spinup_Cycle_Years(spinup_all_met.Date, transient_start);
        n_spinup_years = numel(spinup_cycle_years);
        fprintf('Spin-up forcing: %d complete cycle years (%d-%d), repeated until convergence\n', ...
            numel(spinup_cycle_years), spinup_cycle_years(1), spinup_cycle_years(end));
        if isfinite(handoff_end_dn)
            fprintf('  Handoff extension: %s - %s\n', ...
                datestr(handoff_start_dn,'yyyy-mm-dd HH:MM'), ...
                datestr(handoff_end_dn - 1/24,'yyyy-mm-dd HH:MM'));
        end
        fprintf('  Meteo source: %s\n', spinup_met_info.source);
    else
        met_start = transient_start;
        met_end = transient_end;
        if isnan(met_start) || isnan(met_end)
            error('MOD_PARAM case %d requires transient_start/transient_end for daily meteo.', NN_case_study_ERW_MOD);
        end
        [daily_met, daily_info] = Load_Spinup_Meteo(daily_meteo_path, met_start, met_end, Lat, Lon, DeltaGMT);
        Date = daily_met.Date;    Pr = daily_met.Pr;      Ta = daily_met.Ta;
        Ws = daily_met.Ws;        ea = daily_met.ea;      Pre = daily_met.Pre;
        SAD1 = daily_met.SAD1;    SAD2 = daily_met.SAD2;
        SAB1 = daily_met.SAB1;    SAB2 = daily_met.SAB2;
        N = daily_met.N;          Tdew = daily_met.Tdew;  esat = daily_met.esat;
        PARB = daily_met.PARB;    PARD = daily_met.PARD;
        NN = numel(Date);
        if NN ~= MOD_PARAM_input.NumberOfDays * 24
            error('Daily meteo transient forcing has %d hours, expected %d.', NN, MOD_PARAM_input.NumberOfDays * 24);
        end
        fprintf('Forcing: %d h (%d d), daily meteo, Date(1)=%s, Date(end)=%s\n', ...
            NN, NN/24, datestr(Date(1),'yyyy-mm-dd'), datestr(Date(end),'yyyy-mm-dd'));
        fprintf('  Meteo source: %s (%d daily records)\n', daily_info.source, daily_info.n_days);
    end
else
    % Slice the configured transient window from the latest donor year and
    % year-shift to the scenario year.
    if isnan(transient_start) || isnan(transient_end)
        error('Transient case %d requires transient_start and transient_end in MOD_PARAM.', NN_case_study_ERW_MOD);
    end
    [YS,MS,DS] = datevec(transient_start);
    [YE,ME,DE] = datevec(transient_end);
    [YY,~,~] = datevec(Date(:));
    slice_idx = []; donor_year = NaN;
    for cy = sort(unique(YY), 'descend').'
        d_start = datenum(cy,MS,DS);
        d_end = datenum(cy,ME,DE) + 1;
        if d_end <= d_start
            d_end = datenum(cy+1,ME,DE) + 1;
        end
        mask = (Date(:) >= d_start) & (Date(:) < d_end);
        if sum(mask) == NN
            donor_year = cy; slice_idx = find(mask); break
        end
    end
    if isempty(slice_idx)
        error('No year in scaffold has a full %d-hour configured transient window', NN);
    end
    Date = Shift_Donor_Dates_To_Year(Date(slice_idx), YS); %% donor -> scenario year
    Pr   = Pr(slice_idx);    Ta   = Ta(slice_idx);
    Ws   = Ws(slice_idx);    ea   = ea(slice_idx);    Pre  = Pre(slice_idx);
    SAD1 = SAD1(slice_idx);  SAD2 = SAD2(slice_idx);
    SAB1 = SAB1(slice_idx);  SAB2 = SAB2(slice_idx);
    N    = N(slice_idx);     Tdew = Tdew(slice_idx);  esat = esat(slice_idx);
    PARB = PARB(slice_idx);  PARD = PARD(slice_idx);
    fprintf('Forcing: %d h (%d d), donor %d -> %d, Date(1)=%s, Date(end)=%s\n', ...
        NN, NN/24, donor_year, YS, datestr(Date(1),'yyyy-mm-dd'), datestr(Date(end),'yyyy-mm-dd'));
    clear YS MS DS YE ME DE
end

if ~IS_BG_SPINUP
    Ca = Load_CO2_Record(Date, Input_Directory);
end
ERW_rain_chemistry = Load_Rain_Chemistry(MOD_PARAM_input, Input_Directory, id_location);

%% ========================================================================
%% DERIVE METEOROLOGICAL VARIABLES
%% ========================================================================
t_bef = -0.67; t_aft = 1.67;
Oa = 210000;
Ws(Ws <= 0) = 0.01;
if ~IS_BG_SPINUP
    Ds = esat - ea;
    Ds(Ds < 0) = 0;
    Datam = zeros(NN, 4);
    [YE, MO, DA, HO, MI, SE] = datevec(Date);
    Datam(:, 1:4) = [YE, MO, DA, HO];
    clear YE MO DA HO MI SE
end

%% Rock amendment from scenario
if OPT_ERW == 0
    fprintf('ERW chemistry: off\n');
elseif OPT_ROCK
    fprintf('ERW chemistry: on with workbook rock addition\n');
else
    fprintf('ERW chemistry: on, no rock addition\n');
end

%% ========================================================================
%% RUN MODEL  (transient; spinup cycle loop is in Run_BG_Spinup.m)
%% ========================================================================
cd(Directory);

if OPT_GRID == 1
    topo = Build_Cell_Topology(field_side_m, buffer_width_m, 4, 4);
    if OPT_PARALLEL_GRID == 1
        if isempty(ver('parallel'))
            error('Parallel Computing Toolbox is required when OPT_PARALLEL_GRID = 1.');
        end
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end
        fprintf('Parallel grid-cell routing: %d workers\n', pool.NumWorkers);
    end
else
    topo = [];
end

if IS_BG_SPINUP && isfield(RUN_SETTINGS,'CELL_ID')
    %% Run_BG_Spinup stamps these on the state before its single write.
    RUN_SETTINGS.CELL_ID = validate_cell_id(RUN_SETTINGS.CELL_ID);
    if ~isfield(RUN_SETTINGS,'CELL_AGRONOMIC_IDENTITY') || ...
            isempty(strtrim(char(string(RUN_SETTINGS.CELL_AGRONOMIC_IDENTITY))))
        error('Run_Automatic_Site_Spinup:MissingCellIdentity', ...
            'CELL_AGRONOMIC_IDENTITY is required when CELL_ID is supplied.');
    end
end
if IS_BG_SPINUP
    run('Run_BG_Spinup.m');
else
    % Transient: define parameter paths and run MAIN_FRAME for the full window
    PARAM_IC = fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM), '.m']);
    ERW_PARAM_IC = fullfile(Parameter_Directory, 'ERW_MOD_PARAM_excel.m');
    MAIN_FRAME_BG_ERW_10mins;
    cd(current_directory);
end


% PH_out is pre-allocated as zeros; row 1 (initial state) is never written
% by the model loop, so set it explicitly from the parameter initial value.
if exist('PH_out','var') && exist('PHs','var') && PH_out(1) == 0
    PH_out(1) = PHs;
end

% Save the daily fertilizer schedule actually seen by BIOGEO_UNIT.
FertN_daily = NaN;
FertP_daily = NaN;
FertK_daily = NaN;
if exist('B_IO','var') && exist('jDay','var') && ~istimetable(B_IO.FertN)
    jj = jDay(:);
    jj = jj(jj >= 1 & jj <= numel(B_IO.FertN));
    FertN_daily = B_IO.FertN(jj(:));
    FertP_daily = B_IO.FertP(jj(:));
    FertK_daily = B_IO.FertK(jj(:));
end

%% ========================================================================
%% END-OF-RUN DIAGNOSTICS
%% ========================================================================
if exist('topo','var') && isstruct(topo) && ~isempty(topo)
    domain_area_m2 = sum(topo.cell_area_m2);
end
domain_area_ha = domain_area_m2 / 1e4;
n_years = NN / (365.25 * 24);
n_days  = NN / 24;

if exist('exit_acc', 'var') && isstruct(exit_acc) && isfield(exit_acc, 'NO3_leach_ts_gN_m2') && any(exit_acc.NO3_leach_ts_gN_m2 > 0)
    NO3_leach_kg_ha_yr = mean(exit_acc.NO3_leach_ts_gN_m2, 'omitnan') * 365.25 * 10;
    NO3_source = 'exit_acc';
elseif exist('LEAK_NO3', 'var')
    NO3_leach_kg_ha_yr = mean(LEAK_NO3(:), 'omitnan') * 365.25 * 10;
    NO3_source = 'LEAK_NO3';
else
    NO3_leach_kg_ha_yr = NaN;
    NO3_source = 'unavailable';
end

N2O_kg_ha_yr = 0;
if exist('N2Oflx', 'var') && ~isempty(N2Oflx)
    N2O_kg_ha_yr = mean(N2Oflx(:), 'omitnan') * 365.25 * 10;
end

% Actual mass over the simulated window = annualized rate x window length in years
NO3_leach_kg_ha_window = NO3_leach_kg_ha_yr * n_years;
N2O_kg_ha_window       = N2O_kg_ha_yr       * n_years;

% The summary only describes a run that simulated the configured window.  In
% spin-up mode it does not: MAIN_FRAME's last pass is the short ERW
% charge-balance diagnostic (so the fertilizer, NO3 and N2O series cover a
% couple of days), while NN still describes the transient window that was never
% simulated.  Every line below would therefore be a few days of state reported
% against a months-long window.  Run_BG_Spinup prints its own convergence,
% handoff and budget diagnostics, which are the meaningful spin-up output.
if ~IS_BG_SPINUP
    fprintf('\n=== End-of-Run Summary ===\n');
    fprintf('  Site: %s, Scenario: %s\n', id_location, strtrim(SCENARIO_TEST));
    fprintf('  Domain: %.0f m2 (%.2f ha)\n', domain_area_m2, domain_area_ha);
    fprintf('  Duration: %.2f years (%.0f days)\n', n_years, n_days);
    if ~isscalar(FertN_daily) || ~isnan(FertN_daily)
        fprintf('  Fertilizer N applied: %.2f g N m-2 (%.1f kg N ha-1)\n', ...
            sum(FertN_daily, 'omitnan'), sum(FertN_daily, 'omitnan') * 10);
    end
    fprintf('  NO3 leaching: %.2f kg N ha-1 (window) | %.2f kg N ha-1 yr-1 (annualized) [%s]\n', ...
        NO3_leach_kg_ha_window, NO3_leach_kg_ha_yr, NO3_source);
    fprintf('  N2O emission: %.4f kg N2O-N ha-1 (window) | %.4f kg N2O-N ha-1 yr-1 (annualized)\n', ...
        N2O_kg_ha_window, N2O_kg_ha_yr);
    if n_years < 0.95
        fprintf('  NOTE: annualized rates extrapolate a %.0f-day (%.2f yr) window x%.1f; seasonal bias likely\n', ...
            n_days, n_years, 1/n_years);
    end
    fprintf('================================\n');
end

%% ========================================================================
%% SAVE RESULTS TO SIM/
%% ========================================================================
if SAVE_RESULTS == 1
    fl_save = fullfile(current_directory,'SIM',SCENARIO_TEST);
    if ~exist(fl_save,'dir')
        mkdir(fl_save);
    end
    save_name = fullfile(fl_save, RESULT_FILE_NAME);
    save(save_name);
    fprintf('Saved %s result: %s\n', char(RUN_MODE), save_name);
end

end

function val = table_value(T,name,default)
% Read a scalar from a table column, returning default if the column is absent or empty.
val = default;
if ~ismember(name,T.Properties.VariableNames); return; end
raw = T.(name);
if istable(raw); raw = table2array(raw); end
if iscell(raw); raw = raw{1}; end
if isstring(raw)
    if isempty(raw) || ismissing(raw(1)) || strlength(strtrim(raw(1))) == 0 || strcmpi(strtrim(raw(1)), "none"); return; end
    val = raw(1);
elseif ischar(raw)
    if isempty(strtrim(raw)) || strcmpi(strtrim(raw), 'none'); return; end
    val = raw;
elseif isnumeric(raw)
    if isempty(raw) || (isscalar(raw) && isnan(raw)); return; end
    val = raw(1);
elseif isdatetime(raw)
    if isempty(raw) || isnat(raw(1)); return; end
    val = raw(1);
end
end

function override = validate_parameter_override(override,base,label)
%% The override must carry every workbook column; it may carry more --
%% the optional columns table_value reads with a default (meteo_file,
%% meteo_start, meteo_end, transient_start, transient_end, forcing_file,
%% id_location) are exactly what a
%% caller supplies this way when the workbook row does not.
if ~istable(override) || height(override) ~= 1 || ...
        ~all(ismember(base.Properties.VariableNames,override.Properties.VariableNames))
    error('Run_Automatic_Site_Spinup:BadParameterOverride', ...
        ['%s must be a one-row table carrying every workbook parameter ' ...
         'column (extra optional columns are allowed).'],label);
end
end

function cell_id = validate_cell_id(value)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value == fix(value) && value >= 1)
    error('Run_Automatic_Site_Spinup:BadCellId', ...
        'CELL_ID must be a positive integer.');
end
cell_id = double(value);
end

function S = Load_Rain_Chemistry(MOD_PARAM_input, Input_Directory, id_location)
S = struct('rain_pH_vwm', NaN, ...
    'rain_strong_anion_flux_umolc_m2_d', 0, ...
    'rain_strong_anion_ueq_L', 0, ...
    'rain_strong_anion_deposition_umolc_m2', 0, ...
    'rain_Ca_flux_g_m2_d', 0, ...
    'rain_Mg_flux_g_m2_d', 0, ...
    'rain_K_flux_g_m2_d', 0, ...
    'rain_Na_flux_g_m2_d', 0, ...
    'rain_Na_flux_umol_m2_d', 0, ...
    'source_file', '');
rain_file = char(table_value(MOD_PARAM_input, 'rain_chem_file', ""));
if isempty(strtrim(rain_file))
    return
end
rain_mat = fullfile(Input_Directory, ['Rain_Chem_' char(string(id_location)) '.mat']);
if ~isfile(rain_mat)
    error('TC_SMEW_N_run:MissingRainChemistry', ...
        'Rain chemistry file %s is configured, but %s is missing. Run Extract_NTN_Rain_Chemistry first.', ...
        rain_file, rain_mat);
end
loaded = load(rain_mat, 'Rain_Chem');
if ~isfield(loaded, 'Rain_Chem')
    error('TC_SMEW_N_run:InvalidRainChemistry', '%s does not contain Rain_Chem.', rain_mat);
end
RC = loaded.Rain_Chem;
required = {'pH_vwm','strong_anion_ueq_L','strong_anion_deposition_umolc_m2','strong_anion_flux_umolc_m2_d',...
    'Ca_flux_g_m2_d','Mg_flux_g_m2_d','K_flux_g_m2_d','Na_flux_g_m2_d','Na_flux_umol_m2_d'};
for ii = 1:numel(required)
    if ~isfield(RC, required{ii}) || ~isfinite(RC.(required{ii}))
        error('TC_SMEW_N_run:InvalidRainChemistry', ...
            'Rain_Chem.%s must be finite in %s.', required{ii}, rain_mat);
    end
end
if RC.pH_vwm <= 0
    error('TC_SMEW_N_run:InvalidRainChemistry', ...
        'Rain_Chem.pH_vwm must be positive in %s.', rain_mat);
end
S.rain_pH_vwm = RC.pH_vwm;
S.rain_strong_anion_flux_umolc_m2_d = RC.strong_anion_flux_umolc_m2_d;
S.rain_strong_anion_ueq_L = RC.strong_anion_ueq_L;
S.rain_strong_anion_deposition_umolc_m2 = RC.strong_anion_deposition_umolc_m2;
S.rain_Ca_flux_g_m2_d = RC.Ca_flux_g_m2_d;
S.rain_Mg_flux_g_m2_d = RC.Mg_flux_g_m2_d;
S.rain_K_flux_g_m2_d = RC.K_flux_g_m2_d;
S.rain_Na_flux_g_m2_d = RC.Na_flux_g_m2_d;
S.rain_Na_flux_umol_m2_d = RC.Na_flux_umol_m2_d;
S.source_file = rain_file;
end

function dn = table_datenum(T,name,default)
% Read a date column from a table and convert it to a MATLAB datenum.
dn = default;
val = table_value(T,name,"");
if isdatetime(val)
    dn = datenum(val);
elseif isnumeric(val)
    dn = val;
elseif isstring(val) || ischar(val)
    txt = char(val);
    if isempty(strtrim(txt)); return; end
    try
        dn = datenum(datetime(txt,'InputFormat','yyyy-MM-dd'));
    catch
        dn = datenum(datetime(txt));
    end
end
end

function path_out = scenario_path(name)
% Strip leading path separators from a scenario name for use as a folder path.
name = char(string(name));
path_out = regexprep(strtrim(name),'^[\\/]+','');
end

function file_out = resolve_input_file(Input_Directory,file_name)
% Resolve a forcing file name to a full path, checking Input_Directory as a fallback.
file_name = char(file_name);
if isfile(file_name)
    file_out = file_name;
else
    file_out = fullfile(Input_Directory,file_name);
end
if ~isfile(file_out)
    error('Input file not found: %s', file_out);
end
end

function path_out = resolve_meteo_file(file_name,Data_Directory,Input_Directory)
%% A spin-up meteo file is an absolute path, or a name resolved first in
%% the Data directory (daily station files) and then in Inputs (the
%% hourly forcing .mat files Load_Spinup_Meteo accepts directly).
file_name = char(file_name);
candidates = {file_name, fullfile(Data_Directory,file_name), ...
    fullfile(Input_Directory,file_name)};
hit = find(cellfun(@isfile,candidates),1);
if isempty(hit)
    error('Run_Automatic_Site_Spinup:MeteoFileNotFound', ...
        'Spin-up meteo file %s not found (tried %s and %s).', ...
        file_name,candidates{2},candidates{3});
end
path_out = candidates{hit};
end
