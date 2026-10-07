function S = Load_Rain_Chemistry(MOD_PARAM_input, Input_Directory, id_location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The rain-chemistry record a site's parameter row names (rain_chem_file),
%% as the ERW chemistry takes it: Inputs/Rain_Chem_<id_location>.mat, from
%% Extract_NTN_Rain_Chemistry (long-term NADP NTN wet deposition).  The
%% plot-scale spin-up and the field drivers read it here, so every run of a
%% site sees the same record.
%%
%% Without a record: rain pH NaN (the parameter file's rain pH stands), no
%% Na or strong-anion supply, and an empty source_file, which the drivers
%% read as "no record".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
rain_file = rain_chem_file(MOD_PARAM_input);
if isempty(rain_file)
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

function name = rain_chem_file(T)
%% The row's rain_chem_file, '' when the column is absent, empty or 'none'.
name = '';
if ~istable(T) || ~ismember('rain_chem_file', T.Properties.VariableNames)
    return
end
raw = T.rain_chem_file;
if iscell(raw); raw = raw{1}; end
if ~(isstring(raw) || ischar(raw))
    return
end
raw = strtrim(string(raw));
if isempty(raw) || ismissing(raw(1)) || strlength(raw(1)) == 0 || strcmpi(raw(1), "none")
    return
end
name = char(raw(1));
end
