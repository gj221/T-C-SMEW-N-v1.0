function Rain_Chem = Extract_NTN_Rain_Chemistry(ntn_file, out_file, start_dn, end_dn, ...
    site_precip_mm_yr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract_NTN_Rain_Chemistry
%%
%% Long-term mean wet deposition from an NADP National Trends Network
%% weekly site file (NTN-<site>-w-s-mg.csv: concentrations in mg/L, -9 for
%% missing, precipitation in mm), as the Rain_Chem record
%% Load_Rain_Chemistry reads (Inputs/Rain_Chem_<id_location>.mat).
%%
%% NADP's convention: the precipitation-weighted mean (PWM) concentration
%% of the valid samples (validity code 'w', every ion and the pH present),
%% times all the precipitation the gauge measured in the window, valid
%% sample or not.  Values flagged '<' (below detection) are used as
%% reported.  Rain pH is the pH of the PWM hydrogen-ion concentration.
%%
%% The strong-acid anions are sulphate and chloride.  Nitrate is not among
%% them: its N reaches the model as N deposition (with ammonium), so
%% counting it again in the anion budget would deposit it twice.
%%
%% start_dn, end_dn (optional) restrict the window to whole samples whose
%% collection began and ended inside it.
%%
%% site_precip_mm_yr (optional): the site's own mean precipitation over the
%% same window [mm/yr].  The concentrations are regional and the amount is
%% local -- NADP's deposition maps take the concentration from the network
%% and the precipitation from the place -- so the fluxes become the PWM
%% concentrations times this precipitation, the rain the model receives,
%% instead of the gauge's.  Both precipitations are recorded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = detectImportOptions(ntn_file, 'Delimiter', ',');
opts = setvartype(opts, intersect({'dateOn','dateOff','valcode','invalcode','siteID'}, ...
    opts.VariableNames), 'char');
T = readtable(ntn_file, opts);
ions = {'Ca','Mg','K','Na','NH4','NO3','Cl','SO4'};
need = [ions, {'ph','ppt','dateOn','dateOff','valcode'}];
missing = need(~ismember(need, T.Properties.VariableNames));
if ~isempty(missing)
    error('Extract_NTN_Rain_Chemistry:MissingColumns', ...
        '%s lacks the NTN column(s) %s.', ntn_file, strjoin(missing, ', '));
end
on = datenum(T.dateOn, 'yyyy-mm-dd HH:MM');
off = datenum(T.dateOff, 'yyyy-mm-dd HH:MM');
if nargin < 3 || isempty(start_dn); start_dn = min(on); end
if nargin < 4 || isempty(end_dn); end_dn = max(off); end
in = on >= start_dn & off <= end_dn;
T = T(in,:); on = on(in); off = off(in);
if isempty(T)
    error('Extract_NTN_Rain_Chemistry:EmptyWindow', 'No NTN sample lies in the window.');
end
ppt = T.ppt; ppt(~(ppt >= 0)) = 0;                  % gauge precipitation [mm]
conc = zeros(height(T), numel(ions));
for k = 1:numel(ions); conc(:,k) = T.(ions{k}); end
valid = startsWith(strtrim(string(T.valcode)), "w") & all(conc >= 0, 2) & ...
    T.ph > 0 & ppt > 0;
if ~any(valid)
    error('Extract_NTN_Rain_Chemistry:NoValidSample', 'No valid NTN sample in the window.');
end
w = ppt(valid);
pwm = sum(conc(valid,:).*w, 1) ./ sum(w);           % [mg/L]
H_pwm = sum(10.^(-T.ph(valid)).*w) / sum(w);        % [mol/L]
days = max(off) - min(on);
P_gauge = sum(ppt);                                 % [mm] = [L/m2]
P_total = P_gauge;
amount = 'all gauge precipitation';
if nargin >= 5 && ~isempty(site_precip_mm_yr)
    if ~(isnumeric(site_precip_mm_yr) && isscalar(site_precip_mm_yr) && ...
            isfinite(site_precip_mm_yr) && site_precip_mm_yr > 0)
        error('Extract_NTN_Rain_Chemistry:BadSitePrecipitation', ...
            'site_precip_mm_yr must be a finite positive scalar [mm/yr].');
    end
    P_total = site_precip_mm_yr * days / 365.25;    % the site's, over the window
    amount = sprintf('the site''s %.1f mm/yr', site_precip_mm_yr);
end

mw = struct('Ca',40.078,'Mg',24.305,'K',39.098,'Na',22.990, ...
    'NH4',18.038,'NO3',62.004,'Cl',35.453,'SO4',96.06);
z = struct('Ca',2,'Mg',2,'K',1,'Na',1,'NH4',1,'NO3',1,'Cl',1,'SO4',2);
c = cell2struct(num2cell(pwm), ions, 2);
ueq = @(ion) c.(ion) / mw.(ion) * z.(ion) * 1000;   % [ueq/L]
flux_g = @(ion) c.(ion) * P_total / days / 1000;    % [g ion m-2 d-1]

Rain_Chem = struct();
Rain_Chem.site = char(string(T.siteID(1)));
Rain_Chem.source_file = char(ntn_file);
Rain_Chem.window = [min(on) max(off)];
Rain_Chem.n_samples = height(T);
Rain_Chem.n_valid = nnz(valid);
Rain_Chem.precip_mm_yr = P_total / days * 365.25;
Rain_Chem.gauge_precip_mm_yr = P_gauge / days * 365.25;
Rain_Chem.valid_precip_fraction = sum(w) / P_gauge;
Rain_Chem.pwm_mg_L = c;
Rain_Chem.pH_vwm = -log10(H_pwm);
Rain_Chem.strong_anion_ueq_L = ueq('SO4') + ueq('Cl');
Rain_Chem.strong_anion_deposition_umolc_m2 = Rain_Chem.strong_anion_ueq_L * P_total;
Rain_Chem.strong_anion_flux_umolc_m2_d = Rain_Chem.strong_anion_deposition_umolc_m2 / days;
Rain_Chem.Ca_flux_g_m2_d = flux_g('Ca');
Rain_Chem.Mg_flux_g_m2_d = flux_g('Mg');
Rain_Chem.K_flux_g_m2_d = flux_g('K');
Rain_Chem.Na_flux_g_m2_d = flux_g('Na');
Rain_Chem.Na_flux_umol_m2_d = c.Na / mw.Na * 1000 * P_total / days;
Rain_Chem.NH4_N_flux_g_m2_d = flux_g('NH4') * 14.007 / mw.NH4;
Rain_Chem.NO3_N_flux_g_m2_d = flux_g('NO3') * 14.007 / mw.NO3;
Rain_Chem.N_flux_g_m2_d = Rain_Chem.NH4_N_flux_g_m2_d + Rain_Chem.NO3_N_flux_g_m2_d;
Rain_Chem.S_flux_g_m2_d = flux_g('SO4') * 32.06 / mw.SO4;
Rain_Chem.method = ['NADP precipitation-weighted mean of valid samples x ' ...
    amount '; strong anions = SO4 + Cl (NO3 enters as N deposition)'];
if nargin >= 2 && ~isempty(out_file)
    save(out_file, 'Rain_Chem', '-v7');
end
end
