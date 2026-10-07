function Deposition = Build_Site_Deposition(ntn_file, tdep_dir, forcing_file, lat, lon, years, out_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Site_Deposition
%%
%% A site's event deposition record, for the biogeochemistry's
%% B_IO.Dep_events (Deposition_Day): N, K, Ca and Mg arriving wet with each
%% day's precipitation and dry every day.
%%
%% Wet, per calendar year: the precipitation-weighted mean concentration of
%% the NADP NTN samples collected within the year
%% (Extract_NTN_Rain_Chemistry), which the model multiplies by each day's
%% precipitation.  A year counts when its samples span at least 90 % of it
%% and its valid samples caught at least 75 % of its precipitation (NADP's
%% completeness criteria); any other year takes the pooled concentration of
%% the three nearest complete years.
%% Dry, per calendar year: TDep's dry deposition of the 4 km cell holding
%% the site (Extract_TDep_Dry_Deposition), spread over the year's days; a
%% year without all its grids in tdep_dir takes the mean of the three
%% nearest years that have them.
%% The forcing (a T&C .mat with Date and Pr [mm/h], hourly) gives each
%% year's precipitation, for the summary and the long-term mean rates.
%%
%% Deposition.events     year, wet_g_m2_per_mm, dry_g_m2_d over [N K Ca Mg]
%%                       (B_IO.Dep_events)
%% Deposition.wet_mg_L   per year, as element (N, S) or ion, with its source
%% Deposition.dry_kg_ha_yr  per year, with its source
%% Deposition.annual     per complete forcing year: precipitation, wet, dry
%%                       and total deposition [kg ha-1 yr-1]
%% Deposition.mean_g_m2_d  [N K Ca Mg] over those years (B_IO.DepN ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
years = double(years(:));
n_years = numel(years);
elements = {'N','S','Ca','Mg','K','Na','Cl'};
model = {'N','K','Ca','Mg'};
[~,model_col] = ismember(model,elements);
mw = struct('NH4',18.038,'NO3',62.004,'SO4',96.06,'N',14.007,'S',32.06);

%% --- wet: the year's precipitation-weighted mean concentration -----------
wet = nan(n_years,numel(elements));
valid_mm = zeros(n_years,1);
complete = false(n_years,1);
span = zeros(n_years,1);
valid_fraction = zeros(n_years,1);
for y = 1:n_years
    start_dn = datenum(years(y),1,1);
    end_dn = datenum(years(y)+1,1,1);
    try
        R = Extract_NTN_Rain_Chemistry(ntn_file,'',start_dn,end_dn);
    catch ME
        if any(strcmp(ME.identifier,{'Extract_NTN_Rain_Chemistry:EmptyWindow', ...
                'Extract_NTN_Rain_Chemistry:NoValidSample'}))
            continue
        end
        rethrow(ME)
    end
    days = R.window(2)-R.window(1);
    span(y) = days/(end_dn-start_dn);
    valid_fraction(y) = R.valid_precip_fraction;
    complete(y) = span(y) >= 0.9 && valid_fraction(y) >= 0.75;
    c = R.pwm_mg_L;
    wet(y,:) = [c.NH4/mw.NH4*mw.N + c.NO3/mw.NO3*mw.N, c.SO4/mw.SO4*mw.S, ...
        c.Ca, c.Mg, c.K, c.Na, c.Cl];
    valid_mm(y) = R.valid_precip_fraction*R.gauge_precip_mm_yr*days/365.25;
end
if ~any(complete)
    error('Build_Site_Deposition:NoCompleteYear', ...
        'No year of %s meets the completeness criteria.',ntn_file);
end
wet_source = strings(n_years,1);
for y = 1:n_years
    if complete(y)
        wet_source(y) = sprintf('NTN %d',years(y));
        continue
    end
    near = nearest_three(years(y),years(complete));
    k = ismember(years,near);
    wet(y,:) = sum(wet(k,:).*valid_mm(k),1)/sum(valid_mm(k));
    wet_source(y) = sprintf('NTN %s pooled',strjoin(string(near.'),'+'));
end

%% --- dry: TDep's annual dry deposition of the site's cell ----------------
tdep_species = {'n','s','ca','mg','k','na','cl'};
has_grids = arrayfun(@(y) all(cellfun(@(s) isfile(fullfile(tdep_dir, ...
    sprintf('%s_dw-%d.zip',s,y))),tdep_species)),years);
if ~any(has_grids)
    error('Build_Site_Deposition:NoDryGrids','%s holds no complete year of TDep grids.',tdep_dir);
end
T = Extract_TDep_Dry_Deposition(tdep_dir,lat,lon,years(has_grids),tdep_species);
dry = nan(n_years,numel(elements));
dry(has_grids,:) = T{:,strcat(tdep_species,'_dw_kg_ha')};
dry_source = strings(n_years,1);
dry_source(has_grids) = "TDep " + string(years(has_grids));
for y = find(~has_grids).'
    near = nearest_three(years(y),years(has_grids));
    dry(y,:) = mean(dry(ismember(years,near),:),1);
    dry_source(y) = sprintf('TDep %s mean',strjoin(string(near.'),'+'));
end

%% --- the model's coefficients and the forcing years' totals ------------
days_in_year = datenum(years+1,1,1)-datenum(years,1,1);
events.year = years;
events.wet_g_m2_per_mm = wet(:,model_col)*1e-3;           % mg/L x L/m2/mm
events.dry_g_m2_d = dry(:,model_col)*0.1./days_in_year;   % kg/ha -> g/m2
events.elements = model;

F = load(forcing_file,'Date','Pr');
forcing_year = year(datetime(F.Date,'ConvertFrom','datenum'));
annual_years = years(arrayfun(@(y) nnz(forcing_year == y) >= 24*(datenum(y+1,1,1)-datenum(y,1,1)),years));
precip_mm = arrayfun(@(y) sum(F.Pr(forcing_year == y)),annual_years);
[~,a] = ismember(annual_years,years);
wet_kg = wet(a,:).*precip_mm*1e-2;      % mg/L x mm -> mg/m2 -> kg/ha
total = wet_kg+dry(a,:);
annual = array2table([annual_years precip_mm wet_kg dry(a,:) total], ...
    'VariableNames',[{'year','precip_mm'} strcat('wet_',elements,'_kg_ha') ...
    strcat('dry_',elements,'_kg_ha') strcat('total_',elements,'_kg_ha')]);
mean_g_m2_d = sum(total(:,model_col),1)*0.1/sum(days_in_year(a));

Deposition = struct();
Deposition.site_lat_lon = [lat lon];
Deposition.elements = elements;
Deposition.events = events;
Deposition.wet_mg_L = array2table(wet,'VariableNames',elements);
Deposition.wet_mg_L = addvars(Deposition.wet_mg_L,years,wet_source,complete,span, ...
    valid_fraction,'Before',1,'NewVariableNames', ...
    {'year','source','complete','sample_span','valid_precip_fraction'});
Deposition.dry_kg_ha_yr = array2table(dry,'VariableNames',elements);
Deposition.dry_kg_ha_yr = addvars(Deposition.dry_kg_ha_yr,years,dry_source, ...
    'Before',1,'NewVariableNames',{'year','source'});
Deposition.annual = annual;
Deposition.mean_g_m2_d = mean_g_m2_d;
Deposition.tdep_cell = T.Properties.UserData;
Deposition.sources = struct('ntn_file',char(ntn_file),'tdep_dir',char(tdep_dir), ...
    'tdep_version','TDep MMF 2025.01','forcing_file',char(forcing_file));
Deposition.method = ['wet: the year''s NTN precipitation-weighted mean concentration x ' ...
    'each day''s precipitation; dry: the year''s TDep dry deposition of the site''s ' ...
    '4 km cell, spread evenly over its days (Deposition_Day)'];
if nargin >= 7 && ~isempty(out_file)
    save(out_file,'Deposition','-v7');
end
end

function near = nearest_three(y,candidates)
%% The three candidate years nearest y, the earlier first on a tie.
candidates = candidates(:);
[~,order] = sortrows([abs(candidates-y) candidates]);
near = sort(candidates(order(1:min(3,numel(order)))));
end
