function [met, info] = Load_Spinup_Meteo(meteo_file, start_dn, end_dn, Lat, Lon, DeltaGMT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load_Spinup_Meteo
%%
%% Meteorological forcing for the spin-up phases, from either
%%   * an hourly T&C forcing .mat (Date, Pr, Ta, Ws, ea, Pre, SAD1, SAD2,
%%     SAB1, SAB2, N, Tdew, esat, PARB, PARD -- the variables the
%%     Kantola forcing and the Preprocess_Bruce_Field_ERA5 outputs already
%%     carry), sliced to whole days of the requested window.  The
%%     measured / reanalysis diurnal cycle is used as it is; nothing is
%%     disaggregated and nothing is re-synthesised.
%%   * a daily station file (Load_Daily_Meteo), disaggregated to hourly.
%% Both return the same struct, so the spin-up phases need not know which
%% they were given.  start_dn / end_dn = NaN take the file's first and
%% last complete calendar years, as Load_Daily_Meteo does.
%%
%% Lat, Lon and DeltaGMT are the site's.  An hourly file carries its own,
%% and they must agree: the radiation partition stored in the file was
%% done for the file's location and clock, so a mismatch is reported, not
%% silently accepted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meteo_file = char(meteo_file);
[~,~,ext] = fileparts(meteo_file);
if ~strcmpi(ext,'.mat')
    [met, info] = Load_Daily_Meteo(meteo_file, start_dn, end_dn, Lat, Lon, DeltaGMT);
    info.kind = 'daily station file, disaggregated to hourly';
    return
end
if ~isfile(meteo_file)
    error('Load_Spinup_Meteo:MissingFile', ...
        'Spin-up forcing file not found: %s', meteo_file);
end
if isdatetime(start_dn); start_dn = datenum(start_dn); end %#ok<DATNM>
if isdatetime(end_dn);   end_dn   = datenum(end_dn);   end %#ok<DATNM>

fields = {'Date','Pr','Ta','Ws','ea','Pre','SAD1','SAD2','SAB1','SAB2', ...
    'N','Tdew','esat','PARB','PARD'};
S = load(meteo_file);
missing = fields(~isfield(S,fields));
if ~isempty(missing)
    error('Load_Spinup_Meteo:MissingVariables', ...
        '%s lacks the forcing variable(s) %s.', meteo_file, strjoin(missing,', '));
end
one_second = 1/(24*3600);
Date = S.Date(:);
if numel(Date) < 2 || any(abs(diff(Date) - 1/24) > one_second)
    error('Load_Spinup_Meteo:NotHourly', ...
        '%s must be a continuous hourly series to serve as spin-up forcing.', meteo_file);
end
day = floor(Date + one_second);        % stamps sit on the hour; guard rounding

if isnan(start_dn) || isnan(end_dn)
    [yr,~,~] = datevec(Date);
    years = unique(yr)';
    complete = false(size(years));
    for k = 1:numel(years)
        y0 = datenum(years(k),1,1);
        y1 = datenum(years(k)+1,1,1) - 1/24;
        complete(k) = any(abs(Date - y0) <= one_second) && ...
            any(abs(Date - y1) <= one_second);
    end
    if ~any(complete)
        error('Load_Spinup_Meteo:NoCompleteYear', ...
            'No complete calendar year (Jan 1 00:00 to Dec 31 23:00) in %s.', meteo_file);
    end
    if isnan(start_dn); start_dn = datenum(years(find(complete,1,'first')),1,1);   end
    if isnan(end_dn);   end_dn   = datenum(years(find(complete,1,'last')),12,31);  end
end
target_days = (floor(start_dn):floor(end_dn))';
keep = day >= target_days(1) & day <= target_days(end);
%% Whole days only: every day of the window must carry all 24 hourly
%% stamps.  The first and last days of a forcing file are usually partial
%% (its UTC boundaries fall mid-day in local time: Data_Kantola.mat runs
%% 2008-12-31 18:00 to 2020-12-31 17:00), and a partial day accepted as a
%% whole one would leave a cycle or handoff window silently hours short.
counts = accumarray(day(keep) - target_days(1) + 1, 1, [numel(target_days) 1]);
if any(counts == 0)
    error('Load_Spinup_Meteo:MissingDays', ...
        '%s is missing %d day(s) of the window %s to %s.', meteo_file, nnz(counts == 0), ...
        datestr(target_days(1),'yyyy-mm-dd'), datestr(target_days(end),'yyyy-mm-dd'));
end
if any(counts ~= 24)
    partial = find(counts ~= 24, 1, 'first');
    error('Load_Spinup_Meteo:PartialDays', ...
        ['%s has %d incomplete day(s) in the window %s to %s (first: %s, %d of ' ...
         '24 hours); end the window on the file''s last complete day.'], ...
        meteo_file, nnz(counts ~= 24), datestr(target_days(1),'yyyy-mm-dd'), ...
        datestr(target_days(end),'yyyy-mm-dd'), ...
        datestr(target_days(partial),'yyyy-mm-dd'), counts(partial));
end

met = struct();
for k = 1:numel(fields)
    v = S.(fields{k})(:);
    if numel(v) ~= numel(Date)
        error('Load_Spinup_Meteo:LengthMismatch', ...
            '%s: %s has %d values for %d hours.', meteo_file, fields{k}, numel(v), numel(Date));
    end
    met.(fields{k}) = v(keep);
end

site = {'Lat',Lat,0.01; 'Lon',Lon,0.01; 'DeltaGMT',DeltaGMT,0.01};
for k = 1:size(site,1)
    nm = site{k,1};
    if isfield(S,nm) && isnumeric(S.(nm)) && isscalar(S.(nm)) && ...
            isnumeric(site{k,2}) && isscalar(site{k,2}) && ...
            abs(S.(nm) - site{k,2}) > site{k,3}
        warning('Load_Spinup_Meteo:SiteMismatch', ...
            ['%s in %s is %g but the case row says %g: the forcing was ' ...
             'prepared for another location or clock.'], nm, meteo_file, S.(nm), site{k,2});
    end
end

info = struct('source', meteo_file, 'kind', 'hourly forcing file, used as is', ...
    'n_hours', nnz(keep), 'n_days', numel(target_days), ...
    'start_date', target_days(1), 'end_date', target_days(end));
end
