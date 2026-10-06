function Build_TC_Forcing_From_NLDAS2(nldas_source,output_file,Lat,Lon,Zbas,DeltaGMT)
% Convert NLDAS-2 primary forcing (NLDAS_FORA0125_H, v2.0) for one site to
% hourly T&C forcing, by the route of Build_TC_Forcing_From_ERA5_Land_Timeseries.
%
% nldas_source is a folder holding either
%   - NetCDF files (.nc, .nc4): the hourly NLDAS_FORA0125_H files from GES
%     DISC (direct, OPeNDAP or the subsetter, any window around the site)
%     or per-variable point time series; the grid cell containing Lat/Lon
%     is read, its time axis from the CF units ("<unit> since <origin>"); or
%   - text files (.txt, .asc, .csv), one variable each: point time series,
%     such as the Giovanni time-series service's CSV (through NASA Harmony:
%     a POINT request on the collection, one file per variable) or Data
%     Rods ASCII.  A file may hold any stretch of the record; a variable's
%     files are joined.  A line is data when it begins with an ISO time
%     (yyyy-mm-ddTHH[:MM[:SS]], Z optional, or "yyyy-mm-dd HHZ") followed
%     by a value; the header's keys (key=value, key: value or key,value)
%     name the variable (param_short_name, variable) -- else the file name
%     does -- its units (unit, units) and the grid cell (lat, lon) when
%     they are there.
% Variables (v2.0 names; the GRIB short names of v002 are accepted):
% Tair [K], Qair [kg/kg], PSurf [Pa], Wind_E and Wind_N [m/s], SWdown
% [W/m2], Rainf [kg/m2].  Times are UTC.  Lat and Lon are degrees, Zbas is
% m a.s.l., and DeltaGMT the fixed hour offset of local standard time.
%
% NLDAS-2 time convention (NLDAS-2 README, Table 8): air temperature,
% humidity, pressure, wind and shortwave are instantaneous at the hour
% stamped; precipitation is the total of the hour ending there.  T&C's
% forcing is hour-ending, as in the ERA5 route: the state variables are
% taken at the end of the hour (as ERA5's instantaneous t2m, d2m, sp and
% wind are), precipitation as the hour's total, and shortwave as the
% hour's mean -- the average of the instantaneous values at its two ends,
% where ERA5 gives the accumulation.  The first hour has no start value,
% so the forcing begins one hour after the first record.
%
% Humidity: vapour pressure from specific humidity and pressure,
% e = q p/(0.622 + 0.378 q), capped at saturation (the Magnus form over
% water used by the ERA5 route); the dew point is the Magnus inverse of
% that capped e, so ea, esat, Tdew and Ta stay mutually consistent.

narginchk(6,6);
nldas_source = char(nldas_source);
output_file = char(output_file);
if ~isfolder(nldas_source)
    error('Build_TC_Forcing_From_NLDAS2:MissingInput', ...
        'NLDAS-2 source folder not found: %s',nldas_source);
end
[output_dir,~,output_ext] = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir)
    error('Build_TC_Forcing_From_NLDAS2:MissingOutputDirectory', ...
        'Output directory not found: %s',output_dir);
end
if ~strcmpi(output_ext,'.mat')
    error('Build_TC_Forcing_From_NLDAS2:BadOutputExtension', ...
        'Output file must have a .mat extension.');
end
validateattributes(Lat,{'numeric'},{'scalar','real','finite','>=',-90,'<=',90});
validateattributes(Lon,{'numeric'},{'scalar','real','finite','>=',-180,'<=',180});
validateattributes(Zbas,{'numeric'},{'scalar','real','finite'});
validateattributes(DeltaGMT,{'numeric'},{'scalar','real','finite','>=',-12,'<=',14});

V = variables();
series = read_source(nldas_source,V,Lat,Lon);
[t_utc,X] = common_hours(series,V);

%% hour-ending forcing: drop the first hour (no start value for shortwave)
Tair = X(2:end,1); Qair = X(2:end,2); PSurf = X(2:end,3);
Wind_E = X(2:end,4); Wind_N = X(2:end,5);
SW_mean = (X(1:end-1,6) + X(2:end,6))/2;
Rainf = X(2:end,7);
t_utc = t_utc(2:end);

if any(Tair <= 150 | Tair >= 350)
    error('Build_TC_Forcing_From_NLDAS2:BadTemperature', ...
        'Tair must be in kelvin (a value lies outside 150-350 K).');
end
if any(PSurf <= 3e4 | PSurf >= 1.1e5)
    error('Build_TC_Forcing_From_NLDAS2:BadPressure', ...
        'PSurf must be in pascals (a value lies outside 30-110 kPa).');
end
if any(Qair < 0 | Qair >= 0.1)
    error('Build_TC_Forcing_From_NLDAS2:BadHumidity', ...
        'Qair must be specific humidity in kg/kg.');
end
if any(Rainf < -1e-6)
    error('Build_TC_Forcing_From_NLDAS2:NegativePrecipitation', ...
        'Rainf contains a negative value.');
end
if any(SW_mean < -1e-6 | SW_mean > 1500)
    error('Build_TC_Forcing_From_NLDAS2:BadRadiation', ...
        'SWdown must be in W/m2 (a value lies outside 0-1500).');
end

Ta = Tair - 273.15;
Pre = PSurf/100;
Ws = sqrt(Wind_E.^2 + Wind_N.^2);
Ws(Ws==0) = 0.01;
Pr = max(Rainf,0);
Rsw = max(SW_mean,0);
esat = 610.8*exp(17.27*Ta./(237.3+Ta));
ea = min(Qair.*PSurf./(0.622 + 0.378*Qair),esat);
if any(ea <= 0)
    error('Build_TC_Forcing_From_NLDAS2:BadHumidity', ...
        'Vapour pressure must be positive (Qair is zero somewhere).');
end
g = log(ea/610.8);
Tdew = 237.3*g./(17.27 - g);

local_time = t_utc + hours(DeltaGMT);
local_time.TimeZone = '';
Date = datenum(local_time);
Date = Date(:);
n = numel(Date);

[~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N] = ...
    Automatic_Radiation_Partition_I( ...
    Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,0);
SAD1 = SAD1(:);
SAD2 = SAD2(:);
SAB1 = SAB1(:);
SAB2 = SAB2(:);
PARB = PARB(:);
PARD = PARD(:);
N = N(:);

partition_total = SAD1 + SAD2 + SAB1 + SAB2;
residual = Rsw - partition_total;
nonclosure = abs(residual)>1e-8;
SAD1(nonclosure) = SAD1(nonclosure) + residual(nonclosure)/2;
SAD2(nonclosure) = SAD2(nonclosure) + residual(nonclosure)/2;
if any(~isfinite([SAD1;SAD2;SAB1;SAB2;PARB;PARD;N]))
    error('Build_TC_Forcing_From_NLDAS2:PartitionFailure', ...
        'Radiation partitioning produced NaN or infinite values.');
end
if any([SAD1;SAD2;SAB1;SAB2;PARB;PARD] < -1e-8)
    error('Build_TC_Forcing_From_NLDAS2:NegativePartition', ...
        'Radiation partitioning produced a negative flux.');
end
if max(abs(SAD1+SAD2+SAB1+SAB2-Rsw))>1e-8
    error('Build_TC_Forcing_From_NLDAS2:RadiationNonclosure', ...
        'Partitioned shortwave radiation does not close to NLDAS-2 radiation.');
end

save(output_file,'Date','Pr','Ta','Ws','ea','Pre', ...
    'SAD1','SAD2','SAB1','SAB2','N','Tdew','esat','PARB','PARD', ...
    'Lat','Lon','Zbas','DeltaGMT','-v7');
fprintf('Wrote %s with %d hourly forcing records (%s to %s local standard time).\n', ...
    output_file,n,datestr(Date(1),'yyyy-mm-dd HH:MM'),datestr(Date(end),'yyyy-mm-dd HH:MM'));
end

%% =======================================================================
function V = variables()
%% The seven forcing variables in their working order, with the names they
%% go by and the units they must carry.
V = struct( ...
    'name',  {'Tair','Qair','PSurf','Wind_E','Wind_N','SWdown','Rainf'}, ...
    'alias', {{'tair','tmp2m','tmp'},{'qair','spfh2m','spfh'},{'psurf','pressfc','pres'}, ...
              {'wind_e','ugrd10m','ugrd'},{'wind_n','vgrd10m','vgrd'}, ...
              {'swdown','dswrfsfc','dswrf'},{'rainf','apcpsfc','apcp'}}, ...
    'units', {{'k','kelvin'},{'kgkg-1','kg/kg','1','kgkg**-1'},{'pa'}, ...
              {'ms-1','m/s','ms**-1'},{'ms-1','m/s','ms**-1'}, ...
              {'wm-2','w/m2','wm**-2','w/m^2'},{'kgm-2','kg/m2','mm','kgm**-2'}});
end

function series = read_source(folder,V,Lat,Lon)
nc = [dir(fullfile(folder,'*.nc')); dir(fullfile(folder,'*.nc4'))];
txt = [dir(fullfile(folder,'*.txt')); dir(fullfile(folder,'*.asc')); dir(fullfile(folder,'*.csv'))];
if isempty(nc) && isempty(txt)
    error('Build_TC_Forcing_From_NLDAS2:NoFiles', ...
        'No NetCDF (.nc, .nc4) or text (.txt, .asc, .csv) files in %s.',folder);
end
if ~isempty(nc) && ~isempty(txt)
    error('Build_TC_Forcing_From_NLDAS2:MixedFormats', ...
        '%s holds both NetCDF and text files; keep one format per folder.',folder);
end
parts = repmat({{}},1,numel(V));
if ~isempty(nc)
    [~,order] = sort({nc.name}); nc = nc(order);
    for f = 1:numel(nc)
        file = fullfile(nc(f).folder,nc(f).name);
        info = ncinfo(file);
        names = {info.Variables.Name};
        t = []; cell_ij = [];
        for k = 1:numel(V)
            j = find(ismember(lower(names),V(k).alias),1);
            if isempty(j); continue; end
            if isempty(t)
                t = nc_time(file,info);
                cell_ij = nc_cell(file,info,Lat,Lon);
            end
            require_units(file,names{j},V(k).units,ncreadatt_or_empty(file,names{j},'units'));
            parts{k}{end+1} = [t nc_point(file,info,names{j},cell_ij,numel(t))];
        end
    end
else
    [~,order] = sort({txt.name}); txt = txt(order);
    for f = 1:numel(txt)
        file = fullfile(txt(f).folder,txt(f).name);
        [k,t,x] = read_text_series(file,V,Lat,Lon);
        parts{k}{end+1} = [t x];
    end
end
series = cell(1,numel(V));
for k = 1:numel(V)
    if isempty(parts{k})
        error('Build_TC_Forcing_From_NLDAS2:MissingVariable', ...
            'No file holds %s (or %s).',V(k).name,strjoin(V(k).alias,', '));
    end
    series{k} = vertcat(parts{k}{:});
end
end

function t = nc_time(file,info)
names = {info.Variables.Name};
j = find(strcmpi(names,'time'),1);
if isempty(j)
    error('Build_TC_Forcing_From_NLDAS2:NoTime','%s has no time variable.',file);
end
units = ncreadatt_or_empty(file,names{j},'units');
parts = regexp(strtrim(units),'^(seconds|minutes|hours|days)\s+since\s+(.+)$', ...
    'tokens','once','ignorecase');
if isempty(parts)
    error('Build_TC_Forcing_From_NLDAS2:BadTimeUnits', ...
        '%s: unsupported time units "%s".',file,units);
end
origin_text = regexprep(strtrim(parts{2}),'(\.\d+)?\s*(UTC|Z|[+]00:?00)?$','');
origin = parse_time(regexprep(origin_text,'(\d)T(\d)','$1 $2'),file);
v = double(ncread(file,names{j}));
switch lower(parts{1})
    case 'seconds'; t = origin + seconds(v(:));
    case 'minutes'; t = origin + minutes(v(:));
    case 'hours';   t = origin + hours(v(:));
    case 'days';    t = origin + days(v(:));
end
t = datenum(t);
end

function ij = nc_cell(file,info,Lat,Lon)
%% The grid cell whose centre is nearest the site; it must contain it
%% (NLDAS-2 cells are 0.125 degree).  A file without lat/lon is a point
%% series and has no cell to choose.
names = {info.Variables.Name};
jl = find(ismember(lower(names),{'lat','latitude'}),1);
jo = find(ismember(lower(names),{'lon','longitude'}),1);
ij = [];
if isempty(jl) || isempty(jo); return; end
la = double(ncread(file,names{jl})); lo = double(ncread(file,names{jo}));
lo = mod(lo+180,360)-180;
[dla,i] = min(abs(la(:)-Lat)); [dlo,j] = min(abs(lo(:)-Lon));
if dla > 0.0625+1e-9 || dlo > 0.0625+1e-9
    error('Build_TC_Forcing_From_NLDAS2:SiteOutsideGrid', ...
        '%s does not cover %.4f N %.4f E (nearest cell centre %.4f N %.4f E).', ...
        file,Lat,Lon,la(i),lo(j));
end
if isscalar(la) && isscalar(lo); return; end   % a point series: nothing to choose
ij = struct('lat_name',names{jl},'lon_name',names{jo},'i',i,'j',j);
end

function x = nc_point(file,info,name,ij,nt)
%% The variable at the site's cell, one value per time step.
v = info.Variables(strcmp({info.Variables.Name},name));
dims = {v.Dimensions.Name};
start = ones(1,numel(dims)); count = [v.Dimensions.Length];
if ~isempty(ij)
    a = find(strcmp(dims,ij.lat_name) | strcmpi(dims,'lat') | strcmpi(dims,'latitude'),1);
    b = find(strcmp(dims,ij.lon_name) | strcmpi(dims,'lon') | strcmpi(dims,'longitude'),1);
    if isempty(a) || isempty(b)
        error('Build_TC_Forcing_From_NLDAS2:BadDimensions', ...
            '%s: %s is not laid out on lat and lon.',file,name);
    end
    start(a) = ij.i; count(a) = 1; start(b) = ij.j; count(b) = 1;
end
x = double(ncread(file,name,start,count));
x = x(:);
if numel(x) ~= nt
    error('Build_TC_Forcing_From_NLDAS2:BadDimensions', ...
        '%s: %s does not hold one value per time step at the site.',file,name);
end
x(abs(x) >= 1e19) = NaN;   % the NLDAS fill value, where not already NaN
end

function [k,t,x] = read_text_series(file,V,Lat,Lon)
lines = splitlines(string(fileread(file)));
lines = strtrim(lines(strlength(strtrim(lines)) > 0));
%% (MATLAB captures nothing inside a non-capturing group, so minutes and
%% seconds are plain optional groups, colon included.)
pat = '^(\d{4})-(\d{2})-(\d{2})[T ](\d{2})(:\d{2})?(:\d{2})?Z?[\s,;]+([-+]?[\d.]+(?:[eE][-+]?\d+)?|NaN)\s*$';
tok = regexp(lines,pat,'tokens','once');
is_data = ~cellfun(@isempty,tok);
if ~any(is_data)
    error('Build_TC_Forcing_From_NLDAS2:NoData','%s holds no "time value" lines.',file);
end
header = lines(~is_data);
key = @(names) header_value(header,names);
%% the variable: a header key, else the file name
name = key({'param_short_name','variable','short_name','parameter'});
[~,base] = fileparts(file);
k = match_variable(V,name);
if isempty(k); k = match_variable(V,base); end
if isempty(k)
    error('Build_TC_Forcing_From_NLDAS2:UnknownVariable', ...
        '%s: no header key or file name names one of %s.',file,strjoin({V.name},', '));
end
units = key({'unit','units'});
if ~isempty(units); require_units(file,V(k).name,V(k).units,units); end
la = str2double(key({'lat','latitude'})); lo = str2double(key({'lon','longitude'}));
if isfinite(la) && isfinite(lo) && (abs(la-Lat) > 0.0625+1e-9 || ...
        abs(mod(lo+180,360)-180-Lon) > 0.0625+1e-9)
    error('Build_TC_Forcing_From_NLDAS2:SiteOutsideGrid', ...
        '%s is for %.4f N %.4f E, not the site %.4f N %.4f E.',file,la,lo,Lat,Lon);
end
T = vertcat(tok{is_data});
num = str2double(erase(T(:,1:6),':')); num(isnan(num)) = 0;
t = datenum(num(:,1),num(:,2),num(:,3),num(:,4),num(:,5),num(:,6));
x = str2double(T(:,7));
x(x <= -9998) = NaN;   % the usual missing-value codes
end

function v = header_value(header,names)
v = '';
for n = 1:numel(names)
    hit = regexp(header,['^' names{n} '\s*[=:,]\s*(.+)$'],'tokens','once','ignorecase');
    hit = hit(~cellfun(@isempty,hit));
    if ~isempty(hit)
        v = regexprep(char(hit{1}{1}),'^[\s,;"]+|[\s,;"]+$','');   % "Unit:,K," -> K
        return
    end
end
end

function k = match_variable(V,txt)
k = [];
if isempty(txt); return; end
low = lower(char(txt));
for c = 1:numel(V)
    for a = 1:numel(V(c).alias)
        if ~isempty(regexp(low,['(^|[^a-z0-9])' regexptranslate('escape',V(c).alias{a}) '($|[^a-z0-9])'],'once'))
            k = c; return
        end
    end
end
end

function [t_utc,X] = common_hours(series,V)
%% Every variable on one gap-free hourly UTC axis: the hours all of them
%% cover.  A duplicate hour must agree; a missing hour inside is an error.
lo = -Inf; hi = Inf;
for k = 1:numel(series)
    s = sortrows(series{k},1);
    [tu,first] = unique(round(s(:,1)*24));
    dup = setdiff(1:size(s,1),first);
    if ~isempty(dup)
        again = round(s(dup,1)*24);
        [~,where] = ismember(again,tu);
        if any(abs(s(dup,2) - s(first(where),2)) > 1e-6*max(1,abs(s(dup,2))))
            error('Build_TC_Forcing_From_NLDAS2:DuplicateHours', ...
                '%s has two different values for one hour.',V(k).name);
        end
    end
    series{k} = [tu s(first,2)];
    lo = max(lo,tu(1)); hi = min(hi,tu(end));
end
if hi - lo < 1
    error('Build_TC_Forcing_From_NLDAS2:NoCommonHours', ...
        'The variables share fewer than two hours.');
end
hours_all = (lo:hi).';
X = nan(numel(hours_all),numel(series));
for k = 1:numel(series)
    [ok,where] = ismember(hours_all,series{k}(:,1));
    X(ok,k) = series{k}(where(ok),2);
    gap = find(~ok | isnan(X(:,k)),1);
    if ~isempty(gap)
        error('Build_TC_Forcing_From_NLDAS2:MissingHours', ...
            '%s has no value for %s UTC.',V(k).name,datestr(hours_all(gap)/24,'yyyy-mm-dd HH:MM'));
    end
end
t_utc = datetime(hours_all/24,'ConvertFrom','datenum','TimeZone','UTC');
end

function require_units(file,name,allowed,units)
if isempty(units)
    error('Build_TC_Forcing_From_NLDAS2:MissingUnits', ...
        '%s: variable %s has no units attribute.',file,name);
end
u = lower(regexprep(char(units),'[\s\^\*\{\}]',''));
allowed = regexprep(lower(allowed),'[\s\^\*\{\}]','');
if ~any(strcmp(u,allowed))
    error('Build_TC_Forcing_From_NLDAS2:UnexpectedUnits', ...
        '%s: unexpected units for %s: %s',file,name,char(units));
end
end

function v = ncreadatt_or_empty(file,name,att)
try
    v = char(ncreadatt(file,name,att));
catch
    v = '';
end
end

function t = parse_time(txt,file)
formats = {'yyyy-MM-dd HH:mm:ss','yyyy-MM-dd HH:mm','yyyy-MM-dd'};
for k = 1:numel(formats)
    try
        t = datetime(txt,'InputFormat',formats{k},'TimeZone','UTC');
        return
    catch
    end
end
error('Build_TC_Forcing_From_NLDAS2:BadTimeOrigin', ...
    '%s: could not parse the time origin "%s".',file,txt);
end
