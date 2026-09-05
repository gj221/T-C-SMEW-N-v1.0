function test_spinup_meteo_loader()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_spinup_meteo_loader
%%
%% Load_Spinup_Meteo on a synthetic hourly forcing .mat:
%%   1  a requested window is sliced to whole days, hour-exact
%%   2  NaN bounds resolve to the first/last COMPLETE calendar years
%%   3  the sliced years pass Validate_Spinup_Cycle_Years unchanged
%%   4  a window the file does not cover is an error, not a short series
%%   5  a non-hourly file is refused
%%   6  a missing forcing variable is named in the error
%%   7  a site mismatch (Lat) is warned about, the data still returned
%%   8  a non-.mat name is handed to Load_Daily_Meteo (missing file error
%%      comes from there)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
fprintf('\n=== test_spinup_meteo_loader ===\n');
n_pass = 0; n_fail = 0;
    function check(ok,msg)
        if ok; n_pass = n_pass+1; fprintf('  PASS  %s\n',msg);
        else;  n_fail = n_fail+1; fprintf('  FAIL  %s\n',msg); end
    end

%% synthetic hourly forcing: 2013-07-01 00:00 .. 2016-03-31 23:00 local,
%% so 2014 and 2015 are the complete years
tmp = tempname; mkdir(tmp);
d0 = datenum(2013,7,1); d1 = datenum(2016,4,1) - 1/24;
Date = (d0:1/24:d1)';
nh = numel(Date);
S = struct('Date',Date,'Lat',40.06,'Lon',-88.20,'DeltaGMT',-6);
vars = {'Pr','Ta','Ws','ea','Pre','SAD1','SAD2','SAB1','SAB2','N','Tdew','esat','PARB','PARD'};
for k = 1:numel(vars); S.(vars{k}) = (1:nh)' + k; end   % distinct, index-coded
f = fullfile(tmp,'forcing.mat'); save(f,'-struct','S');

%% 1 window slicing
[met,info] = Load_Spinup_Meteo(f,datenum(2014,3,1),datenum(2014,3,10),40.06,-88.20,-6);
check(numel(met.Date) == 10*24 && abs(met.Date(1)-datenum(2014,3,1)) < 1e-9 && ...
    abs(met.Date(end)-(datenum(2014,3,11)-1/24)) < 1e-9, ...
    'a 10-day window returns 240 hours from 00:00 of the first day to 23:00 of the last');
i0 = find(abs(Date-datenum(2014,3,1)) < 1e-9);
check(isequal(met.Pr, S.Pr(i0:i0+239)) && isequal(met.PARD, S.PARD(i0:i0+239)), ...
    'every variable is sliced to the same rows');
check(info.n_hours == 240 && info.n_days == 10 && strcmp(info.kind,'hourly forcing file, used as is'), ...
    'info reports hours, days and the hourly kind');

%% 2 NaN bounds -> complete years
[met2,~] = Load_Spinup_Meteo(f,NaN,NaN,40.06,-88.20,-6);
[yy,~,~] = datevec(met2.Date);
check(min(yy) == 2014 && max(yy) == 2015 && numel(met2.Date) == 2*365*24, ...
    'NaN bounds resolve to the complete years 2014-2015 only');

%% 3 cycle-year validator accepts the slice
try
    Validate_Spinup_Cycle_Years(met2.Date,[2014 2015]); ok = true;
catch; ok = false; end
check(ok,'Validate_Spinup_Cycle_Years accepts the sliced complete years');

%% 4 uncovered window
ok = false;
try; Load_Spinup_Meteo(f,datenum(2016,3,25),datenum(2016,4,5),40.06,-88.20,-6);
catch ME; ok = strcmp(ME.identifier,'Load_Spinup_Meteo:MissingDays'); end
check(ok,'a window running past the file is an error naming the missing days');

%% 5 non-hourly
S3 = S; S3.Date = Date(1:3:end);
for k = 1:numel(vars); S3.(vars{k}) = S.(vars{k})(1:3:end); end
f3 = fullfile(tmp,'three_hourly.mat'); save(f3,'-struct','S3');
ok = false;
try; Load_Spinup_Meteo(f3,NaN,NaN,40.06,-88.20,-6);
catch ME; ok = strcmp(ME.identifier,'Load_Spinup_Meteo:NotHourly'); end
check(ok,'a 3-hourly file is refused');

%% 6 missing variable
S4 = rmfield(S,'PARD'); f4 = fullfile(tmp,'no_pard.mat'); save(f4,'-struct','S4');
ok = false;
try; Load_Spinup_Meteo(f4,NaN,NaN,40.06,-88.20,-6);
catch ME; ok = strcmp(ME.identifier,'Load_Spinup_Meteo:MissingVariables') && contains(ME.message,'PARD'); end
check(ok,'a missing forcing variable is named');

%% 7 site mismatch warns
lastwarn(''); warning('off','Load_Spinup_Meteo:SiteMismatch');
[met7,~] = Load_Spinup_Meteo(f,NaN,NaN,36.70,-78.50,-5);
[~,wid] = lastwarn; warning('on','Load_Spinup_Meteo:SiteMismatch');
check(strcmp(wid,'Load_Spinup_Meteo:SiteMismatch') && numel(met7.Date) == numel(met2.Date), ...
    'a Lat/Lon/DeltaGMT mismatch is warned about and the data still returned');

%% 8 daily dispatch
ok = false;
try; Load_Spinup_Meteo(fullfile(tmp,'station.txt'),NaN,NaN,40.06,-88.20,-6);
catch ME; ok = strcmp(ME.identifier,'Load_Daily_Meteo:MissingFile'); end
check(ok,'a non-.mat name is dispatched to Load_Daily_Meteo');

rmdir(tmp,'s');
fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_spinup_meteo_loader: %d check(s) failed',n_fail); end
end
