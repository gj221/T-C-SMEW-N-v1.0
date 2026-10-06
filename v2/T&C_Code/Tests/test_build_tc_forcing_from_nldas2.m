function test_build_tc_forcing_from_nldas2
% NLDAS-2 hourly primary forcing must be converted to the T&C forcing
% contract, from hourly NetCDF files on a grid window or from one text time
% series per variable, with NLDAS-2's time convention: state variables and
% shortwave instantaneous at the hour, precipitation the total of the hour
% ending there.
fprintf('\n=== test_build_tc_forcing_from_nldas2 ===\n');
test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
preprocess_dir = fullfile(code_dir,'Inputs','Preprocess');
addpath(code_dir,preprocess_dir);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
    function ok = raises(fn,id)
        ok = false;
        try; fn(); catch ME
            ok = strcmp(ME.identifier,id);
            if ~ok; fprintf('    raised %s: %s\n',ME.identifier,ME.message); end
        end
    end

work_dir = tempname; mkdir(work_dir);
file_cleanup = onCleanup(@() rmdir(work_dir,'s'));
Lat = 41.437; Lon = -72.866; Zbas = 97; DeltaGMT = -5;
lat_c = [41.3125 41.4375]; lon_c = [-72.9375 -72.8125];   % the site is in cell (2,2)

%% --- a day and a bit of hourly NLDAS-2, UTC, at the site cell ---------------
n = 30;
t0 = datenum(2020,6,1,0,0,0);
t = t0 + (0:n-1).'/24;
hod = mod((0:n-1).',24);
local_hour = mod(hod + DeltaGMT,24);
Tair = 288 + 6*sin(2*pi*(local_hour-9)/24);
Qair = 0.008 + 0.001*sin(2*pi*hod/24);
Qair(5) = 0.05;                                   % supersaturated: capped at esat
PSurf = 100500 + 50*cos(2*pi*hod/24);
Wind_E = 3*ones(n,1); Wind_N = 4*ones(n,1);
SWdown = max(800*sin(pi*(local_hour-6)/12),0);
Rainf = zeros(n,1); Rainf([8 9 20]) = [0.4 1.2 0.1];
site = [Tair Qair PSurf Wind_E Wind_N SWdown Rainf];
names = {'Tair','Qair','PSurf','Wind_E','Wind_N','SWdown','Rainf'};
units = {'K','kg kg-1','Pa','m s-1','m s-1','W m-2','kg m-2'};

%% hourly NetCDF files, 2 x 2 cells: the other cells hold other values
nc_dir = fullfile(work_dir,'nc'); mkdir(nc_dir);
for h = 1:n
    f = fullfile(nc_dir,sprintf('NLDAS_FORA0125_H.A%s.020.nc',datestr(t(h),'yyyymmdd.HHMM')));
    nccreate(f,'lon','Dimensions',{'lon',2}); ncwrite(f,'lon',lon_c);
    nccreate(f,'lat','Dimensions',{'lat',2}); ncwrite(f,'lat',lat_c);
    nccreate(f,'time','Dimensions',{'time',1});
    ncwrite(f,'time',(t(h)-datenum(1979,1,1))*24);
    ncwriteatt(f,'time','units','hours since 1979-01-01 00:00:00 UTC');
    for k = 1:numel(names)
        grid = repmat(site(h,k)*1.5 + 7,[2 2]);   % neighbours: not the site
        grid(2,2) = site(h,k);
        nccreate(f,names{k},'Dimensions',{'lon',2,'lat',2,'time',1},'FillValue',9.999e20);
        ncwrite(f,names{k},grid);
        ncwriteatt(f,names{k},'units',units{k});
    end
end
out_nc = fullfile(work_dir,'Data_nc.mat');
Build_TC_Forcing_From_NLDAS2(nc_dir,out_nc,Lat,Lon,Zbas,DeltaGMT);
F = load(out_nc);

required = {'Date','Pr','Ta','Ws','ea','Pre','SAD1','SAD2','SAB1','SAB2', ...
    'N','Tdew','esat','PARB','PARD','Lat','Lon','Zbas','DeltaGMT'};
check('the forcing contract is complete',all(isfield(F,required)));
check('one record per hour after the first, hour-ending local standard time', ...
    numel(F.Date) == n-1 && abs(F.Date(1) - (t(2) + DeltaGMT/24)) < 1e-9 && ...
    all(abs(diff(F.Date) - 1/24) < 1e-9));
check('the site cell is read, not its neighbours', ...
    max(abs(F.Ta - (Tair(2:end)-273.15))) < 1e-9 && max(abs(F.Pre - PSurf(2:end)/100)) < 1e-9);
check('precipitation is the total of the hour ending at the stamp', ...
    isequal(F.Pr,Rainf(2:end)));
Rsw = SAD1_SAD2(F);
check('shortwave is the mean of its instantaneous values at the hour''s ends', ...
    max(abs(Rsw - (SWdown(1:end-1)+SWdown(2:end))/2)) < 1e-8);
check('wind speed combines its components',max(abs(F.Ws - 5)) < 1e-12);
e = Qair(2:end).*PSurf(2:end)./(0.622 + 0.378*Qair(2:end));
es = 610.8*exp(17.27*F.Ta./(237.3+F.Ta));
check('vapour pressure from specific humidity, capped at saturation', ...
    max(abs(F.ea - min(e,es))) < 1e-9 && abs(F.ea(4) - es(4)) < 1e-9 && F.ea(4) < e(4));
check('the dew point is the vapour pressure''s own, never above the air', ...
    max(abs(610.8*exp(17.27*F.Tdew./(237.3+F.Tdew)) - F.ea)) < 1e-6 && all(F.Tdew <= F.Ta + 1e-9));
check('the site is recorded',F.Lat == Lat && F.Lon == Lon && F.Zbas == Zbas && F.DeltaGMT == DeltaGMT);

%% --- the same series as one text file per variable ---------------------------
txt_dir = fullfile(work_dir,'txt'); mkdir(txt_dir);
for k = 1:numel(names)
    write_series(fullfile(txt_dir,sprintf('NLDAS_FORA0125_H_v2.0_%s.txt',names{k})), ...
        names{k},units{k},t,site(:,k),{});
end
out_txt = fullfile(work_dir,'Data_txt.mat');
Build_TC_Forcing_From_NLDAS2(txt_dir,out_txt,Lat,Lon,Zbas,DeltaGMT);
G = load(out_txt);
same = true;
for k = 1:numel(required)
    same = same && isequal(size(F.(required{k})),size(G.(required{k}))) && ...
        max(abs(F.(required{k})(:) - G.(required{k})(:))) < 1e-9;
end
check('text series give the forcing the NetCDF files give',same);

%% --- Giovanni time-series CSV (NASA Harmony point request), each variable in
%% two files that are joined --------------------------------------------------
gio_dir = fullfile(work_dir,'giovanni'); mkdir(gio_dir);
first = 1:15; rest = 16:n;
for k = 1:numel(names)
    write_giovanni(fullfile(gio_dir,sprintf('NLDAS_FORA0125_H_2.0_%s.csv',names{k})), ...
        names{k},units{k},t(first),site(first,k),[41.4375 -72.8125]);
    write_giovanni(fullfile(gio_dir,sprintf('NLDAS_FORA0125_H_2.0_%s_2026.csv',names{k})), ...
        names{k},units{k},t(rest),site(rest,k),[41.4375 -72.8125]);
end
out_gio = fullfile(work_dir,'Data_gio.mat');
Build_TC_Forcing_From_NLDAS2(gio_dir,out_gio,Lat,Lon,Zbas,DeltaGMT);
H = load(out_gio);
same = true;
for k = 1:numel(required)
    same = same && isequal(size(F.(required{k})),size(H.(required{k}))) && ...
        max(abs(F.(required{k})(:) - H.(required{k})(:))) < 1e-9;
end
check('Giovanni CSV series, joined across files, give the NetCDF files'' forcing',same);
gio_far = fullfile(work_dir,'giovanni_far'); copyfile(gio_dir,gio_far);
write_giovanni(fullfile(gio_far,'NLDAS_FORA0125_H_2.0_Qair.csv'),'Qair','kg kg-1', ...
    t(first),Qair(first),[40.0625 -88.1875]);
check('a Giovanni series of another grid cell is refused', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    gio_far,fullfile(work_dir,'x.mat'),Lat,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:SiteOutsideGrid'));
gio_units = fullfile(work_dir,'giovanni_units'); copyfile(gio_dir,gio_units);
write_giovanni(fullfile(gio_units,'NLDAS_FORA0125_H_2.0_PSurf.csv'),'PSurf','hPa', ...
    t(first),PSurf(first)/100,[41.4375 -72.8125]);
check('a Giovanni series in other units is refused', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    gio_units,fullfile(work_dir,'x.mat'),Lat,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:UnexpectedUnits'));

%% --- refusals ---------------------------------------------------------------
gap_dir = fullfile(work_dir,'gap'); copyfile(txt_dir,gap_dir);
write_series(fullfile(gap_dir,'NLDAS_FORA0125_H_v2.0_Rainf.txt'),'Rainf','kg m-2', ...
    t([1:11 13:end]),Rainf([1:11 13:end]),{});
check('a missing hour is refused, not filled', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    gap_dir,fullfile(work_dir,'x.mat'),Lat,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:MissingHours'));
check('a site outside the grid window is refused', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    nc_dir,fullfile(work_dir,'x.mat'),41.7,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:SiteOutsideGrid'));
unit_dir = fullfile(work_dir,'units'); copyfile(txt_dir,unit_dir);
write_series(fullfile(unit_dir,'NLDAS_FORA0125_H_v2.0_Tair.txt'),'Tair','C',t,Tair-273.15,{});
check('a variable in other units is refused', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    unit_dir,fullfile(work_dir,'x.mat'),Lat,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:UnexpectedUnits'));
far_dir = fullfile(work_dir,'far'); copyfile(txt_dir,far_dir);
write_series(fullfile(far_dir,'NLDAS_FORA0125_H_v2.0_Qair.txt'),'Qair','kg kg-1',t,Qair, ...
    {'lat=40.0625','lon=-88.1875'});
check('a series for another place is refused', raises(@() Build_TC_Forcing_From_NLDAS2( ...
    far_dir,fullfile(work_dir,'x.mat'),Lat,Lon,Zbas,DeltaGMT),'Build_TC_Forcing_From_NLDAS2:SiteOutsideGrid'));

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_build_tc_forcing_from_nldas2: %d check(s) failed',n_fail); end
end

function Rsw = SAD1_SAD2(F)
Rsw = F.SAD1 + F.SAD2 + F.SAB1 + F.SAB2;
end

function write_giovanni(file,name,unit,t,x,latlon)
%% The Giovanni time-series CSV as NASA Harmony serves it: key,value header,
%% a blank line, "Timestamp (UTC),Data", then yyyy-mm-dd HH:MM,value.
fid = fopen(file,'w');
fprintf(fid,'prod_name,NLDAS_FORA0125_H.2.0\n');
fprintf(fid,'param_short_name,%s\n',name);
fprintf(fid,'param_name,%s\n',name);
fprintf(fid,'unit,%s\n',unit);
fprintf(fid,'undef,-9999\n');
fprintf(fid,'begin_time,%s\n',datestr(t(1),'yyyy-mm-dd HH:MM:SS'));
fprintf(fid,'end_time,%s\n',datestr(t(end),'yyyy-mm-dd HH:MM:SS'));
fprintf(fid,'lat,%.4f\nlon,%.4f\nlat_resolution,0.125\nlon_resolution,0.125\n',latlon(1),latlon(2));
fprintf(fid,'\nTimestamp (UTC),Data\n');
for h = 1:numel(t)
    fprintf(fid,'%s,%.17g\n',datestr(t(h),'yyyy-mm-dd HH:MM'),x(h));
end
fclose(fid);
end

function write_series(file,name,unit,t,x,extra)
%% A Data Rods-like ASCII series: key=value header, then time<TAB>value.
fid = fopen(file,'w');
fprintf(fid,'prod_name=NLDAS_FORA0125_H_v2.0\n');
fprintf(fid,'param_short_name=%s\n',name);
fprintf(fid,'unit=%s\n',unit);
for k = 1:numel(extra); fprintf(fid,'%s\n',extra{k}); end
fprintf(fid,'Date&Time               Data\n');
for h = 1:numel(t)
    fprintf(fid,'%s\t%.17g\n',datestr(t(h),'yyyy-mm-ddTHH:MM:SS'),x(h));
end
fclose(fid);
end
