function Ca = Load_CO2_Record(Date, Input_Directory, record_file)
%LOAD_CO2_RECORD  Atmospheric CO2 concentration [ppm] from a CO2 record.
%
%   Ca = LOAD_CO2_RECORD(Date, Input_Directory) returns the atmospheric CO2
%   mixing ratio [ppm] interpolated onto the simulation timestamps Date
%   (MATLAB serial date numbers) from the hourly record in Ca_Data.mat
%   (variables: Ca, Date_CO2), v0.9's historical series (1975-2022).
%
%   Ca = LOAD_CO2_RECORD(Date, Input_Directory, record_file) reads another
%   record of the same form in Input_Directory -- a site's, named by its
%   parameter row's co2_file column (e.g. Ca_Data_NOAA_MLO_daily.mat,
%   NOAA Mauna Loa daily means, interpolated linearly between days).
%
%   A timestamp outside the record is an error, not a NaN: an unset CO2
%   would reach photosynthesis silently.  Each record is cached on first
%   use, so repeated spin-up cycles do not reload it.
if nargin < 3 || isempty(record_file)
    record_file = 'Ca_Data.mat';
end
persistent cached_file Date_CO2 Ca_record
key = fullfile(Input_Directory, char(record_file));
if isempty(cached_file) || ~strcmp(cached_file, key)
    if ~isfile(key)
        error('Load_CO2_Record:MissingRecord', 'CO2 record not found: %s', key);
    end
    S = load(key, 'Date_CO2', 'Ca');
    Date_CO2  = S.Date_CO2(:);
    Ca_record = S.Ca(:);
    cached_file = key;
end
if any(Date(:) < Date_CO2(1) | Date(:) > Date_CO2(end))
    error('Load_CO2_Record:NotCovered', ...
        '%s covers %s to %s; the run asks for %s to %s.', char(record_file), ...
        datestr(Date_CO2(1), 'yyyy-mm-dd HH:MM'), datestr(Date_CO2(end), 'yyyy-mm-dd HH:MM'), ...
        datestr(min(Date(:)), 'yyyy-mm-dd HH:MM'), datestr(max(Date(:)), 'yyyy-mm-dd HH:MM'));
end
Ca = interp1(Date_CO2, Ca_record, Date);
end
