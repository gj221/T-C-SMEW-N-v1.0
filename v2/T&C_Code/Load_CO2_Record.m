function Ca = Load_CO2_Record(Date, Input_Directory)
%LOAD_CO2_RECORD  Atmospheric CO2 concentration [ppm] from the Ca_Data.mat record.
%
%   Ca = LOAD_CO2_RECORD(Date, Input_Directory) returns the atmospheric CO2
%   mixing ratio [ppm] interpolated onto the simulation timestamps Date
%   (MATLAB serial date numbers) from the hourly record in Ca_Data.mat
%   (variables: Ca, Date_CO2). Output matches the shape of Date.
%
%   The record covers 1975-2022; timestamps outside that span return NaN
%   (no extrapolation). The record is cached on first call so repeated
%   spin-up cycles do not reload the file.

persistent Date_CO2 Ca_record
if isempty(Date_CO2)
    S = load(fullfile(Input_Directory,'Ca_Data.mat'),'Date_CO2','Ca');
    Date_CO2  = S.Date_CO2(:);
    Ca_record = S.Ca(:);
end

Ca = interp1(Date_CO2, Ca_record, Date);
end
