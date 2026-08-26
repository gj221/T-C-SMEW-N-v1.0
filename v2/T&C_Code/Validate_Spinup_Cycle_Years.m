function Validate_Spinup_Cycle_Years(Date, cycle_years)
% Require complete, strictly hourly meteorological years for cyclic spin-up.

valid_dates = isnumeric(Date) && isreal(Date) && isvector(Date) && ...
    ~isempty(Date) && all(isfinite(Date));
valid_years = isnumeric(cycle_years) && isreal(cycle_years) && ...
    isvector(cycle_years) && ~isempty(cycle_years) && ...
    all(isfinite(cycle_years)) && all(cycle_years == floor(cycle_years)) && ...
    numel(unique(cycle_years)) == numel(cycle_years);
if ~(valid_dates && valid_years)
    error('TC_SMEW_N_run:SpinupYears', ...
        'Spin-up dates and cycle years must be finite numeric vectors with unique integer years.');
end

Date = Date(:);
cycle_years = cycle_years(:).';
[forcing_years, ~, ~] = datevec(Date);
hourly_tolerance_d = 1/(24*60*60); % one second in serial-date units

for year = cycle_years
    year_start = datenum(year, 1, 1, 0, 0, 0);
    next_year = datenum(year + 1, 1, 1, 0, 0, 0);
    expected_count = 24*(next_year - year_start);
    year_dates = Date(forcing_years == year);
    complete_endpoints = numel(year_dates) == expected_count && ...
        abs(year_dates(1) - year_start) <= hourly_tolerance_d && ...
        abs(year_dates(end) - (next_year - 1/24)) <= hourly_tolerance_d;
    if complete_endpoints
        strictly_hourly = all(abs(diff(year_dates) - 1/24) <= hourly_tolerance_d);
    else
        strictly_hourly = false;
    end
    if ~(complete_endpoints && strictly_hourly)
        error('TC_SMEW_N_run:SpinupYears', ...
            ['Meteorological year %d is not a complete hourly year from ' ...
            'January 1 00:00 through December 31 23:00 (%d records expected; %d found).'], ...
            year, expected_count, numel(year_dates));
    end
end
end
