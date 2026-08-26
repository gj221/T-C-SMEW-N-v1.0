function values = Spinup_Fertilizer_Series(source, day_of_year, cycle_year, variable_name)
if ~(isnumeric(cycle_year) && isreal(cycle_year) && isscalar(cycle_year) && ...
        isfinite(cycle_year) && cycle_year == floor(cycle_year))
    error('Spinup_Fertilizer_Series:InvalidCalendar', ...
        'The fertilizer ledger needs one finite integer year.');
end
if ~isnumeric(day_of_year) || ~isreal(day_of_year) || isempty(day_of_year)
    error('Spinup_Fertilizer_Series:InvalidCalendar', ...
        'The fertilizer ledger needs at least one Julian day.');
end
day_of_year = day_of_year(:);
n_days = days(datetime(cycle_year+1,1,1) - datetime(cycle_year,1,1));
if any(~isfinite(day_of_year)) || any(day_of_year < 1) || ...
        any(day_of_year ~= floor(day_of_year)) || any(day_of_year > n_days)
    error('Spinup_Fertilizer_Series:InvalidCalendar', ...
        'The fertilizer ledger contains a Julian day outside the model year.');
end

if ~istimetable(source)
    if ~isnumeric(source) || ~isreal(source) || ~isvector(source) || isempty(source)
        error('Spinup_Fertilizer_Series:InvalidValues', ...
            'The fertilizer schedule must be a numeric vector.');
    end
    if max(day_of_year) > numel(source)
        error('Spinup_Fertilizer_Series:ScheduleCoverage', ...
            'The fertilizer schedule does not cover the simulated year.');
    end
    values = source(day_of_year);
else
    if ~ismember(variable_name, source.Properties.VariableNames)
        error('Spinup_Fertilizer_Series:MissingVariable', ...
            'The fertilizer timetable is missing %s.', variable_name);
    end
    source_days = dateshift(source.Properties.RowTimes, 'start', 'day');
    if numel(unique(source_days)) ~= numel(source_days)
        error('Spinup_Fertilizer_Series:DuplicateDates', ...
            'The fertilizer timetable contains duplicate dates.');
    end
    source_values = source.(variable_name);
    if ~isnumeric(source_values) || ~isreal(source_values) || ...
            ~isvector(source_values) || numel(source_values) ~= height(source)
        error('Spinup_Fertilizer_Series:InvalidValues', ...
            'The fertilizer timetable variable must be a numeric vector.');
    end
    model_days = datetime(cycle_year,1,1) + days(day_of_year-1);
    [matched, location] = ismember(model_days, source_days);
    values = zeros(size(day_of_year));
    values(matched) = source_values(location(matched));
end
values = values(:);
if any(~isfinite(values)) || any(values < 0)
    error('Spinup_Fertilizer_Series:InvalidValues', ...
        'The fertilizer schedule must contain finite nonnegative values.');
end
end
