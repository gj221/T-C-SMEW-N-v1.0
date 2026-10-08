function [B_IO_day,dep,dep_erw] = Deposition_Day(B_IO,Date,precip_mm,n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The atmospheric deposition of one biogeochemistry day [g m-2 d-1].
%%
%% Without B_IO.Dep_events the constant rates B_IO.DepN ... DepSi stand and
%% B_IO is returned unchanged.  With it (Build_Site_Deposition) N, K, Ca and
%% Mg arrive two ways:
%%   wet  the year's precipitation-weighted mean concentration times the
%%        day's precipitation, rain or snow, so it comes with the storms
%%        and wet years deposit more;
%%   dry  the year's dry deposition (gases and particles), spread evenly
%%        over the days of the year, rain or shine.
%% P and Si keep their constant rates.  Date is a datevec row (its year
%% selects the record's row; a year the record lacks is an error, never a
%% silent zero); precip_mm is the day's precipitation [mm].  n_days
%% (default 1) gives the mean daily rate of a span of n_days whose total
%% precipitation is precip_mm.
%%
%% dep = [DepN DepP DepK DepCa DepMg DepSi] of the day.
%% dep_erw = [Na An] of the day for the ERW chemistry [umol, umolc m-2 d-1]:
%% Na and the strong-acid anions (SO4 + Cl), wet and dry the same way, when
%% the record carries them (erw_elements); otherwise [0 0].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_IO_day = B_IO;
dep_erw = [0 0];
if isfield(B_IO,'Dep_events') && ~isempty(B_IO.Dep_events)
    if nargin < 4 || isempty(n_days); n_days = 1; end
    if ~(isnumeric(precip_mm) && isscalar(precip_mm) && isfinite(precip_mm) && ...
            precip_mm >= 0 && isnumeric(n_days) && isscalar(n_days) && ...
            isfinite(n_days) && n_days > 0)
        error('Deposition_Day:BadPrecipitation', ...
            'The precipitation must be a finite nonnegative total over a positive span.');
    end
    events = B_IO.Dep_events;
    row = find(events.year == Date(1),1);
    if isempty(row)
        error('Deposition_Day:NoYear', ...
            'The deposition record has no row for %d (it covers %d-%d).', ...
            Date(1),min(events.year),max(events.year));
    end
    rate = events.dry_g_m2_d(row,:)+events.wet_g_m2_per_mm(row,:)*precip_mm/n_days;
    B_IO_day.DepN = rate(1);
    B_IO_day.DepK = rate(2);
    B_IO_day.DepCa = rate(3);
    B_IO_day.DepMg = rate(4);
    if isfield(events,'erw_elements')
        dep_erw = events.erw_dry_umol_m2_d(row,:)+ ...
            events.erw_wet_umol_m2_per_mm(row,:)*precip_mm/n_days;
    end
end
%% Each rate to double first: the global maps hold P and Si in single, and
%% concatenating would round the whole record to single.
dep = [double(B_IO_day.DepN) double(B_IO_day.DepP) double(B_IO_day.DepK) ...
    double(B_IO_day.DepCa) double(B_IO_day.DepMg) double(B_IO_day.DepSi)];

end
