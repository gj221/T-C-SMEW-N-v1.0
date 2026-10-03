function [cycle_years,handoff_start_dn,handoff_end_dn,avail_years] = ...
    Spinup_Cycle_Years(met_dates,transient_start)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spinup_Cycle_Years
%%
%% The calendar years a spin-up cycles (Run_BG_Spinup's meteorological
%% ring) and the partial-year handoff extension up to the transient start,
%% from the hours of the spin-up forcing.  One definition, shared by the
%% plot-scale spin-up (Run_Automatic_Site_Spinup) and the flux-tracking
%% field run (Run_Field_Lateral_Tracking), whose segments must be exactly
%% the ring the replayed cells cycle.
%%
%%   met_dates        datenums of the spin-up forcing hours
%%   transient_start  datenum of the transient start (NaN: none)
%%
%% A transient year inside the forcing ends the ring at the year before it
%% and adds the extension [1 Jan, transient_start); otherwise every forcing
%% year cycles and there is no extension (NaN bounds).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yr,~,~] = datevec(met_dates(:));
avail_years = unique(yr)';
if isempty(avail_years)
    error('TC_SMEW_N_run:SpinupYears', 'Spin-up forcing contains no meteorological years.');
end
cycle_years = avail_years;
handoff_start_dn = NaN;
handoff_end_dn = NaN;
if ~isnan(transient_start) && isfinite(transient_start)
    tmp = datevec(transient_start - 1/24);
    transient_year = tmp(1);
    if ismember(transient_year, avail_years)
        cycle_years = avail_years(avail_years < transient_year);
        handoff_start_dn = datenum(transient_year, 1, 1);
        handoff_end_dn = transient_start;
        if isempty(cycle_years)
            error('TC_SMEW_N_run:SpinupHandoff', ...
                'No complete spin-up years are available before transient year %d.', transient_year);
        end
    else
        warning('TC_SMEW_N_run:SpinupHandoff', ...
            'Transient year %d is outside spin-up range [%d-%d]; stopping on convergence alone.', ...
            transient_year, avail_years(1), avail_years(end));
    end
end
if any(diff(cycle_years) ~= 1)
    error('TC_SMEW_N_run:SpinupYears', 'Spin-up cycle years must be contiguous complete years.');
end
Validate_Spinup_Cycle_Years(met_dates, cycle_years);
end
