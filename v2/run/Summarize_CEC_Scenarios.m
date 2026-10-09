function [summary,daily] = Summarize_CEC_Scenarios(files,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The CEC scenarios' plot-scale results (Run_Plot_CEC_Scenarios), side by
%% side: daily means of the ERW chemistry's bulk and pore-water pH, base
%% saturation and exchangeable Al and H, and per scenario its start, end
%% and change.
%%
%% files   the scenarios' result files
%% Name/value
%%   'Output'  folder to write CEC_scenarios_summary.csv and
%%             CEC_scenarios_daily.csv to (default: none written)
%% summary  one row per scenario: CEC [cmolc/kg] and the first and last
%%          day's bulk pH, pore pH, base saturation, f_Al and f_H, and the
%%          bulk pH's change over the run
%% daily    one row per scenario and day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('Output','');
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k})
        error('Summarize_CEC_Scenarios:BadOptions','Unknown option %s.',varargin{k});
    end
    opt.(varargin{k}) = varargin{k+1};
end
files = cellstr(files);
rows = cell(numel(files),1);
days_tables = cell(numel(files),1);
for s = 1:numel(files)
    R = load(files{s},'ERW_biogeochem','ERW_PARAM_input','transient_start','t');
    b = R.ERW_biogeochem;
    cec = R.ERW_PARAM_input.CEC_tot_in;
    day = max(1,ceil(R.t(:) - 1e-9));   % t: the end of each ERW step [d]
    series = [b.pH_bulk(:) b.pH(:) b.f_Ca(:)+b.f_Mg(:)+b.f_K(:)+b.f_Na(:) b.f_Al(:) b.f_H(:)];
    nd = max(day);
    d = zeros(nd,size(series,2));
    for c = 1:size(series,2)
        ok = isfinite(series(:,c));   % a series the ERW reports only at some steps
        d(:,c) = accumarray(day(ok),series(ok,c),[nd 1])./accumarray(day(ok),1,[nd 1]);
    end
    date = datetime(R.transient_start,'ConvertFrom','datenum') + days(0:size(d,1)-1)';
    days_tables{s} = table(repmat(cec,size(d,1),1),date,d(:,1),d(:,2),d(:,3),d(:,4),d(:,5), ...
        'VariableNames',{'CEC_cmolc_kg','date','bulk_pH','pore_pH','base_saturation','f_Al','f_H'});
    rows{s} = table(cec,d(1,1),d(end,1),d(end,1)-d(1,1),d(1,2),d(end,2),d(1,3),d(end,3), ...
        d(end,4),d(end,5),size(d,1),'VariableNames',{'CEC_cmolc_kg','bulk_pH_start', ...
        'bulk_pH_end','bulk_pH_change','pore_pH_start','pore_pH_end','base_saturation_start', ...
        'base_saturation_end','f_Al_end','f_H_end','days'});
end
summary = sortrows(vertcat(rows{:}),'CEC_cmolc_kg');
daily = vertcat(days_tables{:});
if ~isempty(opt.Output)
    if ~isfolder(opt.Output), mkdir(opt.Output); end
    writetable(summary,fullfile(opt.Output,'CEC_scenarios_summary.csv'));
    writetable(daily,fullfile(opt.Output,'CEC_scenarios_daily.csv'));
end
end
