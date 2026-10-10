function result_file = Run_Field_Transient(shape_file,agronomic_file,state_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Field_Transient
%%
%% The field transient: every cell of one field polygon goes forward from
%% its spun-up state (Run_Field_Cell_Spinup's CellBGStates.mat) over the
%% site's transient window, in 2-D (MAIN_FRAME_BG_ERW_2D, with lateral flow
%% over the DEM's terrain), with the ERW chemistry on and no rock. For
%% Hamden: bulk soil pH under atmospheric deposition alone, 2020-01-01 to
%% 2026-09-29.
%%
%% CELLS. Each cell is configured as its spin-up configured it
%% (Run_Field_Cell_Spinup 'DryRun': the workbook case, the site's MOD and
%% ERW columns, the cell's own soil and cover) and starts from the state its
%% spin-up handed off at the window's first hour (Load_Field_Spinup_States
%% refuses a state for another instant or other inputs). The forcing, CO2
%% and rain records are the ones the spin-up's handoff check read: the
%% rows' meteo_file (Load_Spinup_Meteo), co2_file (Load_CO2_Record) and
%% rain_chem_file (Load_Rain_Chemistry). A row with rock is refused: this
%% run adds none.
%%
%% SEGMENTS. The 2-D driver holds every column's hourly state, and each ERW
%% column a ten-minute record, for the whole window it is given; 936
%% columns over seven years would need hundreds of GB. The window runs as
%% segments of SegmentDays, each going on from the one before exactly as
%% one run would (Field_Transient_Carry), and each with the whole window's
%% longest day, which the parameter files would otherwise find over the
%% segment's days: segmented and single runs agree bit for bit
%% (test_field_transient).
%%
%% CHECKPOINT AND RESUME. A finished segment is written to
%% <OutputDirectory>/segment_<k>.mat (its records) and carry_<k>.mat (the
%% state the next starts from) under a key that chains the run's identity
%% -- the cells and their rows, the spun-up states, the terrain, the window,
%% the segment plan and each segment's forcing -- through every earlier
%% segment; a rerun reuses it. A segment written by other code is refused
%% unless ForceResume (env TC_CKPT_FORCE_RESUME=1) accepts it. Inside a
%% segment the driver's rolling checkpoint (checkpoints/, every
%% CheckpointMinutes) resumes. A job that runs out of time is simply run
%% again.
%%
%% OUTPUT. <OutputDirectory>/Field_Transient_<from>_<to>.mat holds
%% `result`, one record per ERW day of every cell (rows) and cell
%% (columns), and the field's area-weighted means; the means are also
%% written to a CSV beside it. The model day runs from 01:00 to 01:00 local
%% standard time (the drivers' daily step). Record 1 is the state at 00:00
%% of the first day (the handoff state as the ERW chemistry starts from
%% it); record 2 is 01:00 of the first day; each later record is 01:00 of a
%% day, with the fluxes over the 24 hours before it; the last is the end of
%% the window (00:00 after its last day), whose 23 hours run no daily step
%% of the soil biogeochemistry (its T&C leaching is NaN). States: bulk and
%% pore-water pH, base saturation and the exchanger's composition (charge
%% fractions of the CEC), the CEC, the dissolved Ca, Mg, K, Na, free Al3+
%% and total dissolved Al, the residual strong anion, DIC and the water
%% content of the biogeochemical layer. Fluxes: leaching of Ca, Mg and K
%% (and NO3, NH4) by the soil biogeochemistry, of Na, Al, strong anion and
%% DIC by the ERW chemistry, drainage out of the column and precipitation.
%% result.units names the units.
%%
%% Usage
%%   Run_Field_Transient(shp,agro,states,'OutputDirectory',D, ...
%%       'MODParamColumns',M,'ERWParamColumns',E,'GridCellSizeM',5, ...
%%       'MinCoverage','majority','FieldDEM',dem)
%%   ... ,'NDays',3,'CellIds',[33 34 35 62])   % smoke test
%% On Bouchet: Cluster_Field_Spinup_Task('transient') (v2/cluster).
%%
%% Options
%%   OutputDirectory   required
%%   CaseNumber        workbook case the site borrows (5)
%%   MODParamColumns, ERWParamColumns   structs or one-row CSV files, as the
%%                     spin-up took them
%%   FieldFeatureIndices, FieldPolygonIndex, GridCellSizeM, MinCoverage
%%                     the grid, as the spin-up built it
%%   FieldDEM          DEM directory or GeoTIFF (empty: a flat field)
%%   CEC               a CEC scenario: every cell's CEC at pH 7 [cmolc/kg]
%%                     (empty: each cell's own). The spin-up does not
%%                     depend on it (its ERW chemistry is off), so every
%%                     scenario starts from the one spin-up; each keeps the
%%                     handoff's exchanger fractions at its own capacity,
%%                     re-partitioned by the seed at the pore pH that
%%                     reproduces the bulk pH, as Run_Plot_CEC_Scenarios
%%                     does for one column. The result is tagged _CEC_<value>.
%%   SegmentDays       days per segment (30; at least 2, the days the next
%%                     segment repeats)
%%   NDays             only the window's first NDays (a smoke test)
%%   CellIds           only these cells, as a field of their own (a smoke
%%                     test); the states file may hold every cell
%%   ParallelWorkers   thread workers for the columns' hourly step and daily
%%                     ERW chemistry (0; env TC_PARALLEL_WORKERS)
%%   CheckpointMinutes the driver's rolling checkpoint (60; env
%%                     TC_TRANSIENT_CKPT_MINUTES; 0 disables)
%%   ForceResume       reuse segments and checkpoints of other code (false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('OutputDirectory','','CaseNumber',5,'MODParamColumns',[], ...
    'ERWParamColumns',[],'FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'GridCellSizeM',[],'MinCoverage',[],'FieldDEM','','CEC',[],'SegmentDays',30, ...
    'NDays',[],'CellIds',[],'ParallelWorkers',[],'CheckpointMinutes',[], ...
    'ForceResume',[]);
if mod(numel(varargin),2) ~= 0
    error('Run_Field_Transient:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_Field_Transient:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
shape_file = existing_file(shape_file,'shape_file');
agronomic_file = existing_file(agronomic_file,'agronomic_file');
state_file = existing_file(state_file,'state_file');
if isempty(opt.OutputDirectory)
    error('Run_Field_Transient:BadOptions','OutputDirectory is required.');
end
out_dir = char(opt.OutputDirectory);
if ~exist(out_dir,'dir'); mkdir(out_dir); end
out_dir = canonical(out_dir);
positive_integer(opt.CaseNumber,'CaseNumber');
positive_integer(opt.SegmentDays,'SegmentDays');
overlap = Field_Transient_Carry('overlap_days');
if opt.SegmentDays < overlap
    error('Run_Field_Transient:BadOptions', ...
        'SegmentDays must be at least %d: the next segment repeats that many days.',overlap);
end
if ~isempty(opt.CEC) && ~(isnumeric(opt.CEC) && isscalar(opt.CEC) && isfinite(opt.CEC) && opt.CEC > 0)
    error('Run_Field_Transient:BadOptions','CEC must be one positive value [cmolc/kg].');
end
if ~isempty(opt.NDays); positive_integer(opt.NDays,'NDays'); end
if ~isempty(opt.CellIds)
    ids = opt.CellIds(:);
    if ~(isnumeric(ids) && all(isfinite(ids)) && all(ids == fix(ids)) && ...
            all(ids >= 1) && numel(unique(ids)) == numel(ids))
        error('Run_Field_Transient:BadOptions','CellIds must be distinct cell ids.');
    end
end
workers = opt.ParallelWorkers;
if isempty(workers); workers = env_number('TC_PARALLEL_WORKERS',0); end
ckpt_minutes = opt.CheckpointMinutes;
if isempty(ckpt_minutes); ckpt_minutes = env_number('TC_TRANSIENT_CKPT_MINUTES',60); end
if ~(isnumeric(ckpt_minutes) && isscalar(ckpt_minutes) && isfinite(ckpt_minutes) && ckpt_minutes >= 0)
    error('Run_Field_Transient:BadOptions','CheckpointMinutes must be a nonnegative number.');
end
force_resume = opt.ForceResume;
if isempty(force_resume)
    force_resume = strcmpi(strtrim(getenv('TC_CKPT_FORCE_RESUME')),'1');
end
force_resume = logical(force_resume);

here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
Directory = fullfile(v2,'T&C_Code');
Input_Directory = fullfile(Directory,'Inputs');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory,Input_Directory,Parameter_Directory,fullfile(v2,'spinup'),here);
cols = column_struct(opt.MODParamColumns,'MODParamColumns');
erw_cols = column_struct(opt.ERWParamColumns,'ERWParamColumns');

fprintf('\n=== Run_Field_Transient: %s ===\n',shape_file);
%% ---- the cells, as the spin-up configured them ------------------------
grid_args = {'FieldFeatureIndices',opt.FieldFeatureIndices, ...
    'FieldPolygonIndex',opt.FieldPolygonIndex,'GridCellSizeM',opt.GridCellSizeM, ...
    'MinCoverage',opt.MinCoverage};
[~,CellInputs] = Run_Field_Cell_Spinup(shape_file,agronomic_file,grid_args{:}, ...
    'CaseNumber',opt.CaseNumber,'MODParamColumns',cols,'ERWParamColumns',erw_cols, ...
    'OutputFile',fullfile(out_dir,'CellBGStates_not_written.mat'),'DryRun',true);
setups = Build_Shapefile_Grid_Setups(shape_file,opt.GridCellSizeM, ...
    opt.FieldFeatureIndices,opt.MinCoverage);
poly = opt.FieldPolygonIndex;
if isempty(poly); poly = 1; end
topo = setups(poly);
clear setups
if ~isequal(topo.cell_id(:),[CellInputs.cell_id].')
    error('Run_Field_Transient:GridMismatch', ...
        'The grid and the spin-up configuration hold different cells.');
end
n_field = topo.N_cells;
if ~isempty(opt.CellIds)
    [found,idx] = ismember(opt.CellIds(:),topo.cell_id(:));
    if ~all(found)
        error('Run_Field_Transient:UnknownCell','The grid holds no cell %s.', ...
            mat2str(opt.CellIds(~found).'));
    end
    idx = sort(idx);
    topo = subset_topology(topo,idx);
    CellInputs = CellInputs(idx);
    fprintf('Cells: %d of the field''s %d (%s)\n',topo.N_cells,n_field,mat2str(topo.cell_id(:).'));
end
N_cells = topo.N_cells;
if ~isempty(opt.FieldDEM)
    dem = Load_Field_DEM(char(opt.FieldDEM),topo,'OutResolutionM',dem_resolution(opt.GridCellSizeM,topo));
    topo = Attach_Topography_To_Topology(topo,dem);
    clear dem
    fprintf('Topography: relief %.1f m, slope %.1f-%.1f%%, %d outlet cell(s)\n', ...
        max(topo.cell_z_m)-min(topo.cell_z_m),100*min(topo.cell_slope), ...
        100*max(topo.cell_slope),sum(topo.cell_is_outlet));
else
    fprintf('Topography: none (flat field)\n');
end

%% ---- the rows every cell shares, the window ----------------------------
MOD1 = CellInputs(1).MOD_PARAM;
shared = {'meteo_file','co2_file','rain_chem_file','transient_start','transient_end', ...
    'Lat','Lon','Zs','Zbio','MOD_PARAM','t_res','Location','id_location'};
for c = 2:N_cells
    for k = 1:numel(shared)
        if ~isequaln(row_value(MOD1,shared{k}),row_value(CellInputs(c).MOD_PARAM,shared{k}))
            error('Run_Field_Transient:RowMismatch', ...
                'Cell %d has a different %s from cell %d; the field runs on one forcing.', ...
                CellInputs(c).cell_id,shared{k},CellInputs(1).cell_id);
        end
    end
end
for c = 1:N_cells
    E = CellInputs(c).ERW_PARAM;
    rock = [row_number(E,'M_rock_in') row_number(E,'M_rock_event')];
    if any(rock(isfinite(rock)) > 0)
        error('Run_Field_Transient:Rock', ...
            ['Cell %d''s ERW row adds rock (M_rock_in %g, M_rock_event %g); the field ' ...
             'transient adds none: give the site a no-rock ERW row.'], ...
            CellInputs(c).cell_id,rock(1),rock(2));
    end
end
start_dn = row_datenum(MOD1,'transient_start');
end_dn = row_datenum(MOD1,'transient_end');
if ~(isfinite(start_dn) && isfinite(end_dn) && end_dn >= start_dn)
    error('Run_Field_Transient:NoWindow', ...
        'The rows must carry transient_start and transient_end.');
end
n_window = end_dn - start_dn + 1;
n_days_row = row_number(MOD1,'NumberOfDays');
if isempty(opt.NDays) && isfinite(n_days_row) && n_days_row ~= n_window
    error('Run_Field_Transient:NumberOfDays', ...
        'NumberOfDays %d does not match the window %s to %s (%d days).',n_days_row, ...
        datestr(start_dn,'yyyy-mm-dd'),datestr(end_dn,'yyyy-mm-dd'),n_window);
end
if ~isempty(opt.NDays)
    if opt.NDays > n_window
        error('Run_Field_Transient:BadOptions','NDays %d exceeds the %d-day window.', ...
            opt.NDays,n_window);
    end
    end_dn = start_dn + opt.NDays - 1;
end
days = (start_dn:end_dn)';
n_days = numel(days);
Lat = MOD1.Lat; Lon = MOD1.Lon; DeltaGMT = round(Lon/15);
id_location = char(string(MOD1.Location));
if ~isempty(row_value(MOD1,'id_location')) && strlength(strtrim(string(row_value(MOD1,'id_location')))) > 0
    id_location = char(string(row_value(MOD1,'id_location')));
end
meteo_name = row_value(MOD1,'meteo_file');
if isempty(meteo_name) || strlength(strtrim(string(meteo_name))) == 0
    error('Run_Field_Transient:NoMeteo', ...
        'The rows must name the meteo_file the spin-up''s handoff read.');
end
meteo_path = meteo_file_path(char(string(meteo_name)),fullfile(v2,'Data'),Input_Directory);
co2_file = 'Ca_Data.mat';
if ~isempty(row_value(MOD1,'co2_file')) && strlength(strtrim(string(row_value(MOD1,'co2_file')))) > 0
    co2_file = char(strtrim(string(row_value(MOD1,'co2_file'))));
end
fprintf('Window: %s to %s (%d days), %d segment(s) of up to %d days\n', ...
    datestr(start_dn,'yyyy-mm-dd'),datestr(end_dn,'yyyy-mm-dd'),n_days, ...
    ceil(n_days/opt.SegmentDays),opt.SegmentDays);

%% ---- the spun-up states -------------------------------------------------
states_path = state_file;
if ~isempty(opt.CellIds)
    %% Load_Field_Spinup_States reads a file of exactly the run's cells.
    S = load(state_file,'CellBGStates');
    all_ids = arrayfun(@(s) double(s.cell_id),S.CellBGStates);
    [found,where] = ismember(topo.cell_id(:),all_ids(:));
    if ~all(found)
        error('Run_Field_Transient:StateMissing','%s holds no state for cell(s) %s.', ...
            state_file,mat2str(topo.cell_id(~found).'));
    end
    CellBGStates = S.CellBGStates(where);
    clear S
    %% Rewritten only when it changes: the driver's checkpoints are keyed
    %% on its bytes.
    states_path = fullfile(out_dir,'CellBGStates_cells.mat');
    same = false;
    if isfile(states_path)
        P = load(states_path,'CellBGStates');
        same = isfield(P,'CellBGStates') && strcmp(Field_Transient_Carry('digest', ...
            P.CellBGStates),Field_Transient_Carry('digest',CellBGStates));
        clear P
    end
    if ~same
        save_atomic(states_path,struct('CellBGStates',CellBGStates));
    end
    clear CellBGStates
end
states = Load_Field_Spinup_States(states_path,topo,{CellInputs.agronomic_identity},start_dn);
for c = 1:N_cells
    CellInputs(c).BG_STATE = states(c);
end
fprintf('States: %d cells handed off at %s (%s)\n',N_cells, ...
    datestr(start_dn,'yyyy-mm-dd HH:MM'),state_file);
if ~isempty(opt.CEC)
    %% after the states are matched to the inputs they were spun up with
    for c = 1:N_cells
        CellInputs(c).ERW_PARAM.CEC_tot_in = opt.CEC;
    end
    fprintf(['CEC scenario: %g cmolc/kg in every cell, each with its handoff''s ' ...
        'exchanger fractions\n'],opt.CEC);
end

%% ---- forcing, CO2 and rain records over the window ---------------------
[met,~] = Load_Spinup_Meteo(meteo_path,start_dn,end_dn,Lat,Lon,DeltaGMT);
met_day = floor(met.Date + 1/86400);
Ca_all = Load_CO2_Record(met.Date,Input_Directory,co2_file);
rain_chemistry = Load_Rain_Chemistry(MOD1,Input_Directory,id_location);
%% The run scripts' solar window, and the window's longest day, which every
%% segment takes (MAIN_FRAME_BG_ERW_2D, FIELD_TRANSIENT_LMAX_DAY).
t_bef = -0.67; t_aft = 1.67;
Lmax_window = longest_day(met.Date,Lat,Lon,DeltaGMT,t_bef,t_aft);

%% ---- identity, segment plan ---------------------------------------------
code_identity = code_fingerprint(v2);
run_key = sha256_hex(sprintf(['field_transient|case=%d|cells=%s|states=%s|terrain=%s|' ...
    'window=%.17g:%.17g|segment_days=%d|overlap_days=%d'], opt.CaseNumber, ...
    cell_identity(CellInputs),Field_Transient_Carry('digest',states), ...
    terrain_identity(topo),start_dn,end_dn,opt.SegmentDays,overlap));
seg_first = (1:opt.SegmentDays:n_days)';
seg_last = min(seg_first + opt.SegmentDays - 1,n_days);
n_seg = numel(seg_first);
ckpt_dir = fullfile(out_dir,'checkpoints');
if ~exist(ckpt_dir,'dir'); mkdir(ckpt_dir); end
dependencies = {agronomic_file,states_path};
workers = min(workers,N_cells);
if workers < 2; workers = 0; end
parallel_pool_started = start_thread_pool(workers);
cleanup_pool = onCleanup(@() stop_thread_pool(parallel_pool_started));

ctx = struct('topo',topo,'CellInputs',CellInputs,'Directory',Directory, ...
    'Input_Directory',Input_Directory,'Parameter_Directory',Parameter_Directory, ...
    'case_number',opt.CaseNumber,'rain',rain_chemistry,'Lat',Lat,'Lon',Lon, ...
    'DeltaGMT',DeltaGMT,'t_bef',t_bef,'t_aft',t_aft,'id_location',id_location,'workers',workers, ...
    'ckpt_minutes',ckpt_minutes,'force_resume',force_resume);
segments = cell(n_seg,1);
key = run_key;
t_run = tic;
for s = 1:n_seg
    previous_key = key;
    own = days(seg_first(s):seg_last(s));
    forcing_days = own;
    if s > 1; forcing_days = [own(1)-(overlap:-1:1)'; own]; end
    mask = met_day >= forcing_days(1) & met_day <= forcing_days(end);
    seg_met = slice_met(met,mask);
    seg_Ca = Ca_all(mask);
    if numel(seg_met.Date) ~= 24*numel(forcing_days)
        error('Run_Field_Transient:Forcing','Segment %d has %d forcing hours for %d days.', ...
            s,numel(seg_met.Date),numel(forcing_days));
    end
    key = sha256_hex(sprintf('%s|segment=%d|days=%.17g:%.17g|forcing=%s',key,s, ...
        own(1),own(end),Field_Transient_Carry('digest',{seg_met,seg_Ca})));
    seg_file = fullfile(out_dir,sprintf('segment_%03d.mat',s));
    carry_file = fullfile(out_dir,sprintf('carry_%03d.mat',s));
    if isfile(seg_file)
        P = load(seg_file,'key','code_identity');
        if isfield(P,'key') && strcmp(P.key,key)
            if ~strcmp(P.code_identity,code_identity)
                if ~force_resume
                    error('Run_Field_Transient:CodeChanged', ...
                        ['%s was written by other code. Give a new OutputDirectory, or set ' ...
                         'ForceResume (TC_CKPT_FORCE_RESUME=1) to splice the two deliberately.'], ...
                        seg_file);
                end
                warning('Run_Field_Transient:ForcedResume', ...
                    'Reusing %s, written by other code (ForceResume).',seg_file);
            end
            P = load(seg_file,'records','info');
            segments{s} = P.records;
            fprintf('Segment %d/%d (%s to %s): reusing %s\n',s,n_seg, ...
                datestr(own(1),'yyyy-mm-dd'),datestr(own(end),'yyyy-mm-dd'),seg_file);
            continue
        end
        fprintf('Segment %d/%d: %s belongs to another run; running it afresh\n',s,n_seg,seg_file);
    end
    carry_in = [];
    seg_dependencies = dependencies;
    if s > 1
        previous_carry = fullfile(out_dir,sprintf('carry_%03d.mat',s-1));
        carried_key = '';
        if isfile(previous_carry)
            C = load(previous_carry,'key');
            carried_key = C.key;
        end
        if ~strcmp(carried_key,previous_key)
            error('Run_Field_Transient:CarryMissing', ...
                ['The state segment %d hands on (%s) is missing or belongs to another ' ...
                 'run: delete segment_%03d.mat and the segments after it to run again ' ...
                 'from there.'],s-1,previous_carry,s-1);
        end
        C = load(previous_carry,'carry');
        carry_in = C.carry;
        clear C
        seg_dependencies{end+1} = previous_carry; %#ok<AGROW>
    end
    fprintf('\nSegment %d/%d: %s to %s (%d day(s)%s)\n',s,n_seg, ...
        datestr(own(1),'yyyy-mm-dd'),datestr(own(end),'yyyy-mm-dd'),numel(own), ...
        ternary(s > 1,sprintf(', after %d carried days',overlap),''));
    seg = struct('met',seg_met,'Ca',seg_Ca,'carry',carry_in,'continues',s < n_seg, ...
        'first_record_j',ternary(s == 1,1,overlap+1),'ckpt_file', ...
        ternary(ckpt_minutes > 0,fullfile(ckpt_dir,sprintf('segment_%03d.mat',s)),''), ...
        'dependencies',{seg_dependencies},'Lmax_day',struct('window',Lmax_window, ...
        'segment',longest_day(seg_met.Date,Lat,Lon,DeltaGMT,t_bef,t_aft)));
    t_seg = tic;
    [records,carry,info] = run_segment(ctx,seg);
    info.wall_s = toc(t_seg);
    if s < n_seg
        save_atomic(carry_file,struct('carry',carry,'key',key));
    end
    clear carry carry_in
    save_atomic(seg_file,struct('records',records,'info',info,'key',key, ...
        'code_identity',code_identity,'segment',s,'days',own));
    segments{s} = records;
    if ~isempty(seg.ckpt_file) && isfile(seg.ckpt_file); delete(seg.ckpt_file); end
    old_carry = fullfile(out_dir,sprintf('carry_%03d.mat',s-2));
    if s > 2 && isfile(old_carry); delete(old_carry); end
    remaining = (n_days - seg_last(s))/max(1,numel(own))*info.wall_s;
    fprintf('Segment %d/%d done in %.1f min; about %.1f h to go\n',s,n_seg, ...
        info.wall_s/60,remaining/3600);
end
fprintf('\nAll %d segment(s) in %.1f min\n',n_seg,toc(t_run)/60);

%% ---- the result ---------------------------------------------------------
result = assemble(segments,topo,CellInputs);
result.provenance = struct('shape_file',shape_file,'agronomic_file',agronomic_file, ...
    'state_file',state_file,'case',opt.CaseNumber,'meteo_file',meteo_path, ...
    'co2_file',co2_file,'rain_chem_file',char(rain_chemistry.source_file), ...
    'window',[start_dn end_dn],'segment_days',opt.SegmentDays,'n_segments',n_seg, ...
    'run_key',run_key,'last_segment_key',key,'code_identity',code_identity, ...
    'CEC_scenario',opt.CEC,'cell_ids_requested',opt.CellIds,'n_field_cells',n_field, ...
    'matlab',version,'computer',computer,'slurm_job_id',getenv('SLURM_JOB_ID'), ...
    'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
tag = sprintf('%s_%s',datestr(start_dn,'yyyymmdd'),datestr(end_dn,'yyyymmdd'));
if ~isempty(opt.CEC); tag = sprintf('%s_CEC_%s',tag,strrep(sprintf('%g',opt.CEC),'.','p')); end
if ~isempty(opt.CellIds); tag = sprintf('%s_cells%d',tag,N_cells); end
result_file = fullfile(out_dir,sprintf('Field_Transient_%s.mat',tag));
save_atomic(result_file,struct('result',result));
csv_file = fullfile(out_dir,sprintf('Field_Transient_%s_field_mean.csv',tag));
write_field_mean_csv(result,csv_file);
fprintf('Saved %s\n      %s\n',result_file,csv_file);
fprintf(['Field mean: bulk pH %.3f -> %.3f, pore pH %.3f -> %.3f, base saturation ' ...
    '%.3f -> %.3f\n'],result.field.pH_bulk(1),result.field.pH_bulk(end), ...
    result.field.pH_pore(1),result.field.pH_pore(end), ...
    result.field.base_saturation(1),result.field.base_saturation(end));
clear cleanup_pool
end

%% =======================================================================
function [ft_records,ft_carry,ft_info] = run_segment(ft_ctx,ft_seg)
%% One segment: the workspace MAIN_FRAME_BG_ERW_2D reads, built as
%% Run_One_Subfield and the spin-up's tracking run build it, then the
%% driver, then the segment's records and the state the next starts from.
%% Every name of this function's own starts with ft_: the driver's script
%% shares the workspace.
topo = ft_ctx.topo;
Directory = ft_ctx.Directory; %#ok<NASGU>
Input_Directory = ft_ctx.Input_Directory; %#ok<NASGU>
Parameter_Directory = ft_ctx.Parameter_Directory;
NN_case_study_ERW_MOD = ft_ctx.case_number; %#ok<NASGU>
Lat = ft_ctx.Lat; Lon = ft_ctx.Lon; DeltaGMT = ft_ctx.DeltaGMT; %#ok<NASGU>
id_location = ft_ctx.id_location; %#ok<NASGU>
Date = ft_seg.met.Date;  Pr = ft_seg.met.Pr;  Ta = ft_seg.met.Ta;
Ws = ft_seg.met.Ws;      ea = ft_seg.met.ea;  Pre = ft_seg.met.Pre; %#ok<NASGU>
SAD1 = ft_seg.met.SAD1;  SAD2 = ft_seg.met.SAD2; %#ok<NASGU>
SAB1 = ft_seg.met.SAB1;  SAB2 = ft_seg.met.SAB2; %#ok<NASGU>
N = ft_seg.met.N;        Tdew = ft_seg.met.Tdew;  esat = ft_seg.met.esat; %#ok<NASGU>
PARB = ft_seg.met.PARB;  PARD = ft_seg.met.PARD; %#ok<NASGU>
Ca = ft_seg.Ca;
NN = numel(Date);
NNd = ceil(NN/24) + 1;
Ws(Ws <= 0) = 0.01;
Ds = max(esat - ea,0); %#ok<NASGU>
[ft_YE,ft_MO,ft_DA,ft_HO] = datevec(Date);
Datam = [ft_YE(:) ft_MO(:) ft_DA(:) ft_HO(:)];
clear ft_YE ft_MO ft_DA ft_HO
if Datam(1,4) ~= 0
    error('Run_Field_Transient:SegmentStart','A segment must start at 00:00 (%s).', ...
        datestr(Date(1),'yyyy-mm-dd HH:MM'));
end
dt = 3600; dth = 1; Oa = 210000; t_bef = ft_ctx.t_bef; t_aft = ft_ctx.t_aft; %#ok<NASGU>
MOD_PARAM_input = ft_ctx.CellInputs(1).MOD_PARAM;
ERW_PARAM_input = ft_ctx.CellInputs(1).ERW_PARAM;
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'));
ms = length(Zs) - 1; %#ok<NASGU>
cc = MOD_PARAM_input.cc; %#ok<NASGU>
IrD = zeros(NN,1); %#ok<NASGU>
OPT_SoilBiogeochemistry = 1; %#ok<NASGU>
OPT_ERW = 1; %#ok<NASGU>
parallel_workers = ft_ctx.workers; %#ok<NASGU>
ERW_ckpt_file = ft_seg.ckpt_file; %#ok<NASGU>
ERW_ckpt_minutes = ft_ctx.ckpt_minutes; %#ok<NASGU>
ERW_ckpt_force_resume = ft_ctx.force_resume; %#ok<NASGU>
ERW_cell_input_dependencies = ft_seg.dependencies; %#ok<NASGU>
ft_base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
ColCfg = repmat(struct('PARAM_IC',ft_base_param,'PFT_OVERRIDE','', ...
    'apply_rock',false,'ERW_PARAM',ERW_PARAM_input,'MOD_PARAM',table(), ...
    'BG_STATE',struct(),'agronomic_identity','','appelo_alpha',[], ...
    'appelo_logKH',[],'cec_dpH',[],'KNmaxd',[],'last_application_year',[], ...
    'biogeo_overrides',[],'saturation_exponent',1,'RAIN_CHEM',ft_ctx.rain), ...
    topo.N_cells,1);
for ft_c = 1:topo.N_cells
    ft_in = ft_ctx.CellInputs(ft_c);
    ColCfg(ft_c).MOD_PARAM = ft_in.MOD_PARAM;
    ColCfg(ft_c).ERW_PARAM = ft_in.ERW_PARAM;
    ColCfg(ft_c).BG_STATE = ft_in.BG_STATE;
    ColCfg(ft_c).agronomic_identity = ft_in.agronomic_identity;
    if ismember('pft_override',ft_in.MOD_PARAM.Properties.VariableNames)
        ColCfg(ft_c).PFT_OVERRIDE = char(strtrim(string(ft_in.MOD_PARAM.pft_override)));
    end
end
clear ft_c ft_in ft_base_param
if ~isempty(ft_seg.carry)
    FIELD_TRANSIENT_CARRY = ft_seg.carry; %#ok<NASGU>
    ft_seg.carry = [];       % one copy for the driver's checkpoint to write
end
FIELD_TRANSIENT_CONTINUES = ft_seg.continues; %#ok<NASGU>
FIELD_TRANSIENT_LMAX_DAY = ft_seg.Lmax_day; %#ok<NASGU>
%% The driver's counters, declared: i and j the script alone creates would
%% be read below as the imaginary unit.
i = 1; j = 1; %#ok<NASGU>
ft_old_dir = cd(ft_ctx.Directory);
%% named cleanup: the driver's rolling checkpoint saves every variable but it
cleanup = onCleanup(@() cd(ft_old_dir));
Plant_Carbon_Deficit('reset');
MAIN_FRAME_BG_ERW_2D;
clear cleanup
ft_info = struct('carbon_deficit',Plant_Carbon_Deficit('report'), ...
    'computational_time_s',Computational_Time,'n_hours',NN,'n_columns',topo.N_cells, ...
    'date_first',Date(1),'date_last',Date(end),'parallel_workers',parallel_workers);
%% A segment that another follows leaves its last day's record to it: the
%% next ERW day re-anchors that day's end row (ERW_Daily_Step_2D), as one
%% run would before the record is read.
ft_records = segment_records(ft_seg.first_record_j,~ft_seg.continues,Date,Datam,NN, ...
    j - double(ft_seg.continues), ...
    ERW,PAR,LEAK_Ca,LEAK_Mg,LEAK_K,LEAK_NO3,LEAK_NH4,Lk,Pr,O);
%% The state the next segment starts from, or, at the end of the run, its
%% digest.
ft_names = Field_Transient_Carry('names');
ft_W = struct();
for ft_name = [ft_names.hourly ft_names.daily ft_names.column]
    ft_W.(ft_name{1}) = eval(ft_name{1});
end
ft_W.ERW = ERW; ft_W.lateral_diag = lateral_diag; ft_W.jDay = jDay; ft_W.L_day = L_day;
ft_W.Date = Date; ft_W.NN = NN; ft_W.NNd = NNd; ft_W.j = j;
ft_W.whos = whos;
ft_info.workspace_GB = sum([ft_W.whos(~strcmp({ft_W.whos.name},'ft_W')).bytes])/2^30;
fprintf('Segment workspace: %.2f GB (%d columns, %d hours)\n',ft_info.workspace_GB, ...
    topo.N_cells,NN);
if ft_seg.continues
    ft_carry = Field_Transient_Carry('capture',ft_W);
    ft_end = ft_carry;
else
    ft_carry = [];
    ft_end = Field_Transient_Carry('snapshot',ft_W);
end
ft_info.end_state_digest = Field_Transient_Carry('digest',ft_end);
ft_info.end_state_parts = Field_Transient_Carry('parts',ft_end);
end

%% =======================================================================
function R = segment_records(first_j,final,Date,Datam,NN,j_end,ERW,PAR, ...
    LEAK_Ca,LEAK_Mg,LEAK_K,LEAK_NO3,LEAK_NH4,Lk,Pr,O)
%% One record per ERW day from day row first_j to the segment's last (and,
%% at the end of the run, its final flush): the ERW state at the end of the
%% day's ten-minute window, the ERW leaching summed over it, and the soil
%% biogeochemistry's daily step (which reads the same 24 hours).
N_col = numel(ERW);
daily_rows = find(Datam(:,4) == 1);       % the daily step of day row j is at daily_rows(j-1)
jj = first_j:j_end;
n = numel(jj) + double(final);
res = ERW{1}.res;
for c = 2:N_col
    if ERW{c}.res ~= res
        error('Run_Field_Transient:ERWResolution','The ERW columns run on different steps.');
    end
end
state = {'pH_bulk','pH_pore','base_saturation','f_Ca','f_Mg','f_K','f_Na','f_Al','f_H', ...
    'CEC_umolc_m2','Ca_umol_L','Mg_umol_L','K_umol_L','Na_umol_L','Al3_umol_L', ...
    'Al_dissolved_umol_L','An_umolc_L','DIC_umol_L','theta_bio'};
flux = {'leach_Ca_g_m2','leach_Mg_g_m2','leach_K_g_m2','leach_NO3_gN_m2', ...
    'leach_NH4_gN_m2','leach_Na_umol_m2','leach_Al_umol_m2','leach_An_umolc_m2', ...
    'leach_DIC_umol_m2','drainage_mm'};
R = struct('time',nan(n,1),'hours',nan(n,1),'precipitation_mm',nan(n,1));
for k = [state flux]; R.(k{1}) = nan(n,N_col); end
for r = 1:n
    if r <= numel(jj) && jj(r) == 1
        hour_rows = [];
        R.time(r) = Date(1);
        dm_rows = []; dm_end = 1; t_row = 1;
    elseif r <= numel(jj)
        i_j = daily_rows(jj(r)-1);
        hour_rows = max(1,i_j-24):i_j-1;
        R.time(r) = Date(i_j);
        dm_first = (max(1,i_j-24)-1)*res + 1;
        dm_end = (i_j-1)*res;
        dm_rows = max(dm_first,2):dm_end;
        t_row = jj(r);
    else
        %% the final flush (ERW_Daily_Step_2D at i == NN)
        last_daily = daily_rows(end);
        hour_rows = last_daily:NN;
        R.time(r) = Date(NN) + 1/24;
        dm_first = (last_daily-1)*res + 1;
        dm_end = NN*res;
        dm_rows = max(dm_first,2):dm_end;
        t_row = [];
    end
    R.hours(r) = numel(hour_rows);
    if ~isempty(hour_rows); R.precipitation_mm(r) = sum(Pr(hour_rows)); end
    for c = 1:N_col
        b = ERW{c}.bg;
        scale = ERW{c}.const.Al_CEC_scale;
        f = [b.f_Ca(dm_end) b.f_Mg(dm_end) b.f_K(dm_end) b.f_Na(dm_end) b.f_Al(dm_end) b.f_H(dm_end)];
        R.pH_bulk(r,c) = b.pH_bulk(dm_end);
        R.pH_pore(r,c) = b.pH(dm_end);
        R.base_saturation(r,c) = sum(f(1:4));
        R.f_Ca(r,c) = f(1); R.f_Mg(r,c) = f(2); R.f_K(r,c) = f(3);
        R.f_Na(r,c) = f(4); R.f_Al(r,c) = f(5); R.f_H(r,c) = f(6);
        R.CEC_umolc_m2(r,c) = b.CEC_tot(dm_end);
        R.Ca_umol_L(r,c) = b.Ca(dm_end);
        R.Mg_umol_L(r,c) = b.Mg(dm_end);
        R.K_umol_L(r,c) = b.K(dm_end);
        R.Na_umol_L(r,c) = b.Na(dm_end);
        R.Al3_umol_L(r,c) = b.Al(dm_end)/scale;
        R.Al_dissolved_umol_L(r,c) = b.Al_w(dm_end)/scale;
        R.An_umolc_L(r,c) = b.An(dm_end);
        R.DIC_umol_L(r,c) = b.DIC(dm_end);
        p = PAR{c};
        if isempty(hour_rows)
            R.theta_bio(r,c) = sum(p.Bio_Zs.*O(1,:,c));
        else
            R.theta_bio(r,c) = sum(p.Bio_Zs.*mean(O(hour_rows,:,c),1));
        end
        if ~isempty(dm_rows)
            dtw = ERW{c}.dtw;
            R.leach_Na_umol_m2(r,c) = sum(b.Leak_Na(dm_rows))*dtw;
            R.leach_Al_umol_m2(r,c) = sum(b.Leak_Al(dm_rows))*dtw/scale;
            R.leach_An_umolc_m2(r,c) = sum(b.Leak_An(dm_rows))*dtw;
            R.leach_DIC_umol_m2(r,c) = sum(b.Leak_DIC(dm_rows))*dtw;
            R.drainage_mm(r,c) = sum(Lk(hour_rows,1,c));
        end
        if ~isempty(t_row) && t_row >= 2
            R.leach_Ca_g_m2(r,c) = LEAK_Ca(t_row,1,c);
            R.leach_Mg_g_m2(r,c) = LEAK_Mg(t_row,1,c);
            R.leach_K_g_m2(r,c) = LEAK_K(t_row,1,c);
            R.leach_NO3_gN_m2(r,c) = LEAK_NO3(t_row,1,c);
            R.leach_NH4_gN_m2(r,c) = LEAK_NH4(t_row,1,c);
        end
    end
end
R.state_names = state;
R.flux_names = flux;
end

%% =======================================================================
function result = assemble(segments,topo,CellInputs)
first = segments{1};
names = [first.state_names first.flux_names {'time','hours','precipitation_mm'}];
result = struct();
for k = 1:numel(names)
    parts = cellfun(@(R) R.(names{k}),segments,'UniformOutput',false);
    result.(names{k}) = vertcat(parts{:});
end
if any(diff(result.time) <= 0) || abs(result.time(1) - floor(result.time(1)+1e-9)) > 1e-6
    error('Run_Field_Transient:RecordOrder','The segments'' records do not follow one another.');
end
result.state_names = first.state_names;
result.flux_names = first.flux_names;
result.cell_id = topo.cell_id(:).';
[area,w] = Field_Aggregation_Areas(topo);
result.cell_overlap_area_m2 = area(:).';
result.cell_area_m2 = topo.cell_area_m2(:).';
result.field_weights = w;
A = vertcat(CellInputs.agronomic_input);
cells = table(topo.cell_id(:),'VariableNames',{'cell_id'});
for k = {'cell_x_projected_m','cell_y_projected_m','bulk_pH', ...
        'sand_fraction','clay_fraction','OM_fraction','cover_tree','cover_grass'}
    hit = find(strcmp(A.Properties.VariableNames,k{1}),1);
    if ~isempty(hit); cells.(k{1}) = A.(k{1}); end
end
%% the CEC each cell ran with (a scenario's, or its own)
cells.CEC_cmolc_kg = arrayfun(@(ci) double(ci.ERW_PARAM.CEC_tot_in),CellInputs(:));
if isfield(topo,'topography_attached') && topo.topography_attached
    cells.z_m = topo.cell_z_m(:); cells.slope = topo.cell_slope(:);
    cells.aspect_deg = topo.cell_aspect_deg(:); cells.twi = topo.cell_twi(:);
    cells.is_outlet = topo.cell_is_outlet(:);
end
result.cells = cells;
%% Field means over the cells' in-field areas; pH as the mean hydrogen
%% activity (pH_Area_Weighted).
F = struct('time',result.time,'hours',result.hours,'precipitation_mm',result.precipitation_mm);
for k = [first.state_names first.flux_names]
    x = result.(k{1});
    if any(strcmp(k{1},{'pH_bulk','pH_pore'}))
        F.(k{1}) = pH_Area_Weighted(x,w);
    else
        F.(k{1}) = x*w(:);
    end
end
result.field = F;
result.units = struct('time','datenum, local standard time (end of the record)', ...
    'hours','hours the record integrates', ...
    'pH_bulk','bulk soil pH (the site''s extraction ratio)','pH_pore','pore-water pH', ...
    'base_saturation','(Ca+Mg+K+Na)/CEC, charge fraction', ...
    'f_X','charge fraction of the CEC held by X','CEC_umolc_m2','umolc m-2 of the biogeochemical layer', ...
    'X_umol_L','umol per litre of pore water','Al3_umol_L','free Al3+', ...
    'Al_dissolved_umol_L','total dissolved Al','An_umolc_L','residual strong anion, umolc/L', ...
    'theta_bio','m3 m-3, the biogeochemical layer, mean over the record', ...
    'leach_X_g_m2','g m-2 over the record (soil biogeochemistry; NaN on the last record)', ...
    'leach_X_umol_m2','umol (Al, Na, DIC) or umolc (An) m-2 over the record (ERW)', ...
    'drainage_mm','mm out of the bottom of the soil column over the record', ...
    'precipitation_mm','mm over the record', ...
    'field','area-weighted field means; pH by mean hydrogen activity');
end

function write_field_mean_csv(result,csv_file)
F = result.field;
T = table(string(datestr(F.time,'yyyy-mm-dd HH:MM')),F.hours,F.precipitation_mm, ...
    'VariableNames',{'time','hours','precipitation_mm'});
for k = [result.state_names result.flux_names]
    T.(k{1}) = F.(k{1});
end
tmp = [tempname(fileparts(csv_file)) '.csv'];
writetable(T,tmp);
movefile(tmp,csv_file,'f');
end

%% =======================================================================
function T = subset_topology(topo,idx)
%% A smoke test's cells as a field of their own: the per-cell (cell_*) and
%% per-edge (edge_*) fields of the kept cells, the faces to the dropped
%% cells made part of the perimeter.
idx = idx(:);
N = topo.N_cells;
keep = false(N,1); keep(idx) = true;
new_index = zeros(N,1); new_index(idx) = (1:numel(idx)).';
inside = keep(topo.edge_a(:)) & keep(topo.edge_b(:));
T = topo;
for name = fieldnames(topo).'
    x = topo.(name{1});
    if startsWith(name{1},'cell_') && size(x,1) == N
        T.(name{1}) = rows_of(x,idx);
    elseif startsWith(name{1},'edge_') && size(x,1) == topo.N_edges
        T.(name{1}) = rows_of(x,find(inside));
    end
end
T.edge_a = reshape(new_index(topo.edge_a(inside)),size(T.edge_a));
T.edge_b = reshape(new_index(topo.edge_b(inside)),size(T.edge_b));
T.N_edges = nnz(inside);
T.N_cells = numel(idx);
cut = xor(keep(topo.edge_a(:)),keep(topo.edge_b(:)));
kept_end = topo.edge_a(:);
kept_end(cut & ~keep(topo.edge_a(:))) = topo.edge_b(cut & ~keep(topo.edge_a(:)));
face_width = topo.edge_face_width_m(:);
added = accumarray(kept_end(cut),face_width(cut),[N 1]);
T.outer_face_width_m = reshape(topo.outer_face_width_m(idx),[],1) + added(idx);
T.adj = topo.adj(idx,idx);
T.simulated_area_m2 = sum(T.cell_area_m2);
if isfield(T,'cell_overlap_area_m2')
    T.field_area_m2 = sum(T.cell_overlap_area_m2);
end
if isfield(T,'field_missed_area_m2'); T.field_missed_area_m2 = 0; end
T.subset_of_field_cells = N;
end

function y = rows_of(x,rows)
sz = size(x);
y = reshape(x(rows,:),[numel(rows) sz(2:end)]);
end

%% =======================================================================
function txt = cell_identity(CellInputs)
parts = arrayfun(@(ci) sprintf('%d:%s:%s:%s',ci.cell_id,char(string(ci.agronomic_identity)), ...
    ERW_Runtime_Parameter_Identity(ci.MOD_PARAM),ERW_Runtime_Parameter_Identity(ci.ERW_PARAM)), ...
    CellInputs(:),'UniformOutput',false);
txt = sha256_hex(strjoin(parts,';'));
end

function txt = terrain_identity(topo)
edges = [topo.edge_a(:) topo.edge_b(:)];
if isfield(topo,'topography_attached') && topo.topography_attached
    txt = Field_Transient_Carry('digest',{edges,topo.cell_z_m,topo.cell_slope, ...
        topo.cell_aspect_deg,topo.cell_svf,topo.cell_twi,topo.cell_is_outlet, ...
        topo.cell_area_m2,topo.cell_overlap_area_m2});
else
    txt = ['flat:' Field_Transient_Carry('digest',{edges,topo.cell_area_m2})];
end
end

function fp = code_fingerprint(v2)
%% The files that decide the physics, by content (the driver's checkpoint
%% covers the same), with the repository's location left out.
code = fullfile(v2,'T&C_Code');
txt = ERW_Checkpoint_Fingerprint(0,0,'field_transient',fullfile(code,'*.m'), ...
    fullfile(code,'Parameters','*.m'),fullfile(v2,'run','Run_Field_Transient.m'), ...
    fullfile(v2,'spinup','Run_Field_Cell_Spinup.m'));
txt = strrep(strrep(txt,char(92),'/'),strrep(v2,char(92),'/'),'<V2>');
fp = sha256_hex(txt);
end

%% =======================================================================
function L = longest_day(Date,Lat,Lon,DeltaGMT,t_bef,t_aft)
%% The longest day over the days of Date (from 00:00), found as the
%% parameter files find it (MOD_PARAM_*: Lmax_day = max(L_day), the day
%% length at 01:00 of each day). It scales leaf photosynthesis with day
%% length (RELATIVE_PC), so a segment takes the whole window's.
[YE,MO,DA,HO] = datevec(Date(:));
Datam = [YE MO DA HO];
L_day = zeros(numel(Date),1);
for j = 2:24:numel(Date)
    [~,~,~,~,~,L_day(j)] = SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
L = max(L_day);
end

function met = slice_met(met,mask)
for name = fieldnames(met).'
    met.(name{1}) = met.(name{1})(mask);
end
end

function path_out = meteo_file_path(file_name,Data_Directory,Input_Directory)
%% As the spin-up resolves its meteo file: an absolute path, else the Data
%% directory, else Inputs.
candidates = {file_name,fullfile(Data_Directory,file_name),fullfile(Input_Directory,file_name)};
hit = find(cellfun(@isfile,candidates),1);
if isempty(hit)
    error('Run_Field_Transient:MeteoFileNotFound','Meteo file %s not found (tried %s and %s).', ...
        file_name,candidates{2},candidates{3});
end
path_out = candidates{hit};
end

function v = row_value(T,name)
v = [];
if ismember(name,T.Properties.VariableNames)
    v = T.(name);
    if iscell(v); v = v{1}; end
end
end

function x = row_number(T,name)
v = row_value(T,name);
x = NaN;
if isnumeric(v) && isscalar(v); x = double(v); end
if (isstring(v) || ischar(v)); x = str2double(v); end
end

function dn = row_datenum(T,name)
v = row_value(T,name);
if isempty(v)
    dn = NaN;
elseif isdatetime(v)
    dn = datenum(v);
elseif isnumeric(v)
    dn = double(v);
else
    dn = datenum(datetime(char(string(v)),'InputFormat','yyyy-MM-dd'));
end
end

function cols = column_struct(v,name)
%% A struct of row columns, or a one-row CSV of them (Read_Parameter_Columns).
if isempty(v)
    cols = [];
elseif isstruct(v)
    cols = v;
elseif ischar(v) || (isstring(v) && isscalar(v))
    cols = Read_Parameter_Columns(char(v));
else
    error('Run_Field_Transient:BadOptions','%s must be a struct or a CSV file.',name);
end
if isstruct(cols) && isempty(fieldnames(cols)); cols = []; end
end

function r = dem_resolution(requested,topo)
r = requested;
if isempty(r); r = topo.cell_size_m; end
end

function f = existing_file(f,name)
if ~(ischar(f) || (isstring(f) && isscalar(f))) || ~isfile(f)
    error('Run_Field_Transient:MissingInput','%s not found: %s',name,char(string(f)));
end
f = canonical(char(f));
end

function p = canonical(p)
p = char(java.io.File(char(p)).getCanonicalPath());
end

function positive_integer(v,name)
if ~(isnumeric(v) && isscalar(v) && isfinite(v) && v == fix(v) && v >= 1)
    error('Run_Field_Transient:BadOptions','%s must be a positive integer.',name);
end
end

function x = env_number(name,default)
x = str2double(strtrim(getenv(name)));
if ~isfinite(x); x = default; end
end

function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end

function save_atomic(file_name,S)
tmp = [tempname(fileparts(file_name)) '.mat'];
cleanup = onCleanup(@() delete_if_present(tmp));
save(tmp,'-struct','S','-v7.3');
movefile(tmp,file_name,'f');
clear cleanup
end

function delete_if_present(f)
if isfile(f); delete(f); end
end

function started = start_thread_pool(workers)
started = false;
if workers < 2; return; end
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    error('Run_Field_Transient:MissingParallelToolbox', ...
        'ParallelWorkers needs the Parallel Computing Toolbox.');
end
pool = gcp('nocreate');
if isempty(pool)
    parpool('Threads',workers);
    started = true;
elseif ~isa(pool,'parallel.ThreadPool') || pool.NumWorkers < workers
    error('Run_Field_Transient:WrongParallelPool', ...
        'ParallelWorkers needs a thread pool of at least %d workers.',workers);
end
fprintf('Hourly column steps and daily ERW chemistry: %d thread workers\n',workers);
end

function stop_thread_pool(started)
if started
    pool = gcp('nocreate');
    if ~isempty(pool); delete(pool); end
end
end

function hex = sha256_hex(txt)
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(unicode2native(char(txt),'UTF-8'),'int8'));
raw = typecast(md.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end
