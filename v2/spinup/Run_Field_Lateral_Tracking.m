function out = Run_Field_Lateral_Tracking(topo,CellInputs,cell_ids,seeds, ...
    track_idx,output_directory,case_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Field_Lateral_Tracking
%%
%% The 2-D stage of the flux-tracking field spin-up (Lian et al. 2025,
%% Fig. 3c).  Every cell starts from its class seed and the whole field
%% runs in 2-D (MAIN_FRAME_BG_ERW_2D, topography attached) with the ERW
%% chemistry off, as in every spin-up phase, over the spin-up's
%% meteorological ring and its handoff extension.  Meanwhile the lateral
%% water and solute exchange of the tracked cells is recorded
%% (Lateral_Track_Record) for their plot-scale replays (Run_BG_Spinup with
%% LATERAL_REPLAY_FILE).
%%
%% The run is split exactly as Run_BG_Spinup's Phase 3 splits its ring: one
%% segment per calendar cycle year, then the partial year up to the
%% transient start, every column restarting from the previous segment's end
%% state through the same hand-over (pools, soil moisture, nutrient
%% availability and uptake memory), with the same small-negative pool
%% projection.  A replayed cell therefore receives each year's exchange
%% under the very restart it makes itself.
%%
%% Each finished segment is written to <output>/track/segment_<label>.mat
%% under a key that chains the run's identity through every earlier
%% segment, and is reused on a rerun; the 2-D driver's rolling checkpoint
%% resumes inside a segment.  A cluster step that runs out of time is
%% therefore simply submitted again.
%%
%% INPUTS
%%   topo          topography-attached field topology (N cells)
%%   CellInputs    Configure_Cell_Agronomic_Inputs struct array (N)
%%   cell_ids      [N x 1] cell ids
%%   seeds         [N x 1] struct: each cell's starting state (its class
%%                 seed): B_out, O_out, RexmyI_out, NavlI_out, NupI_L_out
%%   track_idx     indices (into 1..N) of the tracked cells
%%   output_directory, case_number
%%
%% OUTPUT struct
%%   .record_files   {n x 1} record_<id>.mat: the cell's lateral exchange
%%                   over every segment, its terrain and the solar window
%%   .ringend_files  {n x 1} ringend_<id>.mat: the cell's state at the end
%%                   of the last cycle year, the seed of its replay
%%   .handoff_file   handoff_states.mat: every cell's state at the
%%                   transient start
%%   .key            the tracking run's provenance key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spinup_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(spinup_directory);
Directory = fullfile(v2_directory,'T&C_Code');
Input_Directory = fullfile(Directory,'Inputs');
Data_Directory = fullfile(v2_directory,'Data');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory,Input_Directory,Parameter_Directory,spinup_directory);

N = topo.N_cells;
track_idx = track_idx(:);
if numel(CellInputs) ~= N || numel(cell_ids) ~= N || numel(seeds) ~= N
    error('Run_Field_Lateral_Tracking:BadInput', ...
        'CellInputs, cell_ids and seeds must hold one entry per cell (%d).',N);
end
if ~(isfield(topo,'topography_attached') && topo.topography_attached)
    error('Run_Field_Lateral_Tracking:NoTopography', ...
        'The tracked field run needs the DEM attached to the topology.');
end
track_dir = fullfile(output_directory,'track');
if ~exist(track_dir,'dir'); mkdir(track_dir); end

%% ---- the ring and the extension, from the rows every cell shares ------
MOD_PARAM_input = CellInputs(1).MOD_PARAM;
ERW_PARAM_input = CellInputs(1).ERW_PARAM;
window_columns = {'meteo_file','meteo_start','meteo_end','transient_start','Lat','Lon'};
for c = 2:N
    for k = 1:numel(window_columns)
        a = window_value(CellInputs(1).MOD_PARAM,window_columns{k});
        b = window_value(CellInputs(c).MOD_PARAM,window_columns{k});
        if ~isequaln(a,b)
            error('Run_Field_Lateral_Tracking:WindowMismatch', ...
                'Cell %d has a different %s from cell %d; the field runs on one forcing.', ...
                cell_ids(c),window_columns{k},cell_ids(1));
        end
    end
end
Lat = MOD_PARAM_input.Lat; Lon = MOD_PARAM_input.Lon; DeltaGMT = round(Lon/15);
meteo_path = meteo_file_path(window_value(MOD_PARAM_input,'meteo_file'), ...
    Data_Directory,Input_Directory);
[spinup_all_met,~] = Load_Spinup_Meteo(meteo_path, ...
    window_datenum(MOD_PARAM_input,'meteo_start'), ...
    window_datenum(MOD_PARAM_input,'meteo_end'),Lat,Lon,DeltaGMT);
transient_start = window_datenum(MOD_PARAM_input,'transient_start');
[cycle_years,handoff_start_dn,handoff_end_dn] = Spinup_Cycle_Years( ...
    spinup_all_met.Date,transient_start);
[spinup_all_yr,~,~] = datevec(spinup_all_met.Date);
labels = arrayfun(@(y) sprintf('%d',y),cycle_years,'UniformOutput',false);
if isfinite(handoff_end_dn)
    labels{end+1} = 'handoff';
end
n_ring = numel(cycle_years);

%% ---- provenance: cells, inputs, seeds, tracked set, terrain, windows ---
id_parts = arrayfun(@(c) sprintf('%d:%s:%s:%s',cell_ids(c), ...
    char(string(CellInputs(c).agronomic_identity)), ...
    ERW_Runtime_Parameter_Identity(CellInputs(c).MOD_PARAM), ...
    ERW_Runtime_Parameter_Identity(CellInputs(c).ERW_PARAM)),(1:N).', ...
    'UniformOutput',false);
seed_parts = arrayfun(@(c) mat2str([seeds(c).B_out(:); seeds(c).O_out(:); ...
    seeds(c).RexmyI_out(:); seeds(c).NavlI_out(:); seeds(c).NupI_L_out(:)].',17), ...
    (1:N).','UniformOutput',false);
terrain = [topo.cell_z_m(:) topo.cell_slope(:) topo.cell_aspect_deg(:) topo.cell_svf(:)];
run_key = sha256_hex(sprintf('case=%d|cells=%s|seeds=%s|track=%s|terrain=%s|edges=%s|segments=%s|handoff=%.17g', ...
    case_number,strjoin(id_parts,';'),strjoin(seed_parts,';'), ...
    mat2str(cell_ids(track_idx).'),mat2str(terrain,17), ...
    mat2str([topo.edge_a(:) topo.edge_b(:)]),strjoin(labels,','),handoff_end_dn));

%% ---- the segments -----------------------------------------------------
state = seeds(:);
seg_key = run_key;
seg_records = cell(numel(labels),1);
ringend = [];
for s = 1:numel(labels)
    seg_key = sha256_hex([seg_key '|' labels{s}]);
    seg_file = fullfile(track_dir,sprintf('segment_%s.mat',labels{s}));
    reused = false;
    if isfile(seg_file)
        S = load(seg_file,'key');
        if isfield(S,'key') && strcmp(S.key,seg_key)
            S = load(seg_file,'end_states','records');
            reused = true;
            fprintf('  tracking segment %s: reusing %s\n',labels{s},seg_file);
        end
    end
    if ~reused
        if s <= n_ring
            mask = spinup_all_yr == cycle_years(s);
        else
            mask = spinup_all_met.Date >= handoff_start_dn & ...
                spinup_all_met.Date < handoff_end_dn;
        end
        fprintf('  tracking segment %s: %d columns, %d h\n',labels{s},N,nnz(mask));
        [S.end_states,S.records] = run_segment(topo,CellInputs,state,track_idx, ...
            spinup_all_met,mask,MOD_PARAM_input,ERW_PARAM_input,case_number, ...
            Lat,Lon,DeltaGMT,Directory,Input_Directory,Parameter_Directory, ...
            fullfile(track_dir,sprintf('checkpoint_%s.mat',labels{s})));
        key = seg_key; end_states = S.end_states; records = S.records; %#ok<NASGU>
        label = labels{s}; %#ok<NASGU>
        tmp = [tempname(track_dir) '.mat'];
        save(tmp,'key','end_states','records','label','-v7.3');
        movefile(tmp,seg_file,'f');
        ckpt = fullfile(track_dir,sprintf('checkpoint_%s.mat',labels{s}));
        if isfile(ckpt); delete(ckpt); end
    end
    state = S.end_states(:);
    seg_records{s} = S.records;
    if s == n_ring
        ringend = state;
    end
end

%% ---- per-cell records, replay seeds, handoff states -------------------
out = struct('key',run_key,'record_files',{cell(numel(track_idx),1)}, ...
    'ringend_files',{cell(numel(track_idx),1)}, ...
    'handoff_file',fullfile(track_dir,'handoff_states.mat'));
%% The 2-D driver's TOPO_SOLAR_WINDOW default (hour-ending forcing).
solar_window = [1 0];
for k = 1:numel(track_idx)
    c = track_idx(k);
    parts = repmat(seg_records{1}(k),numel(seg_records),1);
    for s = 2:numel(seg_records); parts(s) = seg_records{s}(k); end
    record = struct();
    record.cell_id = cell_ids(c);
    record.hour_dn = vertcat(parts.hour_dn);
    record.q_runon = vertcat(parts.q_runon);
    record.dtheta = vertcat(parts.dtheta);
    record.day_dn = vertcat(parts.day_dn);
    record.retain = vertcat(parts.retain);
    record.input = vertcat(parts.input);
    record.species = parts(1).species;
    record.pool_idx = parts(1).pool_idx;
    record.slope = topo.cell_slope(c);
    record.aspect_deg = topo.cell_aspect_deg(c);
    record.svf = topo.cell_svf(c);
    record.solar_window = solar_window;
    record.cycle_years = cycle_years;
    record.required = true;
    record.tracking_key = run_key;
    %% What a replay is a replay of: its provenance carries this digest.
    record.identity = sha256_arrays({record.hour_dn,record.q_runon,record.dtheta, ...
        record.day_dn,record.retain,record.input,record.pool_idx, ...
        [record.slope record.aspect_deg record.svf solar_window]});
    f = fullfile(track_dir,sprintf('record_%04d.mat',cell_ids(c)));
    tmp = [tempname(track_dir) '.mat'];
    save(tmp,'record','-v7.3');
    movefile(tmp,f,'f');
    out.record_files{k} = f;
    seed = ringend(c);
    seed.cell_id = cell_ids(c);
    seed.tracking_key = run_key;
    f = fullfile(track_dir,sprintf('ringend_%04d.mat',cell_ids(c)));
    tmp = [tempname(track_dir) '.mat'];
    save(tmp,'-struct','seed','-v7.3');
    movefile(tmp,f,'f');
    out.ringend_files{k} = f;
end
handoff_states = state; %#ok<NASGU>
handoff_dn = handoff_end_dn; %#ok<NASGU>
if ~isfinite(handoff_end_dn); handoff_dn = NaN; end %#ok<NASGU>
key = run_key; %#ok<NASGU>
tmp = [tempname(track_dir) '.mat'];
save(tmp,'handoff_states','handoff_dn','cell_ids','key','-v7.3');
movefile(tmp,out.handoff_file,'f');
end

%% =======================================================================
function [end_states,records] = run_segment(topo,CellInputs,state,track_idx, ...
    met,mask,MOD_PARAM_input,ERW_PARAM_input,NN_case_study_ERW_MOD, ...
    Lat,Lon,DeltaGMT,Directory,Input_Directory,Parameter_Directory,ERW_ckpt_file) %#ok<INUSD>
%% The workspace MAIN_FRAME_BG_ERW_2D reads, for one segment, built as
%% Run_BG_Spinup builds each Phase-3 year's: the same forcing slice, the
%% same derived series and the same solar window.
Date = met.Date(mask);  Pr = met.Pr(mask);  Ta = met.Ta(mask);
Ws   = met.Ws(mask);    ea = met.ea(mask);  Pre = met.Pre(mask);
SAD1 = met.SAD1(mask);  SAD2 = met.SAD2(mask);
SAB1 = met.SAB1(mask);  SAB2 = met.SAB2(mask);
N    = met.N(mask);     Tdew = met.Tdew(mask);  esat = met.esat(mask);
PARB = met.PARB(mask);  PARD = met.PARD(mask);
NN   = numel(Date);
NNd  = ceil(NN/24) + 1;
Ws(Ws <= 0) = 0.01;
Ds   = max(esat - ea, 0);
[YE,MO,DA,HO,~,~] = datevec(Date);
Datam = [YE(:) MO(:) DA(:) HO(:)];
if Datam(1,4) ~= 0
    error('Run_Field_Lateral_Tracking:SegmentStart', ...
        'A tracking segment must start at 00:00 (it starts %s).', ...
        datestr(Date(1),'yyyy-mm-dd HH:MM'));
end
co2_record_file = 'Ca_Data.mat';
if ismember('co2_file',MOD_PARAM_input.Properties.VariableNames) && ...
        strlength(strtrim(string(MOD_PARAM_input.co2_file))) > 0
    co2_record_file = char(strtrim(string(MOD_PARAM_input.co2_file)));
end
Ca = Load_CO2_Record(Date,Input_Directory,co2_record_file);
dt = 3600; dth = 1; Oa = 210000;
t_bef = -0.67; t_aft = 1.67;
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'));
ms = length(Zs) - 1;
cc = MOD_PARAM_input.cc;
%% Run_Automatic_Site_Spinup's id_location: the row's own, else Location.
id_location = char(string(MOD_PARAM_input.Location));
if ismember('id_location',MOD_PARAM_input.Properties.VariableNames)
    v = MOD_PARAM_input.id_location;
    if iscell(v); v = v{1}; end
    if ~isempty(v) && ~(isstring(v) && ismissing(v)) && strlength(strtrim(string(v))) > 0
        id_location = char(string(v));
    end
end
IrD = zeros(NN,1);
OPT_SoilBiogeochemistry = 1;
OPT_ERW = 0;                     %% ERW chemistry off, as in every spin-up phase
parallel_workers = 0;
%% The rolling checkpoint holds every column's history: hourly is enough to
%% lose little on a time limit without writing tens of GB every few minutes.
ERW_ckpt_minutes = 60;
LATERAL_TRACK = struct('cells',track_idx);
base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
ColCfg = repmat(struct('PARAM_IC',base_param,'PFT_OVERRIDE','', ...
    'apply_rock',false,'ERW_PARAM',ERW_PARAM_input,'MOD_PARAM',table(), ...
    'BG_STATE',struct(),'agronomic_identity',''),topo.N_cells,1);
for c = 1:topo.N_cells
    ColCfg(c).MOD_PARAM = CellInputs(c).MOD_PARAM;
    ColCfg(c).ERW_PARAM = CellInputs(c).ERW_PARAM;
    %% The vegetation override the cell's row names, as the 1-D driver
    %% applies it in the cell's own spin-ups.
    if ismember('pft_override',CellInputs(c).MOD_PARAM.Properties.VariableNames)
        ColCfg(c).PFT_OVERRIDE = char(strtrim(string(CellInputs(c).MOD_PARAM.pft_override)));
    end
    ColCfg(c).agronomic_identity = CellInputs(c).agronomic_identity;
    ColCfg(c).BG_STATE = struct('B_out',state(c).B_out,'O_out',state(c).O_out, ...
        'RexmyI_out',state(c).RexmyI_out,'NavlI_out',state(c).NavlI_out, ...
        'NupI_L_out',state(c).NupI_L_out);
end
clear YE MO DA HO met mask state c
%% The driver's day counter, declared here: a j the script alone creates
%% would be read below as the imaginary unit.
j = 1;
old = cd(Directory);
cleanup = onCleanup(@() cd(old));
Plant_Carbon_Deficit('reset');
MAIN_FRAME_BG_ERW_2D;
clear cleanup
%% Every column's end state, built and projected as Run_BG_Spinup hands a
%% Phase-3 year to the next.
end_states = repmat(struct('B_out',[],'O_out',[],'RexmyI_out',[],'NavlI_out',[], ...
    'NupI_L_out',[]),topo.N_cells,1);
for c = 1:topo.N_cells
    st = build_bg_state(P(:,:,c),O(:,:,c),RexmyI(:,:,c),NavlI(:,:,c),NupI_L(:,:,:,c), ...
        Nreserve_L(:,:,c),Preserve_L(:,:,c),Kreserve_L(:,:,c),Careserve_L(:,:,c), ...
        Mgreserve_L(:,:,c),Sireserve_L(:,:,c),Nreserve_H(:,:,c),Preserve_H(:,:,c), ...
        Kreserve_H(:,:,c),Careserve_H(:,:,c),Mgreserve_H(:,:,c),Sireserve_H(:,:,c), ...
        id_location);
    end_states(c).B_out = project_small_negative(st.B_out,1e-6);
    end_states(c).O_out = st.O_out;
    end_states(c).RexmyI_out = st.RexmyI_out;
    end_states(c).NavlI_out = st.NavlI_out;
    end_states(c).NupI_L_out = st.NupI_L_out;
end
%% The steps (rows 2:NN) and days (rows 2:NNd) of the tracked cells; row 1
%% is the segment's initial state and carries no exchange.
if any(lateral_track.retain(1,:,:) ~= 1,'all') || any(lateral_track.input(1,:,:) ~= 0,'all')
    error('Run_Field_Lateral_Tracking:ExchangeBeforeFirstDay', ...
        'Lateral exchange was recorded before the first daily block.');
end
records = repmat(struct('hour_dn',[],'q_runon',[],'dtheta',[],'day_dn',[], ...
    'retain',[],'input',[],'species',{{}},'pool_idx',[]),numel(track_idx),1);
days = 2:j;
for k = 1:numel(track_idx)
    rec = Lateral_Track_Record('cell',lateral_track,k);
    records(k).hour_dn = rec.hour_dn(2:NN);
    records(k).q_runon = rec.q_runon(2:NN);
    records(k).dtheta = rec.dtheta(2:NN,:);
    records(k).day_dn = rec.day_dn(days);
    records(k).retain = rec.retain(days,:);
    records(k).input = rec.input(days,:);
    records(k).species = rec.species;
    records(k).pool_idx = rec.pool_idx;
end
end

%% =======================================================================
function B = project_small_negative(B,tolerance)
%% Run_BG_Spinup's projection between Phase-3 years.
large_negative = find(B < -tolerance);
if ~isempty(large_negative)
    error('Run_Field_Lateral_Tracking:NegativeBiogeochemicalPool', ...
        'Biogeochemical pools %s are below the projection limit (minimum %.6g).', ...
        mat2str(large_negative(:).'),min(B(large_negative)));
end
B(B < 0) = 0;
end

function v = window_value(T,name)
v = [];
if ismember(name,T.Properties.VariableNames)
    v = T.(name);
    if iscell(v); v = v{1}; end
    if isstring(v); v = char(v); end
end
end

function dn = window_datenum(T,name)
v = window_value(T,name);
if isempty(v)
    dn = NaN;
elseif isdatetime(v)
    dn = datenum(v);
elseif isnumeric(v)
    dn = double(v);
else
    try
        dn = datenum(datetime(char(v),'InputFormat','yyyy-MM-dd'));
    catch
        dn = datenum(datetime(char(v)));
    end
end
end

function path_out = meteo_file_path(file_name,Data_Directory,Input_Directory)
%% Resolved as Run_Automatic_Site_Spinup resolves a spin-up meteo file:
%% an absolute path, else the Data directory, else Inputs.
file_name = char(file_name);
candidates = {file_name, fullfile(Data_Directory,file_name), ...
    fullfile(Input_Directory,file_name)};
hit = find(cellfun(@isfile,candidates),1);
if isempty(hit)
    error('Run_Field_Lateral_Tracking:MeteoFileNotFound', ...
        'Spin-up meteo file %s not found (tried %s and %s).', ...
        file_name,candidates{2},candidates{3});
end
path_out = candidates{hit};
end

function hex = sha256_arrays(parts)
%% Digest of the binary values of numeric arrays (exact, and far quicker
%% than printing a multi-year hourly record).
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
for k = 1:numel(parts)
    x = double(parts{k});
    digest.update(typecast(int64(size(x)),'int8'));
    digest.update(typecast(x(:),'int8'));
end
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end

function hex = sha256_hex(txt)
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(unicode2native(txt,'UTF-8'),'int8'));
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end
