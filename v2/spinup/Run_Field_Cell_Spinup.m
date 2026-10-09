function [state_file,CellInputs] = Run_Field_Cell_Spinup(shape_file,agronomic_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spin up the soil biogeochemistry of every cell in one field polygon.
%%
%% Method 'flux_tracking' (default): the T&C-BG-2D spin-up of Lian et al.
%% (2025, WRR 61, e2025WR040260, section 2.3 and Fig. 3), adapted for ERW.
%%   1 class seeds    one plot-scale spin-up per land cover x soil class
%%                    (Run_BG_Spinup: long-term mean forcing to steady
%%                    state, then the coupled ring to dynamic steady state)
%%                    on the class's most typical cell, assigned to every
%%                    cell of the class (Fig. 3a-b);
%%   2 tracking       the whole field in 2-D from those seeds over the
%%                    spin-up's meteorological ring and handoff extension,
%%                    recording the lateral water and solute exchange of a
%%                    topographically stratified ~10 % subsample
%%                    (Run_Field_Lateral_Tracking, Fig. 3c);
%%   3 replays        each subsample cell's own plot-scale spin-up with that
%%                    exchange replayed in every phase, to its dynamic
%%                    steady state (Run_BG_Spinup, LATERAL_REPLAY_FILE;
%%                    Fig. 3d);
%%   4 reconstruction a random forest trained on the replayed cells predicts
%%                    every other cell's conserved pool totals from its
%%                    terrain and soil (Reconstruct_Cell_Pools_RF; Fig. 3e).
%% Adapted for ERW: the ERW chemistry is off in every stage, as in every
%% spin-up phase of v0.9.  The exchanger and pore water are never spun up:
%% every cell, replayed or reconstructed, is seeded at the handoff from its
%% own soil chemistry (bulk pH, CEC, measured base saturation) through
%% v0.9's seed path and must pass its own ERW charge-balance handoff gate.
%% The forest predicts biological and nutrient pools only.
%%
%% Method 'per_cell': every cell's own plot-scale spin-up without lateral
%% exchange (the 1-D reference).
%%
%% Stages, one SLURM step or array each (v2/cluster):
%%   'Stage','plan'        spinup_plan.csv and the array-index lists
%%   'Stage','class_seed'  'ClassIndex',k: the seed of class k, with its
%%                         ERW seed's inputs beside it for an exchange
%%                         calibration (classseed_<k>_erw_seed_inputs.mat,
%%                         Fit_Exchange_Affinity)
%%   'Stage','track'       the 2-D tracking run (every class seed on disk)
%%   'Stage','replay'      'CellIndex',c: subsample cell c's replay
%%   'Stage','reconstruct' 'CellIndex',c: cell c from the forest (every
%%                         replay on disk)
%%   'Stage','cell'        per_cell: 'CellIndex',c, cell c's spin-up
%%   'Stage','assemble'    every cellstate_<id>.mat into the validated
%%                         CellBGStates.mat the transient run reads
%% Without a Stage every stage runs here in turn.  'PlanOnly',true and
%% 'AssembleOnly',true are the plan and assemble stages; 'CellIndex' alone
%% is the per_cell stage, or the replay or reconstruct stage of that cell.
%% Finished work (a class seed, a tracking segment, a replay, a
%% reconstruction) is reused under its provenance key, so an interrupted
%% stage or array resumes.
%%
%% 'WarmStartStates', file: start from an earlier field spin-up's
%% CellBGStates.mat on the same grid -- the first spin-up, once site
%% observations replace its provisional soil inputs.  Every cell starts
%% from its own converged state and runs Phase 3 only (Run_BG_Spinup,
%% SPINUP_WARM_START): per_cell, each cell's own spin-up; flux_tracking, no
%% class seeds -- the tracking starts from every cell's earlier state, the
%% replays run Phase 3 with their recorded exchange, and the forest
%% reconstructs the rest as before.  The exchanger and pore water are
%% seeded from the new inputs at the handoff, as always.
%% 'MODParamColumns' / 'ERWParamColumns', struct: columns over the case's
%% MOD and ERW rows for every cell (a site's own settings, the spin-up
%% windows); each cell's soil and land cover come from the agronomic file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'GridCellSizeM',[],'CaseNumber',5,'OutputFile','','DryRun',false, ...
    'Method','flux_tracking','FieldDEM','','SubsampleFraction',0.10, ...
    'MinCoverage',[],'MODParamColumns',[],'CellIndex',[],'ClassIndex',[], ...
    'Stage','','PlanOnly',false,'AssembleOnly',false,'RFOptions',{{}}, ...
    'ClassColumns',{{}},'WarmStartStates','','ERWParamColumns',[]);
if mod(numel(varargin),2) ~= 0
    error('Run_Field_Cell_Spinup:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_Field_Cell_Spinup:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~isfile(shape_file) || ~isfile(agronomic_file)
    error('Run_Field_Cell_Spinup:MissingInput', ...
        'The shapefile and agronomic input file must both exist.');
end
if ~(islogical(opt.DryRun) && isscalar(opt.DryRun))
    error('Run_Field_Cell_Spinup:BadOptions','DryRun must be scalar logical.');
end
if ~(isnumeric(opt.CaseNumber) && isscalar(opt.CaseNumber) && ...
        isfinite(opt.CaseNumber) && opt.CaseNumber == fix(opt.CaseNumber) && ...
        opt.CaseNumber >= 1)
    error('Run_Field_Cell_Spinup:BadOptions','CaseNumber must be a positive integer.');
end
opt.Method = lower(char(opt.Method));
if strcmp(opt.Method,'rf')
    error('Run_Field_Cell_Spinup:BadMethod', ...
        ['Method ''rf'' is superseded by ''flux_tracking'' (Lian et al. 2025: ' ...
         'class seeds, a 2-D tracking run, replays of the tracked cells and the ' ...
         'random forest for the rest).']);
end
if ~any(strcmp(opt.Method,{'flux_tracking','per_cell'}))
    error('Run_Field_Cell_Spinup:BadMethod', ...
        'Method must be ''flux_tracking'' or ''per_cell''.');
end

spinup_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(spinup_directory);
code_directory = fullfile(v2_directory,'T&C_Code');
input_directory = fullfile(code_directory,'Inputs');
addpath(code_directory,fullfile(code_directory,'Parameters'),spinup_directory);
setups = Build_Shapefile_Grid_Setups(shape_file,opt.GridCellSizeM, ...
    opt.FieldFeatureIndices,opt.MinCoverage);
if isempty(opt.FieldPolygonIndex)
    if numel(setups) ~= 1
        error('Run_Field_Cell_Spinup:DisconnectedField', ...
            'Select one of the %d disconnected polygons with FieldPolygonIndex.', ...
            numel(setups));
    end
    opt.FieldPolygonIndex = 1;
elseif ~(isnumeric(opt.FieldPolygonIndex) && ...
        isscalar(opt.FieldPolygonIndex) && isfinite(opt.FieldPolygonIndex) && ...
        opt.FieldPolygonIndex == fix(opt.FieldPolygonIndex) && ...
        opt.FieldPolygonIndex >= 1 && opt.FieldPolygonIndex <= numel(setups))
    error('Run_Field_Cell_Spinup:BadPolygonIndex', ...
        'FieldPolygonIndex must be between 1 and %d.',numel(setups));
end
topo = setups(opt.FieldPolygonIndex);
resolved = Resolve_Cell_Agronomic_Inputs(topo, ...
    readtable(agronomic_file,'VariableNamingRule','preserve'));

workbook = fullfile(input_directory,'ERW_MOD_PARAM_parameters.xlsx');
[MOD, ERW] = Read_Parameter_Workbook(workbook);   % the spin-up's own typing
MOD = MOD(MOD.Number == opt.CaseNumber,:);
ERW = ERW(ERW.Number == opt.CaseNumber,:);
if height(MOD) ~= 1 || height(ERW) ~= 1
    error('Run_Field_Cell_Spinup:CaseSelection', ...
        'CaseNumber %d must select one MOD and one ERW row.',opt.CaseNumber);
end
%% The soil inputs and land cover of every cell come from the agronomic
%% file; a column override may not also set them.
per_cell_MOD = {'PH','Psan','Pcla','Porg'};
if ismember('cover_tree',resolved.Properties.VariableNames) && all(isfinite(resolved.cover_tree))
    per_cell_MOD = [per_cell_MOD {'Ccrown','Cbare','Curb','Cwat','Crock'}];
end
if ~isempty(opt.ERWParamColumns)
    %% Optional ERW row columns for every cell -- a site's own ERW settings
    %% over the case row it borrows -- applied before the cells' own soil
    %% inputs, so a cell's measured or baseline exchanger still wins.
    names = column_names(opt.ERWParamColumns,'ERWParamColumns');
    check_columns(names,ERW.Properties.VariableNames,{'CEC_tot_in'},'ERWParamColumns');
    for k = 1:numel(names)
        ERW.(names{k}) = column_value(opt.ERWParamColumns.(names{k}));
    end
end
CellInputs = Configure_Cell_Agronomic_Inputs(MOD,ERW,resolved);
if ~isempty(opt.MODParamColumns)
    %% Optional MOD_PARAM columns the workbook row does not carry (the
    %% spin-up meteo and transient windows, a forcing file), applied to
    %% every cell's row; the spin-up reads them through table_value.
    names = column_names(opt.MODParamColumns,'MODParamColumns');
    known = [MOD.Properties.VariableNames, {'id_location','forcing_file', ...
        'transient_start','transient_end','scenario_transient','source_data_file', ...
        'meteo_file','rain_chem_file','meteo_start','meteo_end','rain_chem_start', ...
        'rain_chem_end','pft_override','co2_file'}];
    check_columns(names,known,per_cell_MOD,'MODParamColumns');
    for c = 1:numel(CellInputs)
        for k = 1:numel(names)
            CellInputs(c).MOD_PARAM.(names{k}) = column_value(opt.MODParamColumns.(names{k}));
        end
    end
end
%% A row that seeds the exchanger from base saturation (v0.9's
%% f_CEC_inventory convention) needs it in every cell: without it the seed
%% fails at the handoff, hours into the run, so it is refused here.
no_bs = arrayfun(@(ci) lacks_exchanger_seed(ci.ERW_PARAM),CellInputs);
if any(no_bs)
    error('Run_Field_Cell_Spinup:NoBaseSaturation', ...
        ['The ERW row seeds the exchanger from base saturation (f_CEC_inventory) ' ...
         'but %d cell(s) have none (cell id(s) %s): give baseline_bs_* or ' ...
         'measured_bs_* columns in the agronomic file.'], ...
        nnz(no_bs),mat2str(resolved.cell_id(no_bs).'));
end

state_file = char(opt.OutputFile);
if isempty(state_file)
    [~,shape_name] = fileparts(shape_file);
    output_directory = fullfile(input_directory,'Field_Spinup',shape_name, ...
        sprintf('polygon_%02d',opt.FieldPolygonIndex));
    state_file = fullfile(output_directory,'CellBGStates.mat');
else
    output_directory = fileparts(state_file);
    if isempty(output_directory)
        output_directory = pwd;
        state_file = fullfile(output_directory,state_file);
    end
end
if opt.DryRun
    return
end
if ~exist(output_directory,'dir'); mkdir(output_directory); end

stage = resolve_stage(opt);
cluster_step = ~isempty(stage);
N = topo.N_cells;

%% ---- the flux-tracking design: terrain, subsample, classes -----------
track_idx = []; class_of = ones(N,1); class_rep = zeros(0,1);   % per_cell has no classes
if strcmp(opt.Method,'flux_tracking')
    if isempty(opt.FieldDEM)
        error('Run_Field_Cell_Spinup:FluxTrackingNeedsDEM', ...
            ['Method ''flux_tracking'' needs a FieldDEM: the lateral exchange it ' ...
             'tracks and the forest''s covariates come from the terrain. Use ' ...
             'Method ''per_cell'' for a spin-up without lateral exchange.']);
    end
    if any(strcmp(stage,{'','reconstruct'})) && exist('TreeBagger','file') ~= 2
        error('Run_Field_Cell_Spinup:NoRFToolbox', ...
            ['The flux-tracking reconstruction needs TreeBagger (Statistics and ' ...
             'Machine Learning Toolbox) on every node that runs it.']);
    end
    dem = Load_Field_DEM(char(opt.FieldDEM),topo, ...
        'OutResolutionM',dem_resolution(opt.GridCellSizeM,topo));
    topo = Attach_Topography_To_Topology(topo,dem);
    %% Decided once and deterministically (same seed, same grid), so the
    %% plan, every array task and the assembly agree on the design.
    [track_idx,~] = Select_Topographic_Subsample(topo,resolved,opt.SubsampleFraction);
    track_idx = track_idx(:);
    kfold = 5;
    for i = 1:2:numel(opt.RFOptions)
        if strcmpi(char(opt.RFOptions{i}),'KFold'); kfold = opt.RFOptions{i+1}; end
    end
    if numel(track_idx) < 2*kfold
        error('Run_Field_Cell_Spinup:SubsampleTooSmall', ...
            ['SubsampleFraction %.2f gives %d tracked cells on %d; the forest ' ...
             'needs at least 2 x KFold = %d. Raise SubsampleFraction or pass ' ...
             'RFOptions {KFold,k}.'],opt.SubsampleFraction,numel(track_idx),N,2*kfold);
    end
    [class_of,class_rep,class_names] = field_classes(resolved,opt.ClassColumns);
end
if cluster_step
    check_plan_method(output_directory,opt.Method);
end
F = struct('dir',output_directory,'case',opt.CaseNumber,'resolved',resolved, ...
    'inputs',CellInputs,'track_idx',track_idx,'class_of',class_of, ...
    'class_rep',class_rep,'warm',warm_seeds(opt.WarmStartStates,resolved,output_directory));
if strcmp(stage,'cell_auto')
    %% A cell index alone: a tracked cell's replay, any other's reconstruction.
    c = check_index(opt.CellIndex,N,'CellIndex');
    if ismember(c,track_idx); stage = 'replay'; else; stage = 'reconstruct'; end
    if strcmp(stage,'reconstruct') && exist('TreeBagger','file') ~= 2
        error('Run_Field_Cell_Spinup:NoRFToolbox', ...
            ['The flux-tracking reconstruction needs TreeBagger (Statistics and ' ...
             'Machine Learning Toolbox) on every node that runs it.']);
    end
end

switch stage
    case 'plan'
        state_file = write_plan(output_directory,topo,resolved,opt.Method,track_idx, ...
            class_of,class_rep);
        return
    case 'assemble'
        [CellBGStates,report] = assemble_cell_states(output_directory,topo, ...
            resolved,MOD,ERW,opt.Method);
        resolved.spinup_method = repmat(string(opt.Method),N,1);
        save_state_file(state_file,CellBGStates,resolved,report);
        stop_on_failed_cells(report,resolved,state_file);
        fprintf('Saved %d assembled cell spin-up states (%s): %s\n',N,opt.Method,state_file);
        return
    case 'class_seed'
        k = check_index(opt.ClassIndex,numel(class_rep),'ClassIndex');
        state_file = class_seed_path(output_directory,k);
        if F.warm.on
            fprintf('Class %d: no class seed in a warm start (the tracking starts from the earlier states)\n',k);
            state_file = '';
            return
        end
        class_seed(F,k);
        fprintf('Class %d seed (cell id %d): %s\n',k,resolved.cell_id(class_rep(k)),state_file);
        return
    case 'track'
        out = track_field(F,topo);
        state_file = out.handoff_file;
        return
    case 'replay'
        c = check_index(opt.CellIndex,N,'CellIndex');
        if ~ismember(c,track_idx)
            error('Run_Field_Cell_Spinup:NotTracked', ...
                'Cell %d (id %d) is not in the tracked subsample.',c,resolved.cell_id(c));
        end
        st = replay_cell(F,c);
        state_file = write_cell_state(output_directory,c,resolved.cell_id(c),st,"replay");
        fprintf('Saved cell %d (id %d, replay): %s\n',c,resolved.cell_id(c),state_file);
        return
    case 'reconstruct'
        c = check_index(opt.CellIndex,N,'CellIndex');
        [st,method] = reconstruct_cells(F,topo,opt,c);
        state_file = write_cell_state(output_directory,c,resolved.cell_id(c),st,method);
        fprintf('Saved cell %d (id %d, %s): %s\n',c,resolved.cell_id(c),method,state_file);
        if method == "failed"
            error('Run_Field_Cell_Spinup:CellFailed', ...
                ['Cell %d (id %d): the reconstruction and the cold start both ' ...
                 'failed; see the warnings above.'],c,resolved.cell_id(c));
        end
        return
    case 'cell'
        c = check_index(opt.CellIndex,N,'CellIndex');
        [label,extra] = cell_start(F,c);
        st = stamp(full_cell_spinup(spinup_state_path(output_directory, ...
            resolved.cell_id(c)),resolved.cell_id(c),opt.CaseNumber,CellInputs(c), ...
            label,extra),"full",false);
        state_file = write_cell_state(output_directory,c,resolved.cell_id(c),st,"full");
        fprintf('Saved cell %d (id %d, full): %s\n',c,resolved.cell_id(c),state_file);
        return
end

%% ---- every stage here, in turn ----------------------------------------
if strcmp(opt.Method,'flux_tracking')
    for k = 1:numel(class_rep)*~F.warm.on
        class_seed(F,k);
    end
    track_field(F,topo);
    for c = track_idx(:).'
        replay_cell(F,c);
    end
    [CellBGStates,report] = reconstruct_cells(F,topo,opt,[]);
    report.class_names = class_names;
else
    state_cells = cell(N,1);
    for c = 1:N
        [label,extra] = cell_start(F,c);
        state_cells{c} = stamp(full_cell_spinup(spinup_state_path(output_directory, ...
            resolved.cell_id(c)),resolved.cell_id(c),opt.CaseNumber,CellInputs(c), ...
            label,extra),"full",false);
    end
    CellBGStates = harmonise_states(state_cells);
    report = struct('method','per_cell','method_per_cell',repmat("full",N,1), ...
        'n_full',N,'n_failed',0);
end
resolved.spinup_method = repmat(string(opt.Method),N,1);
save_state_file(state_file,CellBGStates,resolved,report);
stop_on_failed_cells(report,resolved,state_file);
fprintf('Saved %d validated cell spin-up states (%s): %s\n',N,opt.Method,state_file);
end

%% =======================================================================
function stage = resolve_stage(opt)
stage = lower(char(opt.Stage));
legacy = {};
if opt.PlanOnly; legacy{end+1} = 'plan'; end
if opt.AssembleOnly; legacy{end+1} = 'assemble'; end
if numel(legacy) > 1 || (~isempty(legacy) && ~isempty(stage) && ~strcmp(stage,legacy{1}))
    error('Run_Field_Cell_Spinup:BadStage','Give one stage only.');
end
if ~isempty(legacy); stage = legacy{1}; end
if isempty(stage) && ~isempty(opt.CellIndex)
    %% A cell index alone names that cell's own stage.
    if strcmp(opt.Method,'per_cell'); stage = 'cell'; else; stage = 'cell_auto'; end
end
known = {'','plan','class_seed','track','replay','reconstruct','cell','assemble','cell_auto'};
if ~any(strcmp(stage,known))
    error('Run_Field_Cell_Spinup:BadStage', ...
        'Stage must be plan, class_seed, track, replay, reconstruct, cell or assemble.');
end
if strcmp(opt.Method,'per_cell') && any(strcmp(stage,{'class_seed','track','replay','reconstruct'}))
    error('Run_Field_Cell_Spinup:BadStage','Stage %s belongs to Method flux_tracking.',stage);
end
if strcmp(opt.Method,'flux_tracking') && strcmp(stage,'cell')
    error('Run_Field_Cell_Spinup:BadStage','Stage cell belongs to Method per_cell.');
end
end

function k = check_index(k,n,name)
if ~(isnumeric(k) && isscalar(k) && isfinite(k) && k == fix(k) && k >= 1 && k <= n)
    error(['Run_Field_Cell_Spinup:Bad' name],'%s must be an integer in 1..%d.',name,n);
end
end

%% =======================================================================
function [class_of,class_rep,class_names] = field_classes(resolved,columns)
%% Land cover x soil classes (Lian et al. 2025: one 1-D seed per class).
%% Categorical columns of the agronomic table -- the named ones, or by
%% default whichever of land_cover, vegetation, pft, soil_class and
%% texture_class it carries -- define the classes; a field with none is one
%% class.  Each class is represented by its most typical cell: the one
%% nearest the class median of the resolved soil properties.
N = height(resolved);
if isempty(columns)
    columns = {'land_cover','vegetation','pft','soil_class','texture_class'};
    columns = columns(ismember(lower(columns),lower(resolved.Properties.VariableNames)));
end
codes = zeros(N,0); parts = cell(N,numel(columns));
for k = 1:numel(columns)
    hit = find(strcmpi(resolved.Properties.VariableNames,columns{k}),1);
    if isempty(hit)
        error('Run_Field_Cell_Spinup:UnknownClassColumn', ...
            'The agronomic table has no column %s.',columns{k});
    end
    v = string(resolved.(resolved.Properties.VariableNames{hit}));
    [~,~,cc] = unique(v); codes = [codes cc(:)]; %#ok<AGROW>
    parts(:,k) = cellstr(v);
end
if isempty(codes)
    class_of = ones(N,1); class_names = {'field'};
else
    [u,~,class_of] = unique(codes,'rows');
    class_names = cell(size(u,1),1);
    for q = 1:size(u,1)
        first = find(class_of == q,1);
        class_names{q} = strjoin(parts(first,:),' / ');
    end
end
soil = {'sand_fraction','clay_fraction','OM_fraction','bulk_pH','CEC_cmolc_kg'};
X = zeros(N,0);
for k = 1:numel(soil)
    hit = find(strcmpi(resolved.Properties.VariableNames,soil{k}),1);
    if ~isempty(hit); X = [X double(resolved.(resolved.Properties.VariableNames{hit}))]; end %#ok<AGROW>
end
n_class = max(class_of);
class_rep = zeros(n_class,1);
for q = 1:n_class
    members = find(class_of == q);
    if isempty(X)
        class_rep(q) = members(1);
        continue
    end
    Xq = X(members,:);
    sd = std(X,0,1); sd(~(sd > 0)) = 1;
    d = sum(((Xq - median(Xq,1))./sd).^2,2);
    [~,best] = min(d);
    class_rep(q) = members(best);
end
end

%% =======================================================================
function st = class_seed(F,k)
%% The class's 1-D spin-up (Lian et al. Fig. 3a): its representative cell's
%% own full spin-up, written to the class seed file.  Its handoff saves the
%% ERW seed's inputs beside it, the input of an exchange calibration.
c = F.class_rep(k);
f = class_seed_path(F.dir,k);
st = full_cell_spinup(f,F.resolved.cell_id(c),F.case,F.inputs(c),'', ...
    struct('ERW_SEED_INPUTS_FILE',strrep(f,'.mat','_erw_seed_inputs.mat')));
end

function v = stand_of(st)
%% A state's perennial stand (Veg_out), or [] when it has none.
v = [];
if isfield(st,'Veg_out') && isstruct(st.Veg_out)
    v = st.Veg_out;
end
end

function f = class_seed_path(output_directory,k)
f = fullfile(output_directory,sprintf('classseed_%02d.mat',k));
end

function out = track_field(F,topo)
%% Lian et al. Fig. 3c: the field in 2-D from its class seeds -- or, in a
%% warm start, from every cell's own earlier state.
if F.warm.on
    seeds = arrayfun(@(s) struct('B_out',s.B_out(end,:),'O_out',s.O_out(end,:), ...
        'RexmyI_out',s.RexmyI_out(end,:),'NavlI_out',s.NavlI_out(end,:), ...
        'NupI_L_out',s.NupI_L_out,'Veg_out',{stand_of(s)}),F.warm.states(:),'UniformOutput',false);
    out = Run_Field_Lateral_Tracking(topo,F.inputs,F.resolved.cell_id,vertcat(seeds{:}), ...
        F.track_idx,F.dir,F.case);
    return
end
n_class = numel(F.class_rep);
seeds_by_class = cell(n_class,1);
for k = 1:n_class
    f = class_seed_path(F.dir,k);
    c = F.class_rep(k);
    if ~isfile(f)
        error('Run_Field_Cell_Spinup:ClassSeedNotReady', ...
            'Class %d seed %s is missing; run the class_seed stage first.',k,f);
    end
    st = load_bg_state(f);
    Validate_BG_State_For_Transient(st,f);
    want = Spinup_Provenance(F.resolved.cell_id(c),F.inputs(c).agronomic_identity, ...
        F.case,F.inputs(c).MOD_PARAM,F.inputs(c).ERW_PARAM,'cold');
    if ~(isfield(st,'spinup_provenance') && strcmp(char(string(st.spinup_provenance)),want))
        error('Run_Field_Cell_Spinup:ClassSeedStale', ...
            '%s was spun up for another run; rerun the class_seed stage.',f);
    end
    seeds_by_class{k} = struct('B_out',st.B_out(end,:),'O_out',st.O_out(end,:), ...
        'RexmyI_out',st.RexmyI_out(end,:),'NavlI_out',st.NavlI_out(end,:), ...
        'NupI_L_out',st.NupI_L_out,'Veg_out',{stand_of(st)});
end
seeds = vertcat(seeds_by_class{F.class_of});
out = Run_Field_Lateral_Tracking(topo,F.inputs,F.resolved.cell_id,seeds, ...
    F.track_idx,F.dir,F.case);
end

function st = replay_cell(F,c)
%% Lian et al. Fig. 3d: the tracked cell's plot-scale spin-up with its
%% recorded lateral exchange replayed, seeded from its own state at the end
%% of the tracking ring.
cell_id = F.resolved.cell_id(c);
track_dir = fullfile(F.dir,'track');
record_file = fullfile(track_dir,sprintf('record_%04d.mat',cell_id));
seed_file = fullfile(track_dir,sprintf('ringend_%04d.mat',cell_id));
if ~isfile(record_file) || ~isfile(seed_file)
    error('Run_Field_Cell_Spinup:TrackingNotReady', ...
        'Cell %d has no tracking record yet (%s); run the track stage first.', ...
        cell_id,record_file);
end
R = load(record_file,'record');
seed = load(seed_file,'tracking_key');
if ~strcmp(R.record.tracking_key,seed.tracking_key)
    error('Run_Field_Cell_Spinup:TrackingMismatch', ...
        'The record and the replay seed of cell %d come from different tracking runs.',cell_id);
end
extra = struct('LATERAL_REPLAY_FILE',record_file,'SPINUP_SEED_FILE',seed_file);
if F.warm.on
    %% the ring-end state of a tracking run that began from the earlier
    %% converged states: Phase 3 starts from it
    extra.SPINUP_WARM_START = true;
end
st = full_cell_spinup(replay_state_path(F.dir,cell_id),cell_id,F.case,F.inputs(c), ...
    replay_label(F,R.record.identity,seed.tracking_key),extra);
st = stamp(st,"replay",false);
end

function label = replay_label(F,record_identity,tracking_key)
%% A replay's provenance label, the same where it is run and where it is read.
label = sprintf('replay:%s;tracking:%s',record_identity,tracking_key);
if F.warm.on; label = ['warm-' label]; end
end

function f = replay_state_path(output_directory,cell_id)
f = fullfile(output_directory,sprintf('replay_%04d_BG.mat',cell_id));
end

function [out,method] = reconstruct_cells(F,topo,opt,c)
%% Lian et al. Fig. 3e-f: the forest over the replayed cells.  With a cell
%% index, that cell alone (a replayed cell is its own replay).
replayed = cell(numel(F.track_idx),1);
for s = 1:numel(F.track_idx)
    replayed{s} = load_replayed_state(F,F.track_idx(s));
end
replayed = harmonise_states(replayed);
handoff_file = fullfile(F.dir,'track','handoff_states.mat');
if ~isfile(handoff_file)
    error('Run_Field_Cell_Spinup:TrackingNotReady', ...
        '%s is missing; run the track stage first.',handoff_file);
end
H = load(handoff_file,'handoff_states','cell_ids');
if ~isequal(H.cell_ids(:),F.resolved.cell_id(:))
    error('Run_Field_Cell_Spinup:TrackingMismatch', ...
        'The tracking run''s handoff states belong to another grid.');
end
args = {'CaseNumber',F.case,'OutputDirectory',F.dir,'RFOptions',opt.RFOptions};
if ~isempty(c); args = [args {'CellIndex',c}]; end
[out,report] = Reconstruct_Cell_Pools_RF(topo,F.resolved,F.inputs,replayed, ...
    F.track_idx,H.handoff_states,args{:});
if isempty(c)
    method = report;
else
    method = report.method_per_cell(1);
end
end

function st = load_replayed_state(F,c)
%% A finished replay of tracked cell c, read and checked, never written.
cell_id = F.resolved.cell_id(c);
f = replay_state_path(F.dir,cell_id);
if ~isfile(f)
    error('Run_Field_Cell_Spinup:ReplayNotReady', ...
        'Tracked cell %d has no replayed state (%s); run the replay array first.',cell_id,f);
end
st = load_bg_state(f);
Validate_BG_State_For_Transient(st,f);
R = load(fullfile(F.dir,'track',sprintf('record_%04d.mat',cell_id)),'record');
seed = load(fullfile(F.dir,'track',sprintf('ringend_%04d.mat',cell_id)),'tracking_key');
want = Spinup_Provenance(cell_id,F.inputs(c).agronomic_identity,F.case, ...
    F.inputs(c).MOD_PARAM,F.inputs(c).ERW_PARAM, ...
    replay_label(F,R.record.identity,seed.tracking_key));
if ~(isfield(st,'spinup_provenance') && strcmp(char(string(st.spinup_provenance)),want))
    error('Run_Field_Cell_Spinup:ReplayStale', ...
        '%s was replayed for another run; rerun the replay array.',f);
end
st = stamp(st,"replay",false);
end

%% =======================================================================
function stop_on_failed_cells(report,resolved,state_file)
%% A cell that could not be spun up holds a neighbour's state flagged
%% unconverged; the transient loader rejects such a file, so the launcher
%% must not announce it as validated.  The file is kept for inspection.
if isfield(report,'n_failed') && report.n_failed > 0
    failed_ids = resolved.cell_id(report.method_per_cell == "failed");
    error('Run_Field_Cell_Spinup:CellsFailed', ...
        ['%d cell(s) could not be spun up (cell id(s) %s); %s holds their ' ...
         'flagged neighbour states and cannot be used by the transient run.'], ...
        report.n_failed,mat2str(failed_ids(:).'),state_file);
end
end

%% =======================================================================
function save_state_file(state_file,CellBGStates,resolved,rf_report) %#ok<INUSD>
output_directory = fileparts(state_file);
temporary_file = [tempname(output_directory) '.mat'];
cleanup = onCleanup(@() delete_if_present(temporary_file));
save(temporary_file,'CellBGStates','resolved','rf_report','-v7.3');
movefile(temporary_file,state_file,'f');
clear cleanup
end

%% =======================================================================
function plan_file = write_plan(output_directory,topo,resolved,method,track_idx,class_of,class_rep)
%% spinup_plan.csv (one row per cell) and the SLURM --array index lists.
N = topo.N_cells;
is_tracked = false(N,1); is_tracked(track_idx) = true;
is_rep = false(N,1); is_rep(class_rep) = true;
plan = table((1:N).',resolved.cell_id(:),is_tracked,class_of(:),is_rep, ...
    'VariableNames',{'cell_index','cell_id','is_tracked','class_index','is_class_seed'});
plan.method = repmat("full",N,1);
if strcmp(method,'flux_tracking')
    plan.method(is_tracked) = "replay";
    plan.method(~is_tracked) = "rf";
end
plan.spinup_method = repmat(string(method),N,1);
plan_file = fullfile(output_directory,'spinup_plan.csv');
writetable(plan,plan_file);
write_index_list(fullfile(output_directory,'array_indices_all.txt'),1:N);
if strcmp(method,'flux_tracking')
    write_index_list(fullfile(output_directory,'array_indices_classes.txt'),1:numel(class_rep));
    write_index_list(fullfile(output_directory,'array_indices_tracked.txt'),find(is_tracked));
    write_index_list(fullfile(output_directory,'array_indices_rf.txt'),find(~is_tracked));
    fprintf(['Run_Field_Cell_Spinup: plan for %d cells (flux_tracking: %d class ' ...
        'seed(s), %d tracked and replayed, %d from the forest) written to %s\n'], ...
        N,numel(class_rep),nnz(is_tracked),nnz(~is_tracked),plan_file);
else
    fprintf('Run_Field_Cell_Spinup: plan for %d cells (per_cell) written to %s\n', ...
        N,plan_file);
end
end

function write_index_list(file_name,idx)
idx = idx(:).';
if isempty(idx)
    txt = '';
elseif isequal(idx,idx(1):idx(end))
    txt = sprintf('%d-%d',idx(1),idx(end));
else
    txt = strjoin(arrayfun(@(i) sprintf('%d',i),idx,'UniformOutput',false),',');
end
fid = fopen(file_name,'w'); fprintf(fid,'%s\n',txt); fclose(fid);
end

%% =======================================================================
function file_name = cell_state_path(output_directory,cell_id)
file_name = fullfile(output_directory,sprintf('cellstate_%04d.mat',cell_id));
end

function file_name = spinup_state_path(output_directory,cell_id)
%% A per_cell spin-up's own output file for a cell (BG_OUTPUT_FILE).
file_name = fullfile(output_directory,sprintf('cell_%04d_BG.mat',cell_id));
end

function file_name = write_cell_state(output_directory,cell_index,cell_id,st,method) %#ok<INUSD>
file_name = cell_state_path(output_directory,cell_id);
temporary_file = [tempname(output_directory) '.mat'];
cleanup = onCleanup(@() delete_if_present(temporary_file));
method = string(method); %#ok<NASGU>
save(temporary_file,'st','method','cell_id','cell_index','-v7.3');
movefile(temporary_file,file_name,'f');
clear cleanup
end

%% =======================================================================
function [CellBGStates,report] = assemble_cell_states(output_directory, ...
    topo,resolved,MOD,ERW,method)
%% Every cell's cellstate_<id>.mat, validated, identity-checked against
%% the resolved inputs, and made a uniform struct array.
N = topo.N_cells;
states = cell(N,1); methods = strings(N,1);
%% completeness first: an incomplete array is the common case and the
%% message should say which cells are missing, not fail on a state
present = arrayfun(@(c) isfile(cell_state_path(output_directory,resolved.cell_id(c))),1:N);
missing = resolved.cell_id(~present);
if ~isempty(missing)
    error('Run_Field_Cell_Spinup:CellStatesMissing', ...
        '%d of %d cell states are missing (cell id(s) %s); the array is not complete.', ...
        numel(missing),N,mat2str(missing(:).'));
end
%% A failed reconstruction is named as such before validation rejects it
%% for the unconverged neighbour state it holds.
failed_ids = [];
for c = 1:N
    S = load(cell_state_path(output_directory,resolved.cell_id(c)),'method');
    if isfield(S,'method') && string(S.method) == "failed"
        failed_ids(end+1) = resolved.cell_id(c); %#ok<AGROW>
    end
end
if ~isempty(failed_ids)
    error('Run_Field_Cell_Spinup:CellsFailed', ...
        ['%d cell(s) could not be spun up (ids %s): their reconstruction and cold ' ...
         'start both failed. Rerun those array tasks after fixing the cause.'], ...
        numel(failed_ids),mat2str(failed_ids));
end
for c = 1:N
    f = cell_state_path(output_directory,resolved.cell_id(c));
    S = load(f,'st','method','cell_id');
    if S.cell_id ~= resolved.cell_id(c) || ~isfield(S.st,'cell_id') || ...
            S.st.cell_id ~= resolved.cell_id(c)
        error('Run_Field_Cell_Spinup:CellStateIdentity', ...
            '%s holds cell %d, expected %d.',f,S.st.cell_id,resolved.cell_id(c));
    end
    Validate_BG_State_For_Transient(S.st,f);
    states{c} = S.st; methods(c) = string(S.method);
end
CellBGStates = harmonise_states(states);
%% The identity check the transient run will make, made now.
Configure_Cell_Agronomic_Inputs(MOD,ERW,resolved,CellBGStates);
failed = arrayfun(@(st) isfield(st,'reconstruction_failed') && ...
    ~isempty(st.reconstruction_failed) && logical(st.reconstruction_failed),CellBGStates);
report = struct('method',method,'assembled',true, ...
    'method_per_cell',methods,'n_full',nnz(methods == "full"), ...
    'n_replay',nnz(methods == "replay"),'n_rf',nnz(methods == "rf"), ...
    'n_fallback',nnz(methods == "fallback"),'n_failed',nnz(failed));
if any(failed)
    warning('Run_Field_Cell_Spinup:AssembledFailures', ...
        '%d assembled cell(s) are flagged reconstruction_failed (ids %s).', ...
        nnz(failed),mat2str(resolved.cell_id(failed).'));
end
fprintf(['Run_Field_Cell_Spinup: assembled %d cells (%d full, %d replayed, %d from ' ...
    'the forest, %d fallback, %d failed)\n'],N,report.n_full,report.n_replay, ...
    report.n_rf,report.n_fallback,report.n_failed);
end

function states = harmonise_states(cells)
%% States written by different paths may differ in optional fields; give
%% every one the union so they form one struct array.
names = {};
for c = 1:numel(cells); names = union(names,fieldnames(cells{c}),'stable'); end
for c = 1:numel(cells)
    for k = 1:numel(names)
        if ~isfield(cells{c},names{k}); cells{c}.(names{k}) = []; end
    end
    cells{c} = orderfields(cells{c},names);
end
states = vertcat(cells{:});
end

%% =======================================================================
function st = stamp(st,method,failed)
st.reconstruction_method = string(method);
st.reconstruction_failed = logical(failed);
end

%% =======================================================================
function st = full_cell_spinup(cell_file,cell_id,case_number,cell_inputs,seed_label,extra)
%% One cell's plot-scale spin-up into cell_file.  What the run is -- the
%% cell, its inputs, the case and its rows (with the forcing, meteo and
%% transient windows), and what it starts from ('cold', or the lateral
%% record and seed it replays) -- is its provenance key: a finished state is
%% reused, and a phase checkpoint resumed, only under the same key, so an
%% interrupted array task resumes and any change starts afresh.
if nargin < 5 || isempty(seed_label); seed_label = 'cold'; end
if nargin < 6; extra = struct(); end
provenance = Spinup_Provenance(cell_id,cell_inputs.agronomic_identity, ...
    case_number,cell_inputs.MOD_PARAM,cell_inputs.ERW_PARAM,seed_label);
if isfile(cell_file)
    try
        st = load_bg_state(cell_file);
        Validate_BG_State_For_Transient(st,cell_file);
        if isfield(st,'cell_id') && isfield(st,'spinup_provenance') && ...
                st.cell_id == cell_id && ...
                strcmp(char(string(st.spinup_provenance)),provenance)
            fprintf('  cell %d: reusing validated spin-up state %s\n',cell_id,cell_file);
            return
        end
        fprintf('  cell %d: existing state was made for another run; spinning up afresh\n',cell_id);
    catch ME
        fprintf('  cell %d: existing state unusable (%s); spinning up afresh\n',cell_id,ME.message);
    end
    delete_if_present(cell_file);
end
settings = struct('NN_case_study_ERW_MOD',case_number, ...
    'RUN_MODE','spinup','RUN_CHEMISTRY','erw','RUN_DOMAIN','single', ...
    'SAVE_RESULTS',0,'MOD_PARAM_override',cell_inputs.MOD_PARAM, ...
    'ERW_PARAM_override',cell_inputs.ERW_PARAM, ...
    'BG_OUTPUT_FILE',cell_file,'CELL_ID',cell_id, ...
    'CELL_AGRONOMIC_IDENTITY',cell_inputs.agronomic_identity, ...
    'SPINUP_PROVENANCE',provenance);
for name = fieldnames(extra).'
    settings.(name{1}) = extra.(name{1});
end
Run_Automatic_Site_Spinup(settings);
st = load_bg_state(cell_file);
Validate_BG_State_For_Transient(st,cell_file);
end

function check_plan_method(output_directory,method)
%% One method per output directory: per_cell and flux_tracking runs would
%% otherwise overwrite each other's cellstate files and the assembly would
%% mix them.
plan_file = fullfile(output_directory,'spinup_plan.csv');
if ~isfile(plan_file)
    return
end
plan = readtable(plan_file,'TextType','string');
if ~ismember('spinup_method',plan.Properties.VariableNames)
    return
end
planned = plan.spinup_method(1);
if ~strcmpi(planned,method)
    error('Run_Field_Cell_Spinup:PlanMethodMismatch', ...
        ['%s plans Method=%s; this step runs Method=%s. Give each method its ' ...
         'own output directory (TC_OUTPUT_DIR).'],plan_file,planned,method);
end
end

function delete_if_present(file_name)
if isfile(file_name)
    delete(file_name)
end
end

%% =======================================================================
function r = dem_resolution(requested,topo)
%% Default the DEM spacing to the process-grid cell size rather than the
%% tile's native 1 m.  Terrain attributes are sampled at cell centres, so
%% a DEM far finer than the cell aliases sub-cell microtopography into
%% slope / curvature / TWI (and slope, a gradient magnitude, is biased
%% one-sided by it).  Load_Field_DEM area-averages to whatever spacing it
%% is given.
r = requested;
if isempty(r) && isfield(topo,'cell_size_m') && ~isempty(topo.cell_size_m)
    r = topo.cell_size_m;
end
end

%% =======================================================================
function W = warm_seeds(states_file,resolved,output_directory)
%% A warm start's seeds: every cell's converged state from an earlier field
%% spin-up (its CellBGStates.mat), matched by cell id and checked to lie at
%% the same place (the identity records the cell centre), so states of
%% another grid are refused.  Kept in memory; a stage writes only its own
%% cell's seed file (cell_start).
W = struct('on',false,'states',struct([]),'dir','');
if isempty(states_file); return; end
states_file = char(states_file);
if ~isfile(states_file)
    error('Run_Field_Cell_Spinup:WarmStatesMissing', ...
        'Warm-start states not found: %s',states_file);
end
%% The warm run writes cell files and a CellBGStates.mat of its own, which
%% would overwrite the run it starts from.
if strcmpi(canonical_path(fileparts(states_file)),canonical_path(output_directory))
    error('Run_Field_Cell_Spinup:WarmStartSameDirectory', ...
        ['A warm start needs its own output directory, not %s, which holds ' ...
         'the states it starts from.'],output_directory);
end
S = load(states_file,'CellBGStates');
prev = S.CellBGStates;
ids = arrayfun(@(s) double(s.cell_id),prev);
order = zeros(height(resolved),1);
for c = 1:height(resolved)
    k = find(ids == resolved.cell_id(c),1);
    if isempty(k)
        error('Run_Field_Cell_Spinup:WarmStateMissing', ...
            'The warm-start states hold no cell %d.',resolved.cell_id(c));
    end
    xy = [];
    if isfield(prev,'cell_agronomic_identity')
        xy = sscanf(regexp(char(string(prev(k).cell_agronomic_identity)), ...
            'x=[^;]*;y=[^;]*','match','once'),'x=%f;y=%f');
    end
    centre = [resolved.cell_x_projected_m(c) resolved.cell_y_projected_m(c)];
    if numel(xy) ~= 2 || any(abs(xy(:).' - centre) > 1e-6)
        error('Run_Field_Cell_Spinup:WarmStateGrid', ...
            ['The warm-start state of cell %d does not lie at this cell''s centre ' ...
             '(another grid, or no cell identity).'],resolved.cell_id(c));
    end
    order(c) = k;
end
W = struct('on',true,'states',prev(order),'dir',fullfile(output_directory,'warm'));
end

function p = canonical_path(p)
if isempty(p); p = pwd; end
p = char(java.io.File(char(p)).getCanonicalPath());
end

function [label,extra] = cell_start(F,c)
%% How a cell's own spin-up starts: cold, or warm from its earlier state.
label = 'cold'; extra = struct();
if ~F.warm.on; return; end
st = F.warm.states(c);
cell_id = F.resolved.cell_id(c);
if ~isfolder(F.warm.dir); mkdir(F.warm.dir); end
f = fullfile(F.warm.dir,sprintf('warm_%04d.mat',cell_id));
tmp = [tempname(F.warm.dir) '.mat'];
save(tmp,'-struct','st','-v7.3');
movefile(tmp,f,'f');
earlier = mat2str(st.B_out(end,:),17);
if isfield(st,'spinup_provenance'); earlier = char(string(st.spinup_provenance)); end
label = ['warm:' earlier];
extra = struct('SPINUP_SEED_FILE',f,'SPINUP_WARM_START',true);
end

%% =======================================================================
function names = column_names(cols,option)
if ~isstruct(cols) || ~isscalar(cols)
    error(['Run_Field_Cell_Spinup:Bad' option], ...
        '%s must be a scalar struct of column values.',option);
end
names = fieldnames(cols);
end

function check_columns(names,known,per_cell,option)
%% A column the spin-up never reads (a misspelling) would otherwise be
%% accepted silently and the row would run without it; one the cells set
%% from the agronomic file would be overwritten, or overwrite them.
unknown = names(~ismember(names,known));
if ~isempty(unknown)
    error(['Run_Field_Cell_Spinup:Unknown' option(1:end-7) 'Column'], ...
        '%s %s: not a workbook column nor a column the spin-up reads.', ...
        option,strjoin(unknown,', '));
end
own = names(ismember(names,per_cell));
if ~isempty(own)
    error('Run_Field_Cell_Spinup:CellColumn', ...
        '%s %s: set per cell from the agronomic file, not for the whole field.', ...
        option,strjoin(own,', '));
end
end

function v = column_value(v)
if ischar(v) || iscellstr(v); v = string(v); end   % text columns are strings
end

function tf = lacks_exchanger_seed(E)
tf = false;
if ~(ismember('fX_conc_convention',E.Properties.VariableNames) && ...
        string(E.fX_conc_convention) == "f_CEC_inventory")
    return
end
f = E.f_CEC_in;
if isstring(f) || ischar(f); f = str2num(char(f)); end %#ok<ST2NM>
tf = ~(isnumeric(f) && numel(f) == 6 && all(isfinite(f)));
end
