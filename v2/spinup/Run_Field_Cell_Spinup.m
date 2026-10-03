function [state_file,CellInputs] = Run_Field_Cell_Spinup(shape_file,agronomic_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one independent automatic spin-up for every cell in one polygon.
%%
%% Cluster modes (one SLURM array task per cell, see v2/cluster):
%%   'PlanOnly',true      write spinup_plan.csv and the array-index lists
%%                        (all cells; the RF subsample; the RF remainder)
%%   'CellIndex',k        spin up ONLY cell k (index into the polygon's
%%                        cells): a full spin-up for Method per_cell or an
%%                        RF-subsample cell, the RF-seeded short spin-up
%%                        otherwise (which needs every subsample cell's
%%                        full state on disk already); writes
%%                        cellstate_<id>.mat and returns its path
%%   'AssembleOnly',true  gather every cellstate_<id>.mat into the
%%                        validated CellBGStates.mat the transient run reads
%% A finished cell_<id>_BG.mat in the output directory is reused rather
%% than recomputed, so an interrupted run or a cluster array resumes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'GridCellSizeM',[],'CaseNumber',5,'OutputFile','','DryRun',false, ...
    'Method','per_cell','FieldDEM','','SubsampleFraction',0.12, ...
    'MinCoverage',[],'MODParamColumns',[],'CellIndex',[],'PlanOnly',false, ...
    'AssembleOnly',false,'RFOptions',{{}});
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
CellInputs = Configure_Cell_Agronomic_Inputs(MOD,ERW,resolved);
if ~isempty(opt.MODParamColumns)
    %% Optional MOD_PARAM columns the workbook row does not carry (the
    %% spin-up meteo and transient windows, a forcing file), applied to
    %% every cell's row; the spin-up reads them through table_value.
    if ~isstruct(opt.MODParamColumns) || ~isscalar(opt.MODParamColumns)
        error('Run_Field_Cell_Spinup:BadMODParamColumns', ...
            'MODParamColumns must be a scalar struct of column values.');
    end
    names = fieldnames(opt.MODParamColumns);
    %% A column the spin-up never reads (a misspelling) would otherwise be
    %% accepted silently and the row would run without it.
    known = [MOD.Properties.VariableNames, {'id_location','forcing_file', ...
        'transient_start','transient_end','scenario_transient','source_data_file', ...
        'meteo_file','rain_chem_file','meteo_start','meteo_end','rain_chem_start', ...
        'rain_chem_end'}];
    unknown = names(~ismember(names,known));
    if ~isempty(unknown)
        error('Run_Field_Cell_Spinup:UnknownMODParamColumn', ...
            'MODParamColumns %s: not a workbook column nor a spin-up column the spin-up reads.', ...
            strjoin(unknown,', '));
    end
    for c = 1:numel(CellInputs)
        for k = 1:numel(names)
            v = opt.MODParamColumns.(names{k});
            if ischar(v) || iscellstr(v); v = string(v); end   % text columns are strings
            CellInputs(c).MOD_PARAM.(names{k}) = v;
        end
    end
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

if ~any(strcmpi(opt.Method,{'per_cell','rf'}))
    error('Run_Field_Cell_Spinup:BadMethod', ...
        'Method must be ''per_cell'' or ''rf''.');
end

run_cell_spinup = @(c) full_cell_spinup(output_directory,resolved.cell_id(c), ...
    opt.CaseNumber,CellInputs(c));

cluster_step = ~isempty(opt.CellIndex) || opt.PlanOnly || opt.AssembleOnly;
if strcmpi(opt.Method,'rf')
    if exist('TreeBagger','file') ~= 2 && cluster_step
        %% An array step cannot switch method on its own: the plan and the
        %% other tasks would still be running the RF path.
        error('Run_Field_Cell_Spinup:NoRFToolbox', ...
            ['Method=''rf'' needs TreeBagger (Statistics and Machine Learning ' ...
             'Toolbox) on every node that runs a step; use Method=''per_cell'' ' ...
             'for the whole array instead.']);
    elseif exist('TreeBagger','file') ~= 2
        warning('Run_Field_Cell_Spinup:NoRFToolbox', ...
            ['TreeBagger unavailable; falling back to Method=''per_cell''. ' ...
             'Install the Statistics and Machine Learning Toolbox for the ' ...
             'RF reconstruction.']);
        opt.Method = 'per_cell';
    elseif isempty(opt.FieldDEM)
        error('Run_Field_Cell_Spinup:RFNeedsDEM', ...
            'Method=''rf'' needs a FieldDEM to attach topography for the RF covariates.');
    end
end

%% The RF subsample is decided once and deterministically (same seed,
%% same grid), so the plan, every array task and the assembly agree on
%% which cells are spun up in full.
sub_idx = [];
if strcmpi(opt.Method,'rf')
    dem = Load_Field_DEM(char(opt.FieldDEM),topo, ...
        'OutResolutionM',dem_resolution(opt.GridCellSizeM,topo));
    topo = Attach_Topography_To_Topology(topo,dem);
    [sub_idx,~] = Select_Topographic_Subsample(topo,resolved,opt.SubsampleFraction);
    sub_idx = sub_idx(:);
    %% The forest's k-fold estimate needs at least 2 x KFold training
    %% cells; find that out now, before any full spin-up is spent.
    kfold = 5;
    for i = 1:2:numel(opt.RFOptions)
        if strcmpi(char(opt.RFOptions{i}),'KFold'); kfold = opt.RFOptions{i+1}; end
    end
    if numel(sub_idx) < 2*kfold
        error('Run_Field_Cell_Spinup:SubsampleTooSmall', ...
            ['SubsampleFraction %.2f gives %d subsample cells on %d; the forest ' ...
             'needs at least 2 x KFold = %d. Raise SubsampleFraction or pass ' ...
             'RFOptions {KFold,k}.'],opt.SubsampleFraction,numel(sub_idx), ...
            topo.N_cells,2*kfold);
    end
end

if cluster_step
    check_plan_method(output_directory,opt.Method);
end
if opt.PlanOnly
    state_file = write_plan(output_directory,topo,resolved,opt.Method,sub_idx);
    return
end

if opt.AssembleOnly
    [CellBGStates,rf_report] = assemble_cell_states(output_directory,topo, ...
        resolved,MOD,ERW,opt.Method);
    resolved.spinup_method = repmat(string(opt.Method),topo.N_cells,1);
    save_state_file(state_file,CellBGStates,resolved,rf_report);
    stop_on_failed_cells(rf_report,resolved,state_file);
    fprintf('Saved %d assembled cell spin-up states (%s): %s\n', ...
        topo.N_cells,opt.Method,state_file);
    return
end

if ~isempty(opt.CellIndex)
    k = opt.CellIndex;
    if ~(isnumeric(k) && isscalar(k) && isfinite(k) && k == fix(k) && ...
            k >= 1 && k <= topo.N_cells)
        error('Run_Field_Cell_Spinup:BadCellIndex', ...
            'CellIndex must be an integer in 1..%d.',topo.N_cells);
    end
    if strcmpi(opt.Method,'rf') && ~ismember(k,sub_idx)
        %% A worker never runs a full spin-up it was not assigned: every
        %% subsample cell's full state must already be on disk (the
        %% subsample array tasks), and is reused from there.
        ready = arrayfun(@(c) isfile(spinup_state_path(output_directory, ...
            resolved.cell_id(c))),sub_idx);
        if ~all(ready)
            error('Run_Field_Cell_Spinup:SubsampleNotReady', ...
                ['RF cell %d needs the full spin-up state of every subsample ' ...
                 'cell; missing cell id(s) %s. Run the subsample array first.'], ...
                k,mat2str(resolved.cell_id(sub_idx(~ready)).'));
        end
        %% Read only: hundreds of RF tasks share these files, and none of
        %% them may re-spin or delete a subsample cell.
        sub_states = repmat(load_finished_state(output_directory,resolved, ...
            CellInputs,sub_idx(1),opt.CaseNumber),numel(sub_idx),1);
        for si = 2:numel(sub_idx)
            sub_states(si) = load_finished_state(output_directory,resolved, ...
                CellInputs,sub_idx(si),opt.CaseNumber);
        end
        [st,rep] = Reconstruct_Cell_Pools_RF(topo,resolved,CellInputs, ...
            sub_states,sub_idx,'CaseNumber',opt.CaseNumber, ...
            'OutputDirectory',output_directory,'CellIndex',k, ...
            'RFOptions',opt.RFOptions);
        method = rep.method_per_cell(1);
    else
        st = stamp(run_cell_spinup(k),"full",false);
        method = "full";
    end
    state_file = write_cell_state(output_directory,k,resolved.cell_id(k),st,method);
    fprintf('Saved cell %d (id %d, %s): %s\n',k,resolved.cell_id(k),method,state_file);
    if method == "failed"
        %% The file holds a flagged neighbour state for the record; the
        %% task itself must fail so the array shows it.
        error('Run_Field_Cell_Spinup:CellFailed', ...
            ['Cell %d (id %d): the RF seed and the cold start both failed; see ' ...
             'the warnings above.'],k,resolved.cell_id(k));
    end
    return
end

if strcmpi(opt.Method,'rf')
    sub_states = repmat(run_cell_spinup(sub_idx(1)),numel(sub_idx),1);
    for si = 2:numel(sub_idx); sub_states(si) = run_cell_spinup(sub_idx(si)); end
    [CellBGStates,rf_report] = Reconstruct_Cell_Pools_RF(topo,resolved, ...
        CellInputs,sub_states,sub_idx,'CaseNumber',opt.CaseNumber, ...
        'OutputDirectory',output_directory,'RFOptions',opt.RFOptions);
    resolved.spinup_method = repmat("rf",topo.N_cells,1);
else
    state_cells = cell(topo.N_cells,1);
    for c = 1:topo.N_cells
        state_cells{c} = stamp(run_cell_spinup(c),"full",false);
    end
    CellBGStates = vertcat(state_cells{:});
    rf_report = struct('method','per_cell');
end
save_state_file(state_file,CellBGStates,resolved,rf_report);
stop_on_failed_cells(rf_report,resolved,state_file);
fprintf('Saved %d validated cell spin-up states (%s): %s\n', ...
    topo.N_cells,opt.Method,state_file);
end

%% =======================================================================
function stop_on_failed_cells(rf_report,resolved,state_file)
%% A cell that could not be spun up holds a neighbour's state flagged
%% unconverged; the transient loader rejects such a file, so the launcher
%% must not announce it as validated.  The file is kept for inspection.
if isfield(rf_report,'n_failed') && rf_report.n_failed > 0
    failed_ids = resolved.cell_id(rf_report.method_per_cell == "failed");
    error('Run_Field_Cell_Spinup:CellsFailed', ...
        ['%d cell(s) could not be spun up (cell id(s) %s); %s holds their ' ...
         'flagged neighbour states and cannot be used by the transient run.'], ...
        rf_report.n_failed,mat2str(failed_ids(:).'),state_file);
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
function plan_file = write_plan(output_directory,topo,resolved,method,sub_idx)
%% spinup_plan.csv (one row per cell) and the SLURM --array index lists:
%% all cells, the RF subsample (full spin-ups), the RF remainder
%% (seeded short spin-ups, run after the subsample array has finished).
N = topo.N_cells;
is_sub = false(N,1); is_sub(sub_idx) = true;
plan = table((1:N).',resolved.cell_id(:),is_sub, ...
    'VariableNames',{'cell_index','cell_id','is_subsample'});
plan.method = repmat("full",N,1);
if strcmpi(method,'rf'); plan.method(~is_sub) = "rf_seed"; end
plan.spinup_method = repmat(string(lower(method)),N,1);
plan_file = fullfile(output_directory,'spinup_plan.csv');
writetable(plan,plan_file);
write_index_list(fullfile(output_directory,'array_indices_all.txt'),1:N);
write_index_list(fullfile(output_directory,'array_indices_subsample.txt'),find(is_sub));
write_index_list(fullfile(output_directory,'array_indices_rf.txt'),find(~is_sub));
fprintf(['Run_Field_Cell_Spinup: plan for %d cells (%s; %d full-spin-up ' ...
    'subsample cells) written to %s\n'],N,method,nnz(is_sub),plan_file);
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
%% The spin-up's own output file for a cell (BG_OUTPUT_FILE).
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
function [CellBGStates,rf_report] = assemble_cell_states(output_directory, ...
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
        ['%d cell(s) could not be spun up (ids %s): their RF seed and cold ' ...
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
rf_report = struct('method',method,'assembled',true, ...
    'method_per_cell',methods,'n_full',nnz(methods == "full"), ...
    'n_rf_seed',nnz(methods == "rf_seed"),'n_fallback',nnz(methods == "fallback"), ...
    'n_failed',nnz(failed));
if any(failed)
    warning('Run_Field_Cell_Spinup:AssembledFailures', ...
        '%d assembled cell(s) are flagged reconstruction_failed (ids %s).', ...
        nnz(failed),mat2str(resolved.cell_id(failed).'));
end
fprintf('Run_Field_Cell_Spinup: assembled %d cells (%d full, %d RF-seeded, %d fallback, %d failed)\n', ...
    N,rf_report.n_full,rf_report.n_rf_seed,rf_report.n_fallback,rf_report.n_failed);
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
function st = full_cell_spinup(output_directory, cell_id, case_number, cell_inputs)
cell_file = spinup_state_path(output_directory,cell_id);
%% What this run is: the cell, its inputs, the case and its rows (with the
%% forcing, meteo and transient windows), from a cold start.  A finished
%% state is reused, and a phase checkpoint resumed, only under the same key,
%% so an interrupted array task resumes and any change starts afresh.
provenance = Spinup_Provenance(cell_id,cell_inputs.agronomic_identity, ...
    case_number,cell_inputs.MOD_PARAM,cell_inputs.ERW_PARAM,'cold');
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
Run_Automatic_Site_Spinup(settings);
st = load_bg_state(cell_file);
Validate_BG_State_For_Transient(st,cell_file);
end

function st = load_finished_state(output_directory,resolved,CellInputs,c,case_number)
%% A finished full spin-up of subsample cell c, read and checked, never
%% written: it must be the cold-start run of this cell's current inputs.
cell_id = resolved.cell_id(c);
cell_file = spinup_state_path(output_directory,cell_id);
provenance = Spinup_Provenance(cell_id,CellInputs(c).agronomic_identity, ...
    case_number,CellInputs(c).MOD_PARAM,CellInputs(c).ERW_PARAM,'cold');
st = load_bg_state(cell_file);
Validate_BG_State_For_Transient(st,cell_file);
if ~(isfield(st,'spinup_provenance') && ...
        strcmp(char(string(st.spinup_provenance)),provenance))
    error('Run_Field_Cell_Spinup:SubsampleStale', ...
        ['%s was spun up for another run (inputs, case or windows); rerun the ' ...
         'subsample array before the RF array.'],cell_file);
end
end

function check_plan_method(output_directory,method)
%% One method per output directory: per_cell and rf runs would otherwise
%% overwrite each other's cellstate files and the assembly would mix them.
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
