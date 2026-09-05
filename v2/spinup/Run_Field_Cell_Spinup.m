function [state_file,CellInputs] = Run_Field_Cell_Spinup(shape_file,agronomic_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one independent automatic spin-up for every cell in one polygon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'GridCellSizeM',[],'CaseNumber',5,'OutputFile','','DryRun',false, ...
    'Method','per_cell','FieldDEM','','SubsampleFraction',0.12, ...
    'MinCoverage',[],'MODParamColumns',[]);
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
if ~exist(output_directory,'dir')
    mkdir(output_directory);
end

if ~any(strcmpi(opt.Method,{'per_cell','rf'}))
    error('Run_Field_Cell_Spinup:BadMethod', ...
        'Method must be ''per_cell'' or ''rf''.');
end

run_cell_spinup = @(c) full_cell_spinup(output_directory,resolved.cell_id(c), ...
    opt.CaseNumber,CellInputs(c));

if strcmpi(opt.Method,'rf')
    if exist('TreeBagger','file') ~= 2
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

if strcmpi(opt.Method,'rf')
    dem = Load_Field_DEM(char(opt.FieldDEM),topo, ...
        'OutResolutionM',dem_resolution(opt.GridCellSizeM,topo));
    topo = Attach_Topography_To_Topology(topo,dem);
    [sub_idx,~] = Select_Topographic_Subsample(topo,resolved,opt.SubsampleFraction);
    sub_states = repmat(run_cell_spinup(sub_idx(1)), numel(sub_idx), 1);
    for s = 2:numel(sub_idx); sub_states(s) = run_cell_spinup(sub_idx(s)); end
    [CellBGStates,rf_report] = Reconstruct_Cell_Pools_RF(topo,resolved, ...
        CellInputs,sub_states,sub_idx,'CaseNumber',opt.CaseNumber, ...
        'OutputDirectory',output_directory);
    resolved.spinup_method = repmat("rf",topo.N_cells,1);
else
    state_cells = cell(topo.N_cells,1);
    for c = 1:topo.N_cells
        state_cells{c} = run_cell_spinup(c);
    end
    CellBGStates = vertcat(state_cells{:});
    rf_report = struct('method','per_cell'); %#ok<NASGU>
end

temporary_file = [tempname(output_directory) '.mat'];
cleanup = onCleanup(@() delete_if_present(temporary_file));
save(temporary_file,'CellBGStates','resolved','rf_report','-v7.3');
movefile(temporary_file,state_file,'f');
clear cleanup
fprintf('Saved %d validated cell spin-up states (%s): %s\n', ...
    topo.N_cells,opt.Method,state_file);
end

%% =======================================================================
function st = full_cell_spinup(output_directory, cell_id, case_number, cell_inputs)
cell_file = fullfile(output_directory,sprintf('cell_%03d_BG.mat',cell_id));
settings = struct('NN_case_study_ERW_MOD',case_number, ...
    'RUN_MODE','spinup','RUN_CHEMISTRY','erw','RUN_DOMAIN','single', ...
    'SAVE_RESULTS',0,'MOD_PARAM_override',cell_inputs.MOD_PARAM, ...
    'ERW_PARAM_override',cell_inputs.ERW_PARAM, ...
    'BG_OUTPUT_FILE',cell_file,'CELL_ID',cell_id, ...
    'CELL_AGRONOMIC_IDENTITY',cell_inputs.agronomic_identity);
Run_Automatic_Site_Spinup(settings);
st = load_bg_state(cell_file);
Validate_BG_State_For_Transient(st,cell_file);
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
