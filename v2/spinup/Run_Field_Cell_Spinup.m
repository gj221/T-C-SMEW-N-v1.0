function [state_file,CellInputs] = Run_Field_Cell_Spinup(shape_file,agronomic_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one independent automatic spin-up for every cell in one polygon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'GridCellSizeM',[],'CaseNumber',5,'OutputFile','','DryRun',false);
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
    opt.FieldFeatureIndices);
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
MOD = readtable(workbook,'Sheet','MOD_PARAM','VariableNamingRule','preserve');
ERW = readtable(workbook,'Sheet','ERW_MOD','VariableNamingRule','preserve');
MOD = MOD(MOD.Number == opt.CaseNumber,:);
ERW = ERW(ERW.Number == opt.CaseNumber,:);
if height(MOD) ~= 1 || height(ERW) ~= 1
    error('Run_Field_Cell_Spinup:CaseSelection', ...
        'CaseNumber %d must select one MOD and one ERW row.',opt.CaseNumber);
end
CellInputs = Configure_Cell_Agronomic_Inputs(MOD,ERW,resolved);

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

state_cells = cell(topo.N_cells,1);
for c = 1:topo.N_cells
    cell_file = fullfile(output_directory,sprintf('cell_%03d_BG.mat', ...
        resolved.cell_id(c)));
    settings = struct('NN_case_study_ERW_MOD',opt.CaseNumber, ...
        'RUN_MODE','spinup','RUN_CHEMISTRY','erw','RUN_DOMAIN','single', ...
        'SAVE_RESULTS',0,'MOD_PARAM_override',CellInputs(c).MOD_PARAM, ...
        'ERW_PARAM_override',CellInputs(c).ERW_PARAM, ...
        'BG_OUTPUT_FILE',cell_file,'CELL_ID',resolved.cell_id(c), ...
        'CELL_AGRONOMIC_IDENTITY',CellInputs(c).agronomic_identity);
    Run_Automatic_Site_Spinup(settings);
    state_cells{c} = load_bg_state(cell_file);
    Validate_BG_State_For_Transient(state_cells{c},cell_file);
end
CellBGStates = vertcat(state_cells{:});

temporary_file = [tempname(output_directory) '.mat'];
cleanup = onCleanup(@() delete_if_present(temporary_file));
save(temporary_file,'CellBGStates','resolved','-v7.3');
movefile(temporary_file,state_file,'f');
clear cleanup
fprintf('Saved %d validated cell spin-up states: %s\n',topo.N_cells,state_file);
end

function delete_if_present(file_name)
if isfile(file_name)
    delete(file_name)
end
end
