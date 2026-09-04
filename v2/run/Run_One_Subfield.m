function result_file = Run_One_Subfield(run_id, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_One_Subfield
%%
%% Run one of the twelve simulations defined by Run_Matrix_12 and save
%% the result. Runs 1 to 9 use the retained sub-field domain. Runs 10 to
%% 12 use one plot-scale column and the original vertical transport only.
%%
%% This is the 2-D counterpart of ERW_run: it performs the same workbook
%% read, spin-up load, forcing load and window slice, then partitions the
%% plot into sub-plot cells and hands the whole set to
%% MAIN_FRAME_BG_ERW_2D instead of the single-column frame.
%%
%% TREATMENT PLACEMENT
%%   Each cell's ColCfg entry carries apply_rock, set from the run's
%%   layout and the cell's region.  Cells that receive no amendment still
%%   run the full carbonate and cation-exchange chemistry, so an
%%   untreated field interior and a treated margin are directly
%%   comparable inside one simulation.
%%
%% Usage
%%   Run_One_Subfield(5)                       % run 5 of the matrix
%%   Run_One_Subfield(11)                      % basalt, one-dimensional
%%   Run_One_Subfield(5,'NDays',3)             % short smoke test
%%   Run_One_Subfield(12,'StartDate','2016-11-09','NDays',3)
%%   Run_One_Subfield(5,'GridNx',1,'GridNy',1) % coarser main grid
%%   Run_One_Subfield(11,'AppeloAlpha',0.8)     % calibration trajectory
%%   Run_One_Subfield(11,'AppeloLogKH',[5.9 4.8 10.47 10.986 16.47])
%%
%% Returns the path of the saved result file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ResultDirectory defaults to the usual SIM/Subfield_12, so behaviour is
%% unchanged when it is not given. It exists because two runs of the same
%% run_id launched concurrently would otherwise race on one output path.
opt = struct('NDays',[],'StartDate',[],'GridNx',[],'GridNy',[], ...
    'FieldShapefile','','GridCellSizeM',[], ...
    'FieldFeatureIndices',[],'FieldPolygonIndex',[], ...
    'CellAgronomicFile','','CellSpinupStateFile','', ...
    'FieldDEM','','DEMResolutionM',[], ...
    'AppeloAlpha',[],'AppeloLogKH',[],'SaveResults',true,'ResultDirectory','');
if mod(numel(varargin),2) ~= 0
    error('Run_One_Subfield:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_One_Subfield:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
validate_optional_integer(opt.NDays,'NDays');
validate_optional_integer(opt.GridNx,'GridNx');
validate_optional_integer(opt.GridNy,'GridNy');
validate_optional_positive_scalar(opt.GridCellSizeM,'GridCellSizeM');
validate_optional_index_vector(opt.FieldFeatureIndices,'FieldFeatureIndices');
validate_optional_integer(opt.FieldPolygonIndex,'FieldPolygonIndex');
validate_optional_nonnegative_scalar(opt.AppeloAlpha,'AppeloAlpha');
if ~isempty(opt.AppeloLogKH) && ~(isnumeric(opt.AppeloLogKH) && ...
        isvector(opt.AppeloLogKH) && numel(opt.AppeloLogKH) == 5 && ...
        all(isfinite(opt.AppeloLogKH)))
    error('Run_One_Subfield:BadOptions', ...
        'AppeloLogKH must be five finite values ordered [Na K Ca Mg Al].');
end
opt.AppeloLogKH = reshape(opt.AppeloLogKH,1,[]);
if ~(islogical(opt.SaveResults) && isscalar(opt.SaveResults))
    error('Run_One_Subfield:BadOptions','SaveResults must be scalar logical.');
end
if ~(ischar(opt.FieldShapefile) || ...
        (isstring(opt.FieldShapefile) && isscalar(opt.FieldShapefile)))
    error('Run_One_Subfield:BadOptions', ...
        'FieldShapefile must be a character vector or scalar string.');
end
opt.FieldShapefile = char(opt.FieldShapefile);
opt.CellAgronomicFile = validate_optional_file_name( ...
    opt.CellAgronomicFile,'CellAgronomicFile');
opt.CellSpinupStateFile = validate_optional_file_name( ...
    opt.CellSpinupStateFile,'CellSpinupStateFile');
if xor(isempty(opt.CellAgronomicFile),isempty(opt.CellSpinupStateFile))
    error('Run_One_Subfield:IncompleteCellInputs', ...
        'CellAgronomicFile and CellSpinupStateFile must be supplied together.');
end
has_field_shape = ~isempty(opt.FieldShapefile);
if ~has_field_shape && (~isempty(opt.GridCellSizeM) || ...
        ~isempty(opt.FieldFeatureIndices) || ~isempty(opt.FieldPolygonIndex))
    error('Run_One_Subfield:IncompleteFieldGrid', ...
        ['GridCellSizeM, FieldFeatureIndices and FieldPolygonIndex ' ...
         'require FieldShapefile.']);
end
if has_field_shape
    if ~isfile(opt.FieldShapefile)
        error('Run_One_Subfield:MissingShapefile', ...
            'FieldShapefile does not exist: %s',opt.FieldShapefile);
    end
    [~,shape_attributes] = fileattrib(opt.FieldShapefile);
    opt.FieldShapefile = shape_attributes.Name;
    if ~isempty(opt.GridNx) || ~isempty(opt.GridNy)
        error('Run_One_Subfield:ConflictingGridOptions', ...
            'GridNx and GridNy cannot be combined with FieldShapefile.');
    end
end
start_datenum = [];
if ~isempty(opt.StartDate)
    try
        start_datenum = datenum(opt.StartDate);
    catch
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
    if ~(isnumeric(start_datenum) && isscalar(start_datenum) && isfinite(start_datenum))
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
end

RM = Run_Matrix_12();
if ~(isnumeric(run_id) && isscalar(run_id) && isfinite(run_id) && ...
        run_id == fix(run_id) && run_id >= 1 && run_id <= numel(RM))
    error('Run_One_Subfield:BadId','run_id must be 1..%d', numel(RM));
end

R = RM(run_id);
if has_field_shape && ~strcmp(R.layout,'whole_field')
    error('Run_One_Subfield:BadFieldLayout', ...
        'FieldShapefile is supported for whole_field runs 1 to 3 only.');
end

%% This function lives in v2/run; the repository root is one level up.
script_directory = fileparts(mfilename('fullpath'));
if isempty(script_directory); script_directory = pwd; end
current_directory = fileparts(script_directory);
Directory           = fullfile(current_directory,'T&C_Code');
Input_Directory     = fullfile(Directory,'Inputs');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory, Input_Directory, Parameter_Directory);

field_setups = [];
if has_field_shape
    field_setups = Build_Shapefile_Grid_Setups( ...
        opt.FieldShapefile,opt.GridCellSizeM,opt.FieldFeatureIndices);
    if isempty(opt.FieldPolygonIndex)
        if numel(field_setups) ~= 1
            error('Run_One_Subfield:DisconnectedField', ...
                ['FieldShapefile contains %d disconnected polygons. ' ...
                 'Run each with FieldPolygonIndex or use Run_All_12_Parallel.'], ...
                numel(field_setups));
        end
        opt.FieldPolygonIndex = 1;
    elseif opt.FieldPolygonIndex > numel(field_setups)
        error('Run_One_Subfield:BadFieldPolygonIndex', ...
            'FieldPolygonIndex must be between 1 and %d.',numel(field_setups));
    end
end

CONTROL_CASE = 5;
NN_case_study_ERW_MOD = CONTROL_CASE;
OPT_ERW = 1;                    %% chemistry on; apply_rock places the rock
OPT_SoilBiogeochemistry = 1;

fprintf('\n=================================================\n');
fprintf(' Run %d: %s\n', R.id, R.name);
fprintf('   treatment %s (workbook case %d)\n', R.treatment, R.case_number);
fprintf('   layout    %s\n', R.layout);
fprintf('=================================================\n');

%% ---------------------------------------------------------------------
%% Workbook
%% ---------------------------------------------------------------------
wb = fullfile(Input_Directory,'ERW_MOD_PARAM_parameters.xlsx');
MOD_opts = detectImportOptions(wb,"Sheet",'MOD_PARAM');
MOD_opts.DataRange = "A2:Y6";
MOD_opts.VariableTypes = ["double","string","string","double","double", ...
    "string","double","double","string","double", ...
    "double","double","double","double","string", ...
    "double","double","string","double","double", ...
    "double","double","double","string","string"];
MOD_PARAM_input = readtable(wb,MOD_opts);

ERW_opts = detectImportOptions(wb,"Sheet",'ERW_MOD');
ERW_opts.DataRange = "A2:AI6";
ERW_opts.VariableTypes = {'double','string','string','double','double',...
    'string','double','double','string','double',...
    'string','string','string','string','double', ...
    'double','double','double','double','double',...
    'double','double','double','string','double',...
    'string','string','double','double','double',...
    'double','double','double','double','double'};
ERW_opts = setvaropts(ERW_opts,{'f_CEC_in','rock_filepath','mineral','rock_f_in'},'FillValue','NaN');
ERW_PARAM_input = readtable(wb,ERW_opts);

[MOD_PARAM_input,ERW_PARAM_input,ERW_PARAM_control] = ...
    Select_Subfield_Parameter_Rows( ...
    MOD_PARAM_input,ERW_PARAM_input,R.case_number);
MOD_PARAM_input.Zbio(:) = R.biogeo_depth_mm;
id_location = char(MOD_PARAM_input.Location);
fprintf('Active biogeochemical depth: %.0f mm\n',R.biogeo_depth_mm);

%% Spin-up pool state
if isempty(opt.CellSpinupStateFile)
    load(fullfile(Input_Directory,['Data_' id_location '_BG.mat']));  %#ok<LOAD>
end

%% ---------------------------------------------------------------------
%% Forcing
%% ---------------------------------------------------------------------
NN  = MOD_PARAM_input.NumberOfDays*24;
dt  = 3600;
dth = 1;
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'));
ms = length(Zs)-1;
cc = MOD_PARAM_input.cc;

load(fullfile(Input_Directory,['Data_' id_location '.mat']));     %#ok<LOAD>

x1 = 61351;                     %% 01-Jan-2016
if ~isempty(start_datenum)
    [date_error,x1] = min(abs(Date-start_datenum));
    if isempty(x1) || date_error > 1/48
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate is outside the available hourly forcing record.');
    end
end
NN = 1826*24;
if ~isempty(opt.NDays)
    NN = opt.NDays*24;
end
x2 = x1+NN-1;
if x2 > numel(Date)
    error('Run_One_Subfield:BadWindow', ...
        'Requested StartDate and NDays exceed the forcing record.');
end

Date=Date(x1:x2); Pr=Pr(x1:x2); Ta=Ta(x1:x2); Ws=Ws(x1:x2); ea=ea(x1:x2);
SAD1=SAD1(x1:x2); SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); SAB2=SAB2(x1:x2);
Pre=Pre(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2); esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD=PARD(x1:x2);

%% Kantola's first amendment event is 10-Nov-2016. A shortened smoke
%% window ending before that date contains no rock application, so those
%% cells must initialise from the control row. The full run retains the
%% treatment row because its future application events are in the window.
treatment_in_window = strcmp(R.treatment,'control') || ...
    Date(end) >= datenum(2016,11,10);

t_bef = -0.67; t_aft = 1.67;
Ds = esat-ea; Ds(Ds<0) = 0;
Ca_file = load(fullfile(Input_Directory,'Ca_Data.mat'));
d1 = find(abs(Ca_file.Date_CO2-Date(1))<1/36);
d2 = find(abs(Ca_file.Date_CO2-Date(end))<1/36);
Ca = Ca_file.Ca(d1:d2);
Oa = 210000;
Ws(Ws<=0) = 0.01;
%% datevec returns outputs shaped like its input, so a row-vector Date
%% would make [YE MO DA HO] a 1 x 4N row instead of an N x 4 matrix.
%% Force columns, as the single-column frame does.
[YE,MO,DA,HO] = datevec(Date);
Datam = [YE(:) MO(:) DA(:) HO(:)];
clear YE MO DA HO
NNd = ceil(NN/24)+1;

%% ---------------------------------------------------------------------
%% Plot partition and per-column configuration
%% ---------------------------------------------------------------------
g = R.geom;
nx = g.main_grid_nx; ny = g.main_grid_ny;
if has_field_shape
    topo = field_setups(opt.FieldPolygonIndex);
    clear field_setups
    nx = topo.main_grid_nx;
    ny = topo.main_grid_ny;
elseif R.use_subdomains
    if ~isempty(opt.GridNx); nx = opt.GridNx; end
    if ~isempty(opt.GridNy); ny = opt.GridNy; end
    topo = Build_Cell_Topology(g.field_side_m,g.margin_width_m,nx,ny);
else
    topo = Build_Single_Cell_Topology(g.field_side_m^2);
end

%% Attach field topography from a DEM so the 2-D driver resolves
%% elevation-driven lateral and vertical flow (Attach_Topography_To_Topology).
%% Without FieldDEM the run stays flat, exactly as before.
if ~isempty(opt.FieldDEM) && topo.N_cells > 1
    dem = Load_Field_DEM(char(opt.FieldDEM),topo, ...
        'OutResolutionM',dem_resolution(opt.DEMResolutionM,topo));
    topo = Attach_Topography_To_Topology(topo,dem);
    fprintf('Topography attached from %s: relief %.1f m, mean slope %.1f%%\n', ...
        char(opt.FieldDEM),max(topo.cell_z_m)-min(topo.cell_z_m), ...
        100*mean(topo.cell_slope));
end

CellAgronomicInputs = table();
CellBGStates = struct([]);
ERW_cell_input_dependencies = {};
if ~isempty(opt.CellAgronomicFile)
    CellAgronomicInputs = Resolve_Cell_Agronomic_Inputs(topo, ...
        readtable(opt.CellAgronomicFile,'VariableNamingRule','preserve'));
    CellBGStates = Load_Field_Spinup_States(opt.CellSpinupStateFile,topo, ...
        CellAgronomicInputs.input_identity);
    ERW_cell_input_dependencies = {opt.CellAgronomicFile,opt.CellSpinupStateFile};
end

base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
is_main = strcmp(topo.cell_region,'main');
ColCfg = repmat(struct('PARAM_IC',base_param,'PFT_OVERRIDE','', ...
                       'apply_rock',false,'ERW_PARAM',ERW_PARAM_control, ...
                       'MOD_PARAM',table(),'BG_STATE',struct(), ...
                       'agronomic_identity','', ...
                       'appelo_alpha',opt.AppeloAlpha, ...
                       'appelo_logKH',opt.AppeloLogKH), ...
                topo.N_cells,1);
for c = 1:topo.N_cells
    if is_main(c)
        ColCfg(c).apply_rock = R.apply_main;
    else
        ColCfg(c).apply_rock = R.apply_margin;
        ColCfg(c).PFT_OVERRIDE = R.margin_pft;   %% '' keeps the crop
    end
    if ColCfg(c).apply_rock && treatment_in_window
        ColCfg(c).ERW_PARAM = ERW_PARAM_input;
    else
        ColCfg(c).ERW_PARAM = ERW_PARAM_control;
    end
end
if ~isempty(CellAgronomicInputs)
    ERW_cell_rows = cell(topo.N_cells,1);
    for c = 1:topo.N_cells
        ERW_cell_rows{c} = ColCfg(c).ERW_PARAM;
    end
    CellInputs = Configure_Cell_Agronomic_Inputs(MOD_PARAM_input, ...
        vertcat(ERW_cell_rows{:}),CellAgronomicInputs,CellBGStates);
    for c = 1:topo.N_cells
        ColCfg(c).MOD_PARAM = CellInputs(c).MOD_PARAM;
        ColCfg(c).ERW_PARAM = CellInputs(c).ERW_PARAM;
        ColCfg(c).BG_STATE = CellInputs(c).BG_STATE;
        ColCfg(c).agronomic_identity = CellInputs(c).agronomic_identity;
    end
    fprintf('Cell soil inputs: %d validated agronomic rows and spin-up states\n', ...
        topo.N_cells);
end
if has_field_shape
    fprintf(['Domain: %d full-grid columns, %.0f m2 field inside ' ...
        '%.0f m2 simulated (%.3g m cells, %.2f deg rotation)\n'], ...
        topo.N_cells,topo.field_area_m2,topo.simulated_area_m2, ...
        topo.cell_size_m,topo.grid_rotation_degrees);
    fprintf('Grid area outside field: %.0f m2 (%.1f%%)\n', ...
        topo.grid_calibration.outside_grid_area_m2, ...
        100*topo.grid_calibration.outside_grid_fraction);
    fprintf('Rock on: %d field-grid cell(s)\n',sum([ColCfg.apply_rock]));
elseif R.use_subdomains
    fprintf('Domain: %d columns, %.0f m2 main + %.0f m2 margin (%g m wide)\n', ...
        topo.N_cells,sum(topo.cell_area_m2(is_main)), ...
        sum(topo.cell_area_m2(~is_main)),g.margin_width_m);
    fprintf('Rock on: %d main cell(s), %d margin cell(s)\n', ...
        sum([ColCfg(is_main).apply_rock]),sum([ColCfg(~is_main).apply_rock]));
else
    fprintf('Domain: one %.0f m2 column; vertical transport only\n', ...
        topo.cell_area_m2);
    fprintf('Rock on: %d column(s)\n',sum([ColCfg.apply_rock]));
end
if ~treatment_in_window && ~strcmp(R.treatment,'control')
    fprintf('No amendment event falls inside this shortened window.\n');
end

%% ---------------------------------------------------------------------
%% Run
%% ---------------------------------------------------------------------
%% Rolling checkpoint for this run.  One file per run, overwritten in
%% place, so a long column that aborts late (or a machine that sleeps)
%% resumes instead of replaying from day 1.  MAIN_FRAME refuses to resume
%% a checkpoint whose code fingerprint has moved on unless
%% ERW_ckpt_force_resume is set, so an edit to the chemistry starts clean.
ckpt_dir = char(opt.ResultDirectory);
if isempty(ckpt_dir)
    ckpt_dir = fullfile(current_directory,'SIM','Subfield_12');
end
ckpt_dir = fullfile(ckpt_dir,'checkpoints');
if ~exist(ckpt_dir,'dir'); mkdir(ckpt_dir); end
ckpt_tag = R.name;
if has_field_shape
    ckpt_tag = sprintf('%s__%s',ckpt_tag, ...
        field_grid_tag(opt.FieldShapefile,opt.GridCellSizeM, ...
        opt.FieldPolygonIndex,opt.FieldFeatureIndices));
end
if ~isempty(opt.CellAgronomicFile)
    ckpt_tag = sprintf('%s__%s',ckpt_tag, ...
        Field_Input_Tag(opt.CellAgronomicFile,opt.CellSpinupStateFile));
end
if ~isempty(opt.AppeloAlpha)
    ckpt_tag = sprintf('%s__alpha%s',ckpt_tag,number_tag(opt.AppeloAlpha));
end
if ~isempty(opt.AppeloLogKH)
    logkh_parts = arrayfun(@number_tag,opt.AppeloLogKH,'UniformOutput',false);
    ckpt_tag = sprintf('%s__logkh%s',ckpt_tag,strjoin(logkh_parts,'_'));
end
ckpt_tag = sprintf('%s__%dd__from%s',ckpt_tag,NN/24,datestr(Date(1),'yyyymmdd'));
ERW_ckpt_file = fullfile(ckpt_dir, ...
    sprintf('ckpt_run%02d_%s.mat',R.id,ckpt_tag));
%% Interval comes from the environment so it reaches the worker
%% processes the launcher spawns (a launcher-workspace variable does
%% not).  TC_CKPT_MINUTES = 0 disables checkpointing entirely.
if ~exist('ERW_ckpt_minutes','var') || isempty(ERW_ckpt_minutes)
    ERW_ckpt_minutes = str2double(getenv('TC_CKPT_MINUTES'));
    if ~isfinite(ERW_ckpt_minutes) || ERW_ckpt_minutes < 0
        ERW_ckpt_minutes = 10;
    end
end
if ERW_ckpt_minutes == 0
    ERW_ckpt_file = '';
end
%% Debugging escape hatch, also via the environment so it reaches the
%% worker.  Resumes a checkpoint whose code fingerprint has moved on:
%% legitimate when iterating on a LATE-timestep failure whose earlier
%% trajectory the edit cannot touch, and never for production results,
%% because the run then contains two model versions spliced at the
%% resume point.  MAIN_FRAME warns and names that iteration.
if ~exist('ERW_ckpt_force_resume','var') || isempty(ERW_ckpt_force_resume)
    ERW_ckpt_force_resume = strcmpi(strtrim(getenv('TC_CKPT_FORCE_RESUME')), '1');
end

old = cd(Directory);
cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_2D;
clear cleanup
cd(old);

%% ---------------------------------------------------------------------
%% Save
%% ---------------------------------------------------------------------
result_file = '';
if opt.SaveResults
    fl_save = char(opt.ResultDirectory);
    if isempty(fl_save)
        fl_save = fullfile(current_directory,'SIM','Subfield_12');
    end
    if ~exist(fl_save,'dir'); mkdir(fl_save); end
    tag = R.name;
    if has_field_shape
        tag = sprintf('%s__%s',tag, ...
            field_grid_tag(opt.FieldShapefile,opt.GridCellSizeM, ...
            opt.FieldPolygonIndex,opt.FieldFeatureIndices));
    elseif R.use_subdomains && (~isempty(opt.GridNx) || ~isempty(opt.GridNy))
        tag = sprintf('%s__grid%dx%d',tag,nx,ny);
    end
    if ~isempty(opt.CellAgronomicFile)
        tag = sprintf('%s__%s',tag, ...
            Field_Input_Tag(opt.CellAgronomicFile,opt.CellSpinupStateFile));
    end
    if ~isempty(opt.AppeloAlpha)
        tag = sprintf('%s__alpha%s',tag,number_tag(opt.AppeloAlpha));
    end
    if ~isempty(opt.AppeloLogKH)
        logkh_parts = arrayfun(@number_tag,opt.AppeloLogKH,'UniformOutput',false);
        tag = sprintf('%s__logkh%s',tag,strjoin(logkh_parts,'_'));
    end
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    if ~isempty(start_datenum); tag = sprintf('%s__from%s',tag,datestr(start_datenum,'yyyymmdd')); end
    result_file = fullfile(fl_save,sprintf('TC_2D_run%02d_%s.mat',R.id,tag));
    tmp_result = [tempname(fl_save) '.mat'];
    try
        save(tmp_result,'-v7.3');
        movefile(tmp_result,result_file,'f');
    catch ME
        if isfile(tmp_result); delete(tmp_result); end
        rethrow(ME)
    end
    fprintf('Saved: %s\n', result_file);
    %% The run finished and its results are on disk, so the rolling
    %% checkpoint is now dead weight.  Drop it: these files are large
    %% and keeping one per completed run would fill the volume.
    if exist('ERW_ckpt_file','var') && ~isempty(ERW_ckpt_file) && isfile(ERW_ckpt_file)
        delete(ERW_ckpt_file);
    end
end

end

function validate_optional_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 1)
    error('Run_One_Subfield:BadOptions','%s must be a positive integer.',name);
end
end

function file_name = validate_optional_file_name(value,name)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a character vector or scalar string.',name);
end
file_name = char(value);
if isempty(file_name)
    return
end
if ~isfile(file_name)
    error('Run_One_Subfield:MissingCellInput', ...
        '%s does not exist: %s',name,file_name);
end
[~,attributes] = fileattrib(file_name);
file_name = attributes.Name;
end

function validate_optional_nonnegative_scalar(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a finite nonnegative scalar.',name);
end
end

function validate_optional_positive_scalar(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a finite positive scalar.',name);
end
end

function tag = field_grid_tag( ...
        shape_file,cell_size_m,polygon_index,feature_indices)
[~,shape_name] = fileparts(shape_file);
shape_name = regexprep(shape_name,'[^A-Za-z0-9]+','_');
if isempty(cell_size_m)
    resolution = 'geometry';
else
    resolution = sprintf('cell%sm',number_tag(cell_size_m));
end
tag = sprintf('field_%s__%s',shape_name,resolution);
tag = sprintf('%s__polygon%02d',tag,polygon_index);
if ~isempty(feature_indices)
    parts = arrayfun(@(x) sprintf('%d',x),feature_indices, ...
        'UniformOutput',false);
    tag = sprintf('%s__features%s',tag,strjoin(parts,'_'));
end
end

function validate_optional_index_vector(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isvector(value) && all(isfinite(value)) && ...
        all(value == fix(value)) && all(value >= 1) && ...
        numel(unique(value)) == numel(value))
    error('Run_One_Subfield:BadOptions', ...
        '%s must contain unique positive integers.',name);
end
end

function tag = number_tag(value)
tag = strrep(sprintf('%.17g',value),'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
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
