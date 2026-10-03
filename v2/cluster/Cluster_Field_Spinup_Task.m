function out = Cluster_Field_Spinup_Task(task, index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster_Field_Spinup_Task
%%
%% The one MATLAB entry point the SLURM scripts call.  Everything about
%% the field comes from environment variables (set by submit_field_spinup.sh
%% or by hand), so the same batch scripts serve every polygon and grid.
%%
%%   task   'plan'         write the plan and the array-index lists the
%%                         submitter reads
%%          'class_seed'   (flux_tracking) the seed of class <index>
%%          'track'        (flux_tracking) the 2-D tracking run; finished
%%                         segments are reused, so a job that runs out of
%%                         time is simply submitted again
%%          'replay'       (flux_tracking) tracked cell <index>'s replay
%%          'reconstruct'  (flux_tracking) cell <index> from the forest
%%          'cell'         (per_cell) cell <index>'s own spin-up
%%          'assemble'     gather every cell into CellBGStates.mat
%%          'validate'     Run_RF_Validation_Experiment on the assembled
%%                         truth states
%%   index  integer for the array tasks (default SLURM_ARRAY_TASK_ID; a
%%          string from the shell is accepted)
%%
%% Environment (TC_ prefix; paths POSIX on the cluster)
%%   TC_SHAPEFILE      field-boundary shapefile
%%   TC_AGRONOMIC_FILE merged per-cell agronomic CSV (Apply_Kriged_Cell_Values)
%%   TC_OUTPUT_DIR     where cell files, the plan and CellBGStates.mat go
%%   TC_METHOD         flux_tracking | per_cell              (default flux_tracking)
%%   TC_FEATURES       polygon features, e.g. "1 2 3"        (default all)
%%   TC_POLYGON        polygon index                         (default 1)
%%   TC_CELL_M         process-grid width [m]                (default: calibrated)
%%   TC_MIN_COVERAGE   coverage rule, e.g. 0.5 (diagnostic)  (default 0.95)
%%   TC_DEM            DEM directory / GeoTIFF (flux_tracking, validate)
%%   TC_SUBSAMPLE      tracked fraction                      (default 0.10)
%%   TC_CLASS_COLUMNS  agronomic columns defining the land cover x soil
%%                     classes, e.g. "land_cover soil_class" (default: those
%%                     of land_cover, vegetation, pft, soil_class,
%%                     texture_class the table carries)
%%   TC_CASE           workbook case number                  (default 5)
%%   TC_FORCING_FILE, TC_METEO_FILE, TC_METEO_START, TC_METEO_END,
%%   TC_TRANSIENT_START, TC_TRANSIENT_END
%%                     the spin-up columns a workbook row may lack
%%                     (Run_Field_Cell_Spinup 'MODParamColumns'); unset ->
%%                     the row must carry them
%%   TC_VALIDATION_FRACTIONS  e.g. "0.10 0.15 0.20 0.30"      (validate)
%%   TC_VALIDATION_SEEDS      e.g. "1 2 3 4 5"                (validate)
%%   TC_RF_STATES_DIR         cell files of a flux_tracking run to score
%%                            against the truth (validate, optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
addpath(fullfile(v2,'spinup'), fullfile(v2,'T&C_Code'), ...
    fullfile(v2,'T&C_Code','Parameters'), fullfile(v2,'T&C_Code','Inputs'), ...
    fullfile(v2,'interpolation'), fullfile(v2,'run'));
setenv('TC_CKPT_MINUTES','0');

shp     = need('TC_SHAPEFILE');
agro    = need('TC_AGRONOMIC_FILE');
outdir  = need('TC_OUTPUT_DIR');
method  = getenv_default('TC_METHOD','flux_tracking');
feats   = str2num(getenv_default('TC_FEATURES','')); %#ok<ST2NM>
poly    = str2double(getenv_default('TC_POLYGON','1'));
cell_m  = str2double(getenv_default('TC_CELL_M','NaN'));
cov     = str2double(getenv_default('TC_MIN_COVERAGE','NaN'));
dem     = getenv_default('TC_DEM','');
frac    = str2double(getenv_default('TC_SUBSAMPLE','0.10'));
classes = strsplit(strtrim(getenv_default('TC_CLASS_COLUMNS','')));
classes = classes(~cellfun(@isempty,classes));
casenum = str2double(getenv_default('TC_CASE','5'));
if ~exist(outdir,'dir'); mkdir(outdir); end
state_file = fullfile(outdir,'CellBGStates.mat');

cols = struct();
map = {'TC_FORCING_FILE','forcing_file'; 'TC_METEO_FILE','meteo_file'; ...
    'TC_METEO_START','meteo_start'; 'TC_METEO_END','meteo_end'; ...
    'TC_TRANSIENT_START','transient_start'; 'TC_TRANSIENT_END','transient_end'};
for k = 1:size(map,1)
    v = getenv(map{k,1});
    if ~isempty(v); cols.(map{k,2}) = v; end
end
if isempty(fieldnames(cols)); cols = []; end

common = {'FieldFeatureIndices',feats,'FieldPolygonIndex',poly, ...
    'CaseNumber',casenum,'OutputFile',state_file,'Method',method, ...
    'SubsampleFraction',frac,'MODParamColumns',cols,'ClassColumns',classes};
if isfinite(cell_m); common = [common, {'GridCellSizeM',cell_m}]; end
if isfinite(cov);    common = [common, {'MinCoverage',cov}]; end
if ~isempty(dem);    common = [common, {'FieldDEM',dem}]; end

fprintf('Cluster_Field_Spinup_Task %s: %s, polygon %d, method %s, output %s\n', ...
    task, shp, poly, method, outdir);
switch lower(task)
    case {'plan','track','assemble'}
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'Stage',lower(task));
    case 'class_seed'
        if nargin < 2; index = []; end
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'Stage','class_seed', ...
            'ClassIndex',array_index(index,task));
    case {'replay','reconstruct','cell'}
        if nargin < 2; index = []; end
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'Stage',lower(task), ...
            'CellIndex',array_index(index,task));
    case 'validate'
        fr = str2num(getenv_default('TC_VALIDATION_FRACTIONS','0.10 0.15 0.20 0.30')); %#ok<ST2NM>
        sd = str2num(getenv_default('TC_VALIDATION_SEEDS','1 2 3 4 5')); %#ok<ST2NM>
        rfdir = getenv_default('TC_RF_STATES_DIR','');
        vopts = {'FieldFeatureIndices',feats,'FieldPolygonIndex',poly, ...
            'CaseNumber',casenum,'TruthStateFile',state_file,'FieldDEM',dem, ...
            'Fractions',fr,'Seeds',sd,'OutputDirectory',outdir,'MODParamColumns',cols};
        if isfinite(cell_m); vopts = [vopts, {'GridCellSizeM',cell_m}]; end
        if isfinite(cov);    vopts = [vopts, {'MinCoverage',cov}]; end
        if ~isempty(rfdir);  vopts = [vopts, {'RFStatesDirectory',rfdir}]; end
        out = Run_RF_Validation_Experiment(shp,agro,vopts{:});
    otherwise
        error('Cluster_Field_Spinup_Task:BadTask', ...
            ['task must be plan | class_seed | track | replay | reconstruct | ' ...
             'cell | assemble | validate (got %s).'], task);
end
end

function index = array_index(index,task)
%% The array task's index: the argument, else SLURM_ARRAY_TASK_ID.
if isempty(index); index = getenv('SLURM_ARRAY_TASK_ID'); end
if ischar(index) || isstring(index); index = str2double(index); end
if ~(isnumeric(index) && isscalar(index) && isfinite(index) && index >= 1)
    error('Cluster_Field_Spinup_Task:NoCellIndex', ...
        'task ''%s'' needs an index (argument or SLURM_ARRAY_TASK_ID).',task);
end
end

function v = need(name)
v = getenv(name);
if isempty(v)
    error('Cluster_Field_Spinup_Task:MissingEnv','%s must be set.', name);
end
end

function v = getenv_default(name, default)
v = getenv(name);
if isempty(v); v = default; end
end
