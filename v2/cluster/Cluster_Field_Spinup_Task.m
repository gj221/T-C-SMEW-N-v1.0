function out = Cluster_Field_Spinup_Task(task, cell_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster_Field_Spinup_Task
%%
%% The one MATLAB entry point the SLURM scripts call.  Everything about
%% the field comes from environment variables (set by submit_field_spinup.sh
%% or by hand), so the same batch scripts serve every polygon and grid.
%%
%%   task        'plan'      write the cell plan (indices, subsample) and
%%                           the array-index lists the submitter reads
%%               'cell'      spin up ONE cell (cell_index = SLURM array id):
%%                           a full spin-up for Method per_cell or for a
%%                           subsample cell, the RF-seeded short spin-up
%%                           otherwise
%%               'assemble'  gather every cell into CellBGStates.mat
%%               'validate'  Run_RF_Validation_Experiment on the
%%                           assembled truth states
%%   cell_index  integer, 'cell' only (a string from the shell is accepted)
%%
%% Environment (TC_ prefix; paths POSIX on the cluster)
%%   TC_SHAPEFILE      field-boundary shapefile
%%   TC_AGRONOMIC_FILE merged per-cell agronomic CSV (Apply_Kriged_Cell_Values)
%%   TC_OUTPUT_DIR     where cell files, the plan and CellBGStates.mat go
%%   TC_METHOD         per_cell | rf                        (default per_cell)
%%   TC_FEATURES       polygon features, e.g. "1 2 3"        (default all)
%%   TC_POLYGON        polygon index                         (default 1)
%%   TC_CELL_M         process-grid width [m]                (default: calibrated)
%%   TC_MIN_COVERAGE   coverage rule, e.g. 0.5 (diagnostic)  (default 0.95)
%%   TC_DEM            DEM directory / GeoTIFF (required for rf, validate)
%%   TC_SUBSAMPLE      RF subsample fraction                 (default 0.12)
%%   TC_CASE           workbook case number                  (default 5)
%%   TC_FORCING_FILE, TC_METEO_FILE, TC_METEO_START, TC_METEO_END,
%%   TC_TRANSIENT_START, TC_TRANSIENT_END
%%                     the spin-up columns a workbook row may lack
%%                     (Run_Field_Cell_Spinup 'MODParamColumns'); unset ->
%%                     the row must carry them
%%   TC_VALIDATION_FRACTIONS  e.g. "0.10 0.15 0.20 0.30"      (validate)
%%   TC_VALIDATION_SEEDS      e.g. "1 2 3 4 5"                (validate)
%%   TC_RF_STATES_DIR         cell files of an RF-path run to score
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
method  = getenv_default('TC_METHOD','per_cell');
feats   = str2num(getenv_default('TC_FEATURES','')); %#ok<ST2NM>
poly    = str2double(getenv_default('TC_POLYGON','1'));
cell_m  = str2double(getenv_default('TC_CELL_M','NaN'));
cov     = str2double(getenv_default('TC_MIN_COVERAGE','NaN'));
dem     = getenv_default('TC_DEM','');
frac    = str2double(getenv_default('TC_SUBSAMPLE','0.12'));
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
    'SubsampleFraction',frac,'MODParamColumns',cols};
if isfinite(cell_m); common = [common, {'GridCellSizeM',cell_m}]; end
if isfinite(cov);    common = [common, {'MinCoverage',cov}]; end
if ~isempty(dem);    common = [common, {'FieldDEM',dem}]; end

fprintf('Cluster_Field_Spinup_Task %s: %s, polygon %d, method %s, output %s\n', ...
    task, shp, poly, method, outdir);
switch lower(task)
    case 'plan'
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'PlanOnly',true);
    case 'cell'
        if nargin < 2; cell_index = getenv('SLURM_ARRAY_TASK_ID'); end
        if ischar(cell_index) || isstring(cell_index); cell_index = str2double(cell_index); end
        if ~(isfinite(cell_index) && cell_index >= 1)
            error('Cluster_Field_Spinup_Task:NoCellIndex', ...
                'task ''cell'' needs a cell index (argument or SLURM_ARRAY_TASK_ID).');
        end
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'CellIndex',cell_index);
    case 'assemble'
        out = Run_Field_Cell_Spinup(shp,agro,common{:},'AssembleOnly',true);
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
            'task must be plan | cell | assemble | validate (got %s).', task);
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
