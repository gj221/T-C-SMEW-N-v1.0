function [reports,files] = Krige_Cell_Agronomic_Values( ...
    shape_file, sample_file, attribute_map, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Krige_Cell_Agronomic_Values
%%
%% Krige measured soil attributes onto every simulation cell centre of a
%% field's calibrated grid, and write a per-cell REVIEW table (prediction,
%% uncertainty and extrapolation diagnostics) -- not a merged parameter
%% file. Apply_Kriged_Cell_Values does the deliberate, threshold-gated
%% merge into the Cell_Agronomic_*.csv template that Run_Field_Cell_Spinup
%% actually reads; keeping the two steps separate means an estimate you
%% have not reviewed can never silently become a model input.
%%
%% INPUTS
%%   shape_file     field-boundary shapefile (e.g. the control fields)
%%   sample_file    soil-sample workbook (xlsx) or delimited file
%%   attribute_map  struct array, one entry per attribute to krige:
%%                    .column        workbook column name, e.g. 'pH'
%%                    .template      matching Cell_Agronomic_*.csv column
%%                                   stem, e.g. 'bulk_pH' (the template
%%                                   carries baseline_<template> and
%%                                   measured_<template>)
%%                    .scale         multiplier applied after reading,
%%                                   e.g. 0.01 for Organic Matter [%] ->
%%                                   OM_fraction [-].  Default 1.
%%                  Bruce Farm control-field default (pH, Organic Matter,
%%                  Cation Exchange Capacity) is built if attribute_map is
%%                  left [].  Soil TEXTURE (sand/clay) is not in this
%%                  workbook and is out of scope here.
%%
%% Name/value
%%   'FeatureIndices'  polygon features (default [1 2], the control fields)
%%   'CellSizeM'        process-grid width (default: calibrated coverage
%%                      width per polygon, matching Write_Field_Agronomic_
%%                      Templates / Run_Field_Cell_Spinup)
%%   'SampleSheet'      workbook sheet (default 2, per Run_Bruce_Interpolation)
%%   'Stage'             'Sample Type' value to use (default 'BL')
%%   'XColumn','YColumn' sample coordinate columns (default 'lng (X)'/'lat (Y)')
%%   'VariogramModel'    default 'exponential'
%%   'OutputDirectory'   where the review CSVs go (default v2/figures)
%%
%% OUTPUTS
%%   reports  struct array, one per polygon: .table (the review table),
%%            .polygon_index, .setup
%%   files    the written review CSV paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'));

opt = struct('FeatureIndices',[1 2],'CellSizeM',[],'SampleSheet',2, ...
    'Stage','BL','XColumn','lng (X)','YColumn','lat (Y)', ...
    'VariogramModel','exponential','OutputDirectory',fullfile(v2,'figures'));
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if ~exist(opt.OutputDirectory,'dir'); mkdir(opt.OutputDirectory); end

if nargin < 3 || isempty(attribute_map)
    attribute_map = struct( ...
        'column',   {'pH','Organic Matter','Cation Exchange Capacity'}, ...
        'template', {'bulk_pH','OM_fraction','CEC_cmolc_kg'}, ...
        'scale',    {1, 0.01, 1});
end
for a = 1:numel(attribute_map)
    if ~isfield(attribute_map(a),'scale') || isempty(attribute_map(a).scale)
        attribute_map(a).scale = 1;
    end
end

points = Read_Sample_Points(sample_file,opt.SampleSheet,opt.XColumn,opt.YColumn);
stage_col = 'Sample Type';
if ~ismember(stage_col,points.Properties.VariableNames)
    error('Krige_Cell_Agronomic_Values:NoStageColumn', ...
        'Sample table has no "%s" column.',stage_col);
end
points = points(strcmp(string(points.(stage_col)),opt.Stage),:);
if isempty(points)
    error('Krige_Cell_Agronomic_Values:NoStageRows', ...
        'No sample carries Sample Type = %s.',opt.Stage);
end

setups = Build_Shapefile_Grid_Setups(shape_file,opt.CellSizeM,opt.FeatureIndices);
[~,shape_name] = fileparts(shape_file);
reports = repmat(struct('table',table(),'polygon_index',0,'setup',struct()), ...
    numel(setups),1);
files = cell(numel(setups),1);

for p = 1:numel(setups)
    topo = setups(p);
    T = table(topo.cell_id,repmat(p,topo.N_cells,1), ...
        'VariableNames',{'cell_id','polygon_index'});
    for a = 1:numel(attribute_map)
        spec = attribute_map(a);
        col = spec.column;
        [pred,sd,nsamp,near_m,hull_ok,n_used] = krige_one_attribute( ...
            points,topo,col,spec.scale,opt);
        T.(sprintf('%s_predicted',spec.template)) = pred;
        T.(sprintf('%s_std',spec.template)) = sd;
        T.(sprintf('%s_nearest_sample_m',spec.template)) = near_m;
        T.(sprintf('%s_in_sample_hull',spec.template)) = hull_ok;
        T.(sprintf('%s_n_samples',spec.template)) = repmat(n_used,topo.N_cells,1);
        fprintf(['Krige_Cell_Agronomic_Values: polygon %d, %s <- "%s" ' ...
            '(n=%d, %d/%d cells inside the sample hull)\n'], ...
            p,spec.template,col,nsamp,nnz(hull_ok),topo.N_cells);
    end
    reports(p).table = T;
    reports(p).polygon_index = p;
    reports(p).setup = topo;
    files{p} = fullfile(opt.OutputDirectory,sprintf( ...
        'Cell_Kriged_Values_%s_polygon_%02d.csv',shape_name,p));
    writetable(T,files{p});
    fprintf('  wrote %s\n',files{p});
end
end

%% =======================================================================
function [pred,sd,nsamp,near_m,hull_ok,n_used] = krige_one_attribute( ...
    points,topo,column,scale,opt)
if ~ismember(column,points.Properties.VariableNames)
    error('Krige_Cell_Agronomic_Values:MissingColumn', ...
        'Sample table has no "%s" column.',column);
end
values = Numeric_Sample_Column(points,column);
located = isfinite(values);
[samples,~] = Project_Sample_Points( ...
    points(located,:),topo,opt.XColumn,opt.YColumn);
v = Numeric_Sample_Column(samples,column);
finite = isfinite(v);
samples = samples(finite,:); v = v(finite);
n_used = height(samples);
if n_used < 3
    warning('Krige_Cell_Agronomic_Values:TooFewSamples', ...
        ['Only %d located finite sample(s) for "%s" in this polygon; ' ...
         'leaving it unkriged (NaN).'],n_used,column);
    pred = nan(topo.N_cells,1); sd = pred; near_m = pred;
    hull_ok = false(topo.N_cells,1); nsamp = n_used;
    return
end
model = Fit_Ordinary_Kriging(samples.x_grid_m,samples.y_grid_m,v, ...
    opt.VariogramModel);
[pred,var_,near_m,hull_ok] = Krige_Attribute_At_Points( ...
    model,topo.cell_xc,topo.cell_yc);
pred = pred*scale;
sd = sqrt(max(var_,0))*abs(scale);
nsamp = n_used;
end
