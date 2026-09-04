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
%% UNCERTAINTY.  A control field's ten-odd samples cannot resolve a
%% semivariogram (Fit_Empirical_Variogram_Range .resolved = false), so the
%% variogram range is not fitted here -- it is treated as unknown within
%% the interval the sample geometry supports (Kriging_Range_Prior_Bounds)
%% and integrated out by Monte Carlo (Krige_Attribute_Range_Ensemble).
%% <stem>_std is therefore the TOTAL standard deviation, interpolation
%% error plus range uncertainty by the law of total variance, and it is
%% what Apply_Kriged_Cell_Values gates on; the two components are written
%% alongside it.  Each attribute is also leave-one-out cross-validated
%% (Krige_Attribute_Cross_Validate) so the review carries an empirical
%% calibration check (rms standardised error ~1) rather than only the
%% model's own opinion of itself.
%%
%% Name/value
%%   'FeatureIndices'  polygon features (default []: every feature in the
%%                      shapefile -- AdamWilbourne_Control_Fields.shp has
%%                      three, Control 2 / Control 1 / Control 3)
%%   'CellSizeM'        process-grid width (default: calibrated coverage
%%                      width per polygon, matching Write_Field_Agronomic_
%%                      Templates / Run_Field_Cell_Spinup)
%%   'SampleSheet'      workbook sheet (default 2, per Run_Bruce_Interpolation)
%%   'Stage'             'Sample Type' value to use (default 'BL')
%%   'XColumn','YColumn' sample coordinate columns (default: auto --
%%                       'lng (X)'/'lat (Y)' or 'sme_lng'/'sme_lat')
%%   'VariogramModel'    default 'exponential'
%%   'RangeDraws'        Monte Carlo range draws per attribute (default
%%                       200).  0 = classical single fit at Fit_Ordinary_
%%                       Kriging's default 80 m range, reported with
%%                       range_assumed = true.
%%   'RangeBoundsM'      [lower upper] range prior (default: from the
%%                       samples, per polygon and attribute)
%%   'Seed'              draw seed (default 0, reproducible review tables)
%%   'CrossValidate'     leave-one-out check per attribute (default true;
%%                       needs >= 4 samples)
%%   'OutputDirectory'   where the review CSVs go (default v2/figures)
%%
%% OUTPUTS
%%   reports  struct array, one per polygon: .table (the review table),
%%            .polygon_index, .setup, .attributes (per attribute: n_samples,
%%            range_bounds_m, n_draws, range_assumed, sill, nugget, loocv)
%%   files    the written review CSV paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'));

opt = struct('FeatureIndices',[],'CellSizeM',[],'SampleSheet',2, ...
    'Stage','BL','XColumn','','YColumn','', ...
    'VariogramModel','exponential','RangeDraws',200,'RangeBoundsM',[], ...
    'Seed',0,'CrossValidate',true,'OutputDirectory',fullfile(v2,'figures'));
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if ~exist(opt.OutputDirectory,'dir'); mkdir(opt.OutputDirectory); end

%% Coordinate columns: the agronomic workbook carries 'lng (X)'/'lat (Y)',
%% the complete sample-catalogue export carries 'sme_lng'/'sme_lat'.
%% Resolve whichever is present unless the caller named the columns.
[points,opt.XColumn,opt.YColumn] = read_points_any_frame( ...
    sample_file,opt.SampleSheet,opt.XColumn,opt.YColumn);

if nargin < 3 || isempty(attribute_map)
    attribute_map = default_attribute_map(points);
end
for a = 1:numel(attribute_map)
    if ~isfield(attribute_map(a),'scale') || isempty(attribute_map(a).scale)
        attribute_map(a).scale = 1;
    end
end
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
reports = repmat(struct('table',table(),'polygon_index',0,'setup',struct(), ...
    'attributes',struct()),numel(setups),1);
files = cell(numel(setups),1);

for p = 1:numel(setups)
    topo = setups(p);
    T = table(topo.cell_id,repmat(p,topo.N_cells,1), ...
        'VariableNames',{'cell_id','polygon_index'});
    %% Projected centres let Apply_Kriged_Cell_Values verify the merge
    %% geometrically (cell_id alone is only unique within one grid).
    if isfield(topo,'cell_xc_projected_m') && isfield(topo,'cell_yc_projected_m')
        T.cell_x_projected_m = topo.cell_xc_projected_m(:);
        T.cell_y_projected_m = topo.cell_yc_projected_m(:);
    end
    attr_info = cell(numel(attribute_map),1);
    for a = 1:numel(attribute_map)
        spec = attribute_map(a);
        col = spec.column;
        [r,info] = krige_one_attribute(points,topo,col,spec.scale,opt);
        stem = spec.template;
        T.(sprintf('%s_predicted',stem))          = r.prediction;
        T.(sprintf('%s_std',stem))                = r.std;
        T.(sprintf('%s_std_interpolation',stem))  = r.std_interpolation;
        T.(sprintf('%s_std_range',stem))          = r.std_range;
        T.(sprintf('%s_nearest_sample_m',stem))   = r.nearest_sample_m;
        T.(sprintf('%s_in_sample_hull',stem))     = r.in_hull;
        T.(sprintf('%s_n_samples',stem))          = repmat(info.n_samples,topo.N_cells,1);
        T.(sprintf('%s_range_assumed',stem))      = repmat(info.range_assumed,topo.N_cells,1);
        info.template = stem; info.column = col;
        attr_info{a} = info;
        report_attribute(p,stem,col,info,r,topo.N_cells);
    end
    reports(p).table = T;
    reports(p).polygon_index = p;
    reports(p).setup = topo;
    reports(p).attributes = [attr_info{:}];
    files{p} = fullfile(opt.OutputDirectory,sprintf( ...
        'Cell_Kriged_Values_%s_polygon_%02d.csv',shape_name,p));
    writetable(T,files{p});
    fprintf('  wrote %s\n',files{p});
end
end

%% =======================================================================
function [points,xcol,ycol] = read_points_any_frame(sample_file,sheet,xcol,ycol)
if ~isempty(xcol) || ~isempty(ycol)
    points = Read_Sample_Points(sample_file,sheet,xcol,ycol);
    return
end
candidates = {'lng (X)','lat (Y)'; 'sme_lng','sme_lat'; 'long','lat'; 'lng','lat'};
last = [];
for k = 1:size(candidates,1)
    try
        points = Read_Sample_Points(sample_file,sheet,candidates{k,1},candidates{k,2});
        xcol = candidates{k,1}; ycol = candidates{k,2};
        return
    catch ME
        %% Read_Sample_Points raises its MissingColumn tag through a
        %% one-argument error() call, so it arrives as the message, not
        %% the identifier -- accept either.
        missing = endsWith(ME.identifier,':MissingColumn') || ...
            contains(ME.message,'MissingColumn');
        if ~missing; rethrow(ME); end
        last = ME;
    end
end
error('Krige_Cell_Agronomic_Values:NoCoordinateColumns', ...
    ['%s has none of the known coordinate column pairs (%s); pass ' ...
     'XColumn/YColumn.  Last error: %s'], char(string(sample_file)), ...
    strjoin(cellfun(@(a,b) sprintf('%s/%s',a,b),candidates(:,1), ...
    candidates(:,2),'UniformOutput',false),', '), last.message);
end

%% =======================================================================
function map = default_attribute_map(points)
%% Bruce Farm default: pH, organic matter, CEC under either naming scheme
%% (agronomic workbook vs. the complete sample-catalogue export).  Soil
%% TEXTURE is added when the file carries it (the catalogue export's
%% 'Sand Mass (%)' / 'Clay Mass (%)'); the workbook does not.
names = points.Properties.VariableNames;
pick = @(alts) alts(ismember(alts,names));
cands = { ...
    pick({'pH','pH Value'}),                                  'bulk_pH',       1; ...
    pick({'Organic Matter','Organic Matter (%)'}),            'OM_fraction',   0.01; ...
    pick({'Cation Exchange Capacity','CEC Value'}),           'CEC_cmolc_kg',  1; ...
    pick({'Sand Mass (%)','Percent Sand','sand_pct'}),        'sand_fraction', 0.01; ...
    pick({'Clay Mass (%)','Percent Clay','clay_pct'}),        'clay_fraction', 0.01};
keep = ~cellfun(@isempty,cands(:,1));
cands = cands(keep,:);
if isempty(cands)
    error('Krige_Cell_Agronomic_Values:NoKnownAttributes', ...
        'No default attribute column found; pass attribute_map explicitly.');
end
map = struct('column',cellfun(@(c) c{1},cands(:,1),'UniformOutput',false).', ...
    'template',cands(:,2).','scale',cands(:,3).');
end

%% =======================================================================
function [r,info] = krige_one_attribute(points,topo,column,scale,opt)
%% r: per-cell prediction / std components / hull diagnostics, already in
%% template units (scale applied).  info: what the estimate rests on.
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
%% Coincident sample locations (replicates) would each zero the kriging
%% right-hand side at that point and collapse the variance there to ~0
%% even though the replicates disagree.  Average them, and carry their
%% scatter forward as the floor of the nugget the estimator reports.
[x0,y0,v0,replicate_var,n_dup] = merge_coincident(samples.x_grid_m, ...
    samples.y_grid_m,v);
if n_dup > 0
    warning('Krige_Cell_Agronomic_Values:CoincidentSamples', ...
        ['"%s": %d sample(s) share a location with another; replicates ' ...
         'averaged, their scatter (variance %.3g) kept as a nugget floor.'], ...
        column,n_dup,replicate_var);
end
samples = table(x0,y0,'VariableNames',{'x_grid_m','y_grid_m'}); v = v0;
n_used = height(samples);

info = struct('n_samples',n_used,'range_bounds_m',[NaN NaN],'n_draws',0, ...
    'range_assumed',true,'sill',NaN,'nugget',NaN,'loocv',[]);
nanv = nan(topo.N_cells,1);
r = struct('prediction',nanv,'std',nanv,'std_interpolation',nanv, ...
    'std_range',nanv,'nearest_sample_m',nanv,'in_hull',false(topo.N_cells,1));
if n_used < 3
    warning('Krige_Cell_Agronomic_Values:TooFewSamples', ...
        ['Only %d located finite sample(s) for "%s" in this polygon; ' ...
         'leaving it unkriged (NaN).'],n_used,column);
    return
end

x = samples.x_grid_m; y = samples.y_grid_m;
vp_floor = struct();
if replicate_var > 0
    %% Nugget floor from the replicate scatter (sill stays the default).
    vp_floor.nugget = max(0.05*var(v,1), replicate_var);
    vp_floor.sill = max(var(v,1), vp_floor.nugget);
end
ens = Krige_Attribute_Range_Ensemble(x,y,v,topo.cell_xc,topo.cell_yc, ...
    'VariogramModel',opt.VariogramModel,'VariogramParameters',vp_floor, ...
    'RangeDraws',opt.RangeDraws,'RangeBoundsM',opt.RangeBoundsM,'Seed',opt.Seed);
s = abs(scale);
r.prediction        = ens.prediction*scale;
r.std               = ens.std*s;
r.std_interpolation = ens.std_interpolation*s;
r.std_range         = ens.std_range*s;
r.nearest_sample_m  = ens.nearest_sample_m;
r.in_hull           = ens.in_hull;

info.range_bounds_m = ens.range_bounds_m;
info.n_draws        = ens.n_draws;
info.range_assumed  = ens.range_assumed;
info.sill           = ens.sill*s^2;
info.nugget         = ens.nugget*s^2;
if opt.CrossValidate && n_used >= 5
    %% Each fold re-derives sill, nugget and range prior from its own
    %% n-1 samples -- the estimator as deployed on that many samples --
    %% so the held-out value can never leak into the fold's variance.
    %% Only an assumed range (RangeDraws = 0) is shared, being external.
    cv_vp = vp_floor;
    if ens.range_assumed; cv_vp.range = ens.range_draws_m(1); end
    cv = Krige_Attribute_Cross_Validate(x,y,v, ...
        'VariogramModel',opt.VariogramModel,'VariogramParameters',cv_vp, ...
        'RangeDraws',opt.RangeDraws,'RangeBoundsM',opt.RangeBoundsM, ...
        'Seed',opt.Seed);
    cv.rmse = cv.rmse*s; cv.mean_error = cv.mean_error*s;
    cv.observed = cv.observed*scale; cv.predicted = cv.predicted*scale;
    cv.error = cv.error*scale; cv.std = cv.std*s;
    cv.std_interpolation = cv.std_interpolation*s;
    info.loocv = cv;
end
end

%% =======================================================================
function [x,y,v,replicate_var,n_dup] = merge_coincident(x,y,v)
%% Average samples that share a location (to 1e-6 m); return the pooled
%% within-location variance of the replicates as a nugget floor.
x = x(:); y = y(:); v = v(:);
replicate_var = 0; n_dup = 0;
if isempty(v); return; end            % a polygon with no located samples
key = round([x y]*1e6);
[~,~,g] = unique(key,'rows','stable');
n_groups = max(g);
xm = accumarray(g,x,[n_groups 1],@mean);
ym = accumarray(g,y,[n_groups 1],@mean);
vm = accumarray(g,v,[n_groups 1],@mean);
counts = accumarray(g,1,[n_groups 1]);
n_dup = sum(counts(counts > 1));
if n_dup > 0
    dev = v - vm(g);
    replicate_var = sum(dev.^2)/max(1,sum(counts(counts > 1)-1));
end
x = xm; y = ym; v = vm;
end

%% =======================================================================
function report_attribute(p,stem,column,info,r,n_cells)
if info.n_samples < 3; return; end
if info.range_assumed
    range_txt = sprintf('range %.0f m ASSUMED (not fitted, %d draws)', ...
        info.range_bounds_m(1),info.n_draws);
else
    range_txt = sprintf('range prior %.0f-%.0f m, %d draws', ...
        info.range_bounds_m(1),info.range_bounds_m(2),info.n_draws);
end
fprintf(['Krige_Cell_Agronomic_Values: polygon %d, %s <- "%s" (n=%d; %s; ' ...
    '%d/%d cells inside the sample hull)\n'], ...
    p,stem,column,info.n_samples,range_txt,nnz(r.in_hull),n_cells);
in = r.in_hull & isfinite(r.prediction);
if any(in)
    fprintf(['    in-hull cells: median std %.3g (interpolation %.3g, ' ...
        'range %.3g) on a median prediction of %.3g\n'], ...
        median(r.std(in)),median(r.std_interpolation(in)), ...
        median(r.std_range(in)),median(r.prediction(in)));
end
if ~isempty(info.loocv)
    cv = info.loocv;
    fprintf(['    leave-one-out: RMSE %.3g, bias %+.3g, rms standardised ' ...
        'error %.2f (%.2f against interpolation-only std; ~1 = calibrated, ' ...
        '>1 = overconfident)\n'], cv.rmse,cv.mean_error,cv.rms_z, ...
        cv.rms_z_interpolation);
end
end
