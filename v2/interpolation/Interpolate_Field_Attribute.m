function result = Interpolate_Field_Attribute(points,field_setup, ...
        attribute,sample_stage,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate one measured attribute at one sample stage over one field.
%%
%% MATLAB form of plot_kriged_attribute from bruce_interpolation.ipynb, with
%% the figure left to the plotting function. The notebook passed its
%% variogram_model argument to the title only and always fitted an
%% exponential model; here the argument reaches the fit, and the default
%% keeps the notebook result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parser = inputParser;
parser.addParameter('StageColumn','Sample Type');
parser.addParameter('XColumn','lng (X)');
parser.addParameter('YColumn','lat (Y)');
parser.addParameter('VariogramModel','exponential');
parser.addParameter('VariogramParameters',struct());
parser.addParameter('GridResolution',200);
parser.parse(varargin{:});
settings = parser.Results;

attribute = char(attribute);
sample_stage = char(sample_stage);
required_columns = {attribute,settings.StageColumn, ...
    settings.XColumn,settings.YColumn};
missing = required_columns(~ismember(required_columns, ...
    points.Properties.VariableNames));
if ~isempty(missing)
    error('Interpolate_Field_Attribute:MissingColumns', ...
        'Missing required columns: %s',strjoin(missing,', '));
end

stage = string(points.(settings.StageColumn));
subset = points(strcmp(stage,sample_stage),:);
if isempty(subset)
    error('Interpolate_Field_Attribute:UnknownStage', ...
        'No sample carries stage %s.',sample_stage);
end
measured = Numeric_Sample_Column(subset,attribute);
located = isfinite(Numeric_Sample_Column(subset,settings.XColumn)) & ...
    isfinite(Numeric_Sample_Column(subset,settings.YColumn)) & ...
    isfinite(measured);
subset = subset(located,:);
if height(subset) < 3
    error('Interpolate_Field_Attribute:TooFewSamples', ...
        ['Only %d valid samples found for %s. ' ...
         'Kriging requires more points'],height(subset),sample_stage);
end

[samples,projection_report] = Project_Sample_Points( ...
    subset,field_setup,settings.XColumn,settings.YColumn);
if height(samples) < 3
    error('Interpolate_Field_Attribute:TooFewSamplesInside', ...
        ['Fewer than 3 samples remain inside the field boundary ' ...
         'for %s at %s.'],attribute,sample_stage);
end

x_m = samples.x_grid_m;
y_m = samples.y_grid_m;
values = Numeric_Sample_Column(samples,attribute);
kriging_model = Fit_Ordinary_Kriging(x_m,y_m,values, ...
    settings.VariogramModel,settings.VariogramParameters);

boundary_polygon = Field_Setup_Polygon(field_setup);
vertices = boundary_polygon.Vertices;
vertices = vertices(all(isfinite(vertices),2),:);
grid_resolution = double(settings.GridResolution);
if ~(isscalar(grid_resolution) && isfinite(grid_resolution) && ...
        grid_resolution >= 2)
    error('Interpolate_Field_Attribute:BadGridResolution', ...
        'GridResolution must be a finite scalar of at least 2.');
end
grid_x = linspace(min(vertices(:,1)),max(vertices(:,1)),grid_resolution);
grid_y = linspace(min(vertices(:,2)),max(vertices(:,2)),grid_resolution);
[prediction,variance] = Krige_Attribute_Grid( ...
    kriging_model,grid_x,grid_y,boundary_polygon);

result.attribute = attribute;
result.sample_stage = sample_stage;
result.grid_x = grid_x;
result.grid_y = grid_y;
result.prediction = prediction;
result.variance = variance;
result.samples = samples;
result.sample_values = values;
result.kriging_model = kriging_model;
result.boundary_polygon = boundary_polygon;
result.stage_sample_count = height(subset);
result.inside_sample_count = height(samples);
result.projected_crs_epsg = projection_report.projected_crs_epsg;
result.value_limits = [min(values) max(values)];

end
