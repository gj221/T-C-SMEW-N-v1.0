function [report,figure_files] = Run_Field_Interpolation_Maps( ...
        shape_file,sample_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Krige one attribute at every sampling timepoint over every polygon.
%%
%% Each disconnected polygon keeps its own grid, its own samples and its own
%% surface, and each sampling stage is interpolated separately so a
%% timepoint is never pooled with another. Stages that leave a polygon with
%% fewer than three usable samples inside the boundary are reported rather
%% than silently dropped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(script_directory,fullfile(v2_directory,'T&C_Code'));

parser = inputParser;
parser.addParameter('SampleSheet',1);
parser.addParameter('FeatureIndices',[]);
parser.addParameter('CellSizeM',[]);
parser.addParameter('Attribute','Cation Exchange Capacity');
parser.addParameter('SampleStages',[]);
parser.addParameter('StageColumn','Sample Type');
parser.addParameter('XColumn','lng (X)');
parser.addParameter('YColumn','lat (Y)');
parser.addParameter('VariogramModel','exponential');
parser.addParameter('VariogramParameters',struct());
parser.addParameter('GridResolution',200);
parser.addParameter('SummaryName','Field_Interpolation_Timepoints');
parser.parse(varargin{:});
settings = parser.Results;
attribute = char(settings.Attribute);

points = Read_Sample_Points(sample_file,settings.SampleSheet, ...
    settings.XColumn,settings.YColumn);
setups = Build_Shapefile_Grid_Setups( ...
    shape_file,settings.CellSizeM,settings.FeatureIndices);
features = shaperead(shape_file);
field_name = strings(numel(setups),1);
for k = 1:numel(setups)
    field_name(k) = Field_Setup_Label( ...
        features,setups(k).source_feature_indices,k);
end

if ~ismember(settings.StageColumn,points.Properties.VariableNames)
    error('Run_Field_Interpolation_Maps:MissingStageColumn', ...
        'The sample table has no column named %s.',settings.StageColumn);
end

polygon_index = [];
polygon_name = strings(0,1);
stage_name = strings(0,1);
inside_count = [];
usable_count = [];
for k = 1:numel(setups)
    inside_samples = Project_Sample_Points(points,setups(k), ...
        settings.XColumn,settings.YColumn);
    if isempty(inside_samples); continue; end
    inside_stage = string(inside_samples.(settings.StageColumn));
    usable = isfinite(Numeric_Sample_Column(inside_samples,attribute));
    for stage = reshape(unique(inside_stage),1,[])
        selected = inside_stage == stage;
        polygon_index(end+1,1) = k; %#ok<AGROW>
        polygon_name(end+1,1) = field_name(k); %#ok<AGROW>
        stage_name(end+1,1) = stage; %#ok<AGROW>
        inside_count(end+1,1) = sum(selected); %#ok<AGROW>
        usable_count(end+1,1) = sum(selected & usable); %#ok<AGROW>
    end
end
if isempty(stage_name)
    error('Run_Field_Interpolation_Maps:NoSamplesInside', ...
        'No sample in %s falls inside any selected polygon.',sample_file);
end

if isempty(settings.SampleStages)
    stages = unique(stage_name(usable_count >= 3),'stable');
else
    stages = reshape(string(settings.SampleStages),[],1);
end
if isempty(stages)
    error('Run_Field_Interpolation_Maps:NoUsableStage', ...
        ['No sampling stage leaves 3 or more finite %s values inside ' ...
         'a polygon.'],attribute);
end
stages = sort(stages);

status = repmat("not attempted",numel(stage_name),1);
cell_count = NaN(numel(stage_name),1);
cell_size_m = NaN(numel(stage_name),1);
value_min = NaN(numel(stage_name),1);
value_max = NaN(numel(stage_name),1);
value_mean = NaN(numel(stage_name),1);
surface_min = NaN(numel(stage_name),1);
surface_max = NaN(numel(stage_name),1);
figure_files = strings(0,1);

for stage = reshape(stages,1,[])
    fprintf('\n--- %s at %s ---\n',attribute,stage);
    output_name = Field_Grid_Output_Name(shape_file,settings.CellSizeM, ...
        settings.FeatureIndices,attribute,stage);
    [~,~,figure_handle,interpolations] = Plot_Field_Grid_Interpolation( ...
        shape_file,settings.CellSizeM,settings.FeatureIndices, ...
        'Points',points,'Attribute',attribute,'SampleStage',char(stage), ...
        'StageColumn',settings.StageColumn,'XColumn',settings.XColumn, ...
        'YColumn',settings.YColumn, ...
        'VariogramModel',settings.VariogramModel, ...
        'VariogramParameters',settings.VariogramParameters, ...
        'GridResolution',settings.GridResolution, ...
        'OutputName',output_name);
    close(figure_handle);
    figure_files(end+1,1) = string(output_name); %#ok<AGROW>
    for k = 1:numel(setups)
        row = find(polygon_index == k & stage_name == stage,1);
        if isempty(row); continue; end
        cell_count(row) = setups(k).N_cells;
        cell_size_m(row) = setups(k).cell_size_m;
        if isempty(interpolations{k})
            status(row) = "skipped, fewer than 3 usable samples";
            continue
        end
        status(row) = "interpolated";
        value_min(row) = interpolations{k}.value_limits(1);
        value_max(row) = interpolations{k}.value_limits(2);
        value_mean(row) = mean(interpolations{k}.sample_values);
        surface = interpolations{k}.prediction(:);
        surface = surface(isfinite(surface));
        surface_min(row) = min(surface);
        surface_max(row) = max(surface);
        fprintf(['  %-12s n=%2d  cells=%5d at %.3f m  ' ...
            'samples %.4g to %.4g  surface %.4g to %.4g\n'], ...
            polygon_name(row),interpolations{k}.inside_sample_count, ...
            setups(k).N_cells,setups(k).cell_size_m, ...
            value_min(row),value_max(row),surface_min(row),surface_max(row));
    end
end

report = table(polygon_index,polygon_name,stage_name, ...
    repmat(string(attribute),numel(stage_name),1),inside_count, ...
    usable_count,status,cell_count,cell_size_m,value_min,value_max, ...
    value_mean,surface_min,surface_max, ...
    'VariableNames',{'Polygon','Field','SampleStage','Attribute', ...
    'SamplesInside','UsableSamples','Status','CellCount','CellSizeM', ...
    'SampleMin','SampleMax','SampleMean','SurfaceMin','SurfaceMax'});
report = sortrows(report,{'Polygon','SampleStage'});

figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
writetable(report,fullfile(figure_directory, ...
    [char(settings.SummaryName) '.csv']));

end

