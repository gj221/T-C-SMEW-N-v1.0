function [results,site] = Run_Notebook_Interpolation(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The whole-site figure bruce_interpolation.ipynb produced, in MATLAB.
%%
%% The notebook selected the AW Bruce fields from AdamWilbourne_Fields.shp,
%% took their union as one boundary, fitted a single ordinary kriging model
%% to every sample inside it and drew one surface in Easting and Northing.
%% That is what this reproduces, one figure per sampling stage.
%%
%% Records 65, 77, 115, 116, 117 and 118 of AdamWilbourne_Fields.shp are AW
%% Bruce 4, 2, 1, 6, 13 and 12, the notebook's six fields. Sheet 1 of the
%% agronomic workbook is the notebook's bruce.csv: 841 of its 1871 rows fall
%% inside these fields. If the shapefile is missing this falls back to the
%% Yield Field applied strips and says so; pass ShapeFile to point elsewhere.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(script_directory,fullfile(v2_directory,'T&C_Code'));

parser = inputParser;
parser.addParameter('ShapeFile','');
parser.addParameter('AttributeFilters',[]);
parser.addParameter('SampleFile','');
parser.addParameter('SampleSheet',[]);
parser.addParameter('Attribute','Cation Exchange Capacity');
parser.addParameter('SampleStages',[]);
parser.addParameter('StageColumn','Sample Type');
parser.addParameter('VariogramModel','exponential');
parser.addParameter('VariogramParameters',struct());
parser.addParameter('GridResolution',200);
parser.parse(varargin{:});
settings = parser.Results;
attribute = char(settings.Attribute);

shape_file = settings.ShapeFile;
attribute_filters = settings.AttributeFilters;
sample_file = settings.SampleFile;
sample_sheet = settings.SampleSheet;
if isempty(shape_file)
    % Prefer the copy tracked in this repository, so a clone reproduces the
    % figure without the survey folder; fall back to the survey folder for
    % a shapefile update the repository has not picked up yet.
    field_file = '';
    tracked_field_file = fullfile(v2_directory,'Spatial inputs', ...
        'Bruce Farm Virginia','Whole farm','AdamWilbourne_Fields.shp');
    survey_field_file = fullfile(Bruce_Farm_Survey_Root,'Whole farm', ...
        'Map files','AdamWilbourne_Fields.shp');
    if isfile(tracked_field_file)
        field_file = tracked_field_file;
    elseif isfile(survey_field_file)
        field_file = survey_field_file;
    end
    if ~isempty(field_file)
        shape_file = field_file;
        if isempty(attribute_filters)
            attribute_filters = struct('Farm','AW Bruce','Identifier', ...
                {{'AW Bruce 1','AW Bruce 2','AW Bruce 4','AW Bruce 6', ...
                'AW Bruce 12','AW Bruce 13'}});
        end
        if isempty(sample_sheet); sample_sheet = 1; end
    else
        fprintf(['AdamWilbourne_Fields.shp is not tracked in this ' ...
            'repository and not found at %s, so the notebook ' ...
            'selection cannot be made.\n'],survey_field_file);
        fprintf(['Falling back to the Yield Field applied strips, the ' ...
            'only compact well-sampled boundary here.\n']);
        fprintf('Pass ShapeFile to use the field boundaries instead.\n');
        shape_file = Bruce_Farm_Input_File(fullfile('Yield field', ...
            'Map files','applied_yield_shp.shp'), ...
            fullfile('Spatial inputs','Bruce Farm Virginia', ...
            'applied_yield_shp.shp'));
        attribute_filters = [];
        if isempty(sample_file)
            sample_file = Bruce_Farm_Input_File(fullfile('Yield field', ...
                'Data','Yield Field - BL, R1, R2, R3 - Soil Data.xlsx'));
        end
        if isempty(sample_sheet); sample_sheet = 1; end
    end
end
if isempty(sample_sheet); sample_sheet = 1; end
if isempty(sample_file)
    sample_file = Bruce_Farm_Input_File(fullfile('Whole farm','Data', ...
        'RJ_msft_yale_geochem - agronomic.xlsx'));
end

if isempty(attribute_filters)
    feature_indices = [];
else
    feature_indices = Select_Shapefile_Features(shape_file, ...
        attribute_filters);
end
site = Build_Site_Boundary(shape_file,feature_indices);
fprintf('%s: %d records, %d disjoint fields, %.1f m2 in %s\n', ...
    shape_file,numel(site.source_feature_indices),site.region_count, ...
    site.field_area_m2,site.projected_crs_name);

points = Read_Sample_Points(sample_file,sample_sheet);
inside_samples = Project_Sample_Points(points,site);
if isempty(inside_samples)
    error('Run_Notebook_Interpolation:NoSamplesInside', ...
        ['None of the %d located samples in %s fall inside the ' ...
         'selected fields.'],height(points),sample_file);
end
fprintf('%d of %d located samples fall inside the site\n', ...
    height(inside_samples),height(points));

usable = isfinite(Numeric_Sample_Column(inside_samples,attribute));
stage_of_sample = string(inside_samples.(settings.StageColumn));
if isempty(settings.SampleStages)
    stages = unique(stage_of_sample(usable));
    counts = arrayfun(@(s) sum(usable & stage_of_sample == s),stages);
    stages = stages(counts >= 3);
else
    stages = reshape(string(settings.SampleStages),[],1);
end
if isempty(stages)
    error('Run_Notebook_Interpolation:NoUsableStage', ...
        'No stage leaves 3 or more finite %s values inside the site.', ...
        attribute);
end

results = cell(numel(stages),1);
for k = 1:numel(stages)
    [result,figure_handle] = Plot_Site_Interpolation(site,points, ...
        attribute,char(stages(k)),'StageColumn',settings.StageColumn, ...
        'VariogramModel',settings.VariogramModel, ...
        'VariogramParameters',settings.VariogramParameters, ...
        'GridResolution',settings.GridResolution);
    close(figure_handle);
    results{k} = result;
    surface_values = result.prediction(isfinite(result.prediction));
    fprintf(['  %-5s n=%3d  samples %.4g to %.4g  ' ...
        'surface %.4g to %.4g  -> %s.png\n'],stages(k), ...
        result.inside_sample_count,result.value_limits(1), ...
        result.value_limits(2),min(surface_values),max(surface_values), ...
        result.output_name);
end

end

function survey_root = Bruce_Farm_Survey_Root()
survey_root = getenv('BRUCE_FARM_ROOT');
if isempty(survey_root)
    survey_root = fullfile(getenv('USERPROFILE'), ...
        'OneDrive - Yale University','Documents','Bruce Farm_Virginia');
end
end
