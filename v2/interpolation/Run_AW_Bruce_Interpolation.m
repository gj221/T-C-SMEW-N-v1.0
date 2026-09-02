function [report,figure_files] = Run_AW_Bruce_Interpolation(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map one soil attribute at every sampling timepoint over the six AW Bruce
%% fields, one calibrated grid per field.
%%
%% This is the simulation view of the notebook's site: AW Bruce 1, 2, 4, 6,
%% 12 and 13 from AdamWilbourne_Fields.shp, each gridded and interpolated on
%% its own, the way a 2D run would see it. Samples come from sheet 1 of the
%% agronomic workbook, the notebook's bruce.csv. For the single unioned
%% surface the notebook actually drew, use Run_Notebook_Interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(script_directory,fullfile(v2_directory,'T&C_Code'));

shape_file = Bruce_Farm_Input_File( ...
    fullfile('Whole farm','Map files','AdamWilbourne_Fields.shp'), ...
    fullfile('Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Fields.shp'));
sample_file = Bruce_Farm_Input_File(fullfile('Whole farm','Data', ...
    'RJ_msft_yale_geochem - agronomic.xlsx'));
attribute_filters = struct('Farm','AW Bruce','Identifier', ...
    {{'AW Bruce 1','AW Bruce 2','AW Bruce 4','AW Bruce 6', ...
    'AW Bruce 12','AW Bruce 13'}});
feature_indices = Select_Shapefile_Features(shape_file,attribute_filters);

[report,figure_files] = Run_Field_Interpolation_Maps( ...
    shape_file,sample_file,'SampleSheet',1, ...
    'FeatureIndices',feature_indices, ...
    'SummaryName','AW_Bruce_Interpolation_Timepoints',varargin{:});

end
