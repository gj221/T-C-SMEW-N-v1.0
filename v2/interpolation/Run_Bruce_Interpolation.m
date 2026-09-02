function [report,figure_files] = Run_Bruce_Interpolation(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map one soil attribute at every sampling timepoint over the whole-site
%% control polygons.
%%
%% MATLAB form of bruce_interpolation.ipynb. The boundary comes from the
%% whole-farm shapefile and the samples from the agronomic workbook. Only
%% its second sheet reaches the control polygons: the Bruce sheet covers the
%% Backside, Frontside and Yield Field block, more than a kilometre away.
%% Source feature 3, Control 3, belongs to Field 7 and is left out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape_file = Bruce_Farm_Input_File( ...
    fullfile('Whole farm','Map files','AdamWilbourne_Control_Fields.shp'), ...
    fullfile('Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp'));
sample_file = Bruce_Farm_Input_File(fullfile('Whole farm','Data', ...
    'RJ_msft_yale_geochem - agronomic.xlsx'));

[report,figure_files] = Run_Field_Interpolation_Maps( ...
    shape_file,sample_file,'SampleSheet',2,'FeatureIndices',[1 2], ...
    'SummaryName','Bruce_Control_Interpolation_Timepoints',varargin{:});

end
