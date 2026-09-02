function [report,figure_files] = Run_Yield_Field_Interpolation(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map one soil attribute at every sampling timepoint over the Yield Field
%% applied strips.
%%
%% The shapefile holds each strip four times over; the connected-region
%% split collapses those to the three distinct strips, which are gridded and
%% interpolated separately. Samples come from the Yield Field soil workbook,
%% which carries BL, BLP, R1, R2 and R3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape_file = Bruce_Farm_Input_File( ...
    fullfile('Yield field','Map files','applied_yield_shp.shp'), ...
    fullfile('Spatial inputs','Bruce Farm Virginia','applied_yield_shp.shp'));
sample_file = Bruce_Farm_Input_File(fullfile('Yield field','Data', ...
    'Yield Field - BL, R1, R2, R3 - Soil Data.xlsx'));

[report,figure_files] = Run_Field_Interpolation_Maps( ...
    shape_file,sample_file,'SampleSheet',1,'FeatureIndices',[1 2 3], ...
    'SummaryName','Yield_Field_Interpolation_Timepoints',varargin{:});

end
