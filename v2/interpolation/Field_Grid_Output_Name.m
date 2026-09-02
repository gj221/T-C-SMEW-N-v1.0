function output_name = Field_Grid_Output_Name(shape_file,cell_size_m, ...
        feature_indices,attribute,sample_stage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the figure and summary file stem for one grid or interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4; attribute = ''; end
if nargin < 5; sample_stage = ''; end

[~,shape_name] = fileparts(shape_file);
shape_name = regexprep(shape_name,'[^A-Za-z0-9]+','_');
feature_tag = '';
if ~isempty(feature_indices)
    parts = arrayfun(@(x) sprintf('%d',x),feature_indices, ...
        'UniformOutput',false);
    feature_tag = sprintf('_features%s',strjoin(parts,'_'));
end
if isempty(cell_size_m)
    resolution = 'geometry';
else
    resolution = strrep(sprintf('cell%gm',cell_size_m),'.','p');
end
attribute_tag = '';
if ~isempty(attribute) && ~isempty(sample_stage)
    attribute_tag = sprintf('_%s_%s', ...
        regexprep(char(attribute),'[^A-Za-z0-9]+','_'), ...
        regexprep(char(sample_stage),'[^A-Za-z0-9]+','_'));
end
output_name = sprintf('Field_Grid_Interpolation_%s%s_%s%s', ...
    shape_name,feature_tag,resolution,attribute_tag);

end
