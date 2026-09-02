function [feature_indices,features] = Select_Shapefile_Features( ...
        shape_file,attribute_filters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolve shapefile records from attribute values, the way the notebook
%% selected its fields by Farm and Identifier.
%%
%% attribute_filters is a struct whose fields name shapefile attributes and
%% whose values are the accepted entries, so
%% struct('Farm','AW Bruce','Identifier',{{'AW Bruce 1','AW Bruce 2'}})
%% keeps the records matching both.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(attribute_filters)
    attribute_filters = struct();
end
if ~isfile(shape_file)
    error('Select_Shapefile_Features:MissingShapefile', ...
        '%s does not exist.',shape_file);
end

features = shaperead(shape_file);
keep = true(numel(features),1);
filter_names = fieldnames(attribute_filters);
for k = 1:numel(filter_names)
    name = filter_names{k};
    if ~isfield(features,name)
        error('Select_Shapefile_Features:MissingAttribute', ...
            '%s has no attribute named %s. It carries: %s.', ...
            shape_file,name,strjoin(fieldnames(features),', '));
    end
    accepted = string(attribute_filters.(name));
    present = string({features.(name)}).';
    keep = keep & ismember(present,accepted);
end

feature_indices = reshape(find(keep),1,[]);
if isempty(feature_indices)
    error('Select_Shapefile_Features:NoMatchingFeature', ...
        'No record of %s matches the requested attributes.',shape_file);
end

end
