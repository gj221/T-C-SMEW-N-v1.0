function [prediction,variance,grid_x,grid_y] = Krige_Attribute_Grid( ...
        model,grid_x,grid_y,boundary_polygon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate a fitted ordinary kriging model over a regular grid.
%%
%% This is the grid form of pykrige's execute. A target that coincides with
%% a sample takes zero semivariance against it, so the surface honours every
%% measurement exactly instead of returning it with the nugget attached.
%% Cells outside the boundary polygon are returned as NaN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4; boundary_polygon = []; end
if ~isstruct(model) || ~isfield(model,'system_matrix') || ...
        ~isfield(model,'values')
    error('Krige_Attribute_Grid:BadModel', ...
        'model must come from Fit_Ordinary_Kriging.');
end
grid_x = reshape(double(grid_x),1,[]);
grid_y = reshape(double(grid_y),1,[]);
if isempty(grid_x) || isempty(grid_y) || ...
        ~all(isfinite(grid_x)) || ~all(isfinite(grid_y))
    error('Krige_Attribute_Grid:BadGrid', ...
        'grid_x and grid_y must be non-empty finite coordinate vectors.');
end

[target_x,target_y] = meshgrid(grid_x,grid_y);
target_count = numel(target_x);
n = model.sample_count;
prediction = zeros(target_count,1);
variance = zeros(target_count,1);

% One solve per block keeps the right-hand side inside memory on fine grids.
block_size = max(1,floor(4e6/(n+1)));
for first = 1:block_size:target_count
    last = min(target_count,first+block_size-1);
    block = first:last;
    separation_m = hypot(model.x_m-reshape(target_x(block),1,[]), ...
        model.y_m-reshape(target_y(block),1,[]));
    right_hand_side = ones(n+1,numel(block));
    right_hand_side(1:n,:) = Variogram_Semivariance( ...
        model.variogram_model,model.model_parameters,separation_m);
    right_hand_side([separation_m <= model.exact_tolerance_m; ...
        false(1,numel(block))]) = 0;
    weights = model.system_matrix\right_hand_side;
    prediction(block) = weights(1:n,:).'*model.values;
    variance(block) = sum(weights.*right_hand_side,1).';
end

prediction = reshape(prediction,size(target_x));
variance = reshape(variance,size(target_x));
if ~isempty(boundary_polygon)
    if ~(isa(boundary_polygon,'polyshape') && isscalar(boundary_polygon))
        error('Krige_Attribute_Grid:BadBoundary', ...
            'boundary_polygon must be a scalar polyshape or empty.');
    end
    outside = ~isinterior(boundary_polygon,target_x(:),target_y(:));
    prediction(reshape(outside,size(prediction))) = NaN;
    variance(reshape(outside,size(variance))) = NaN;
end

end
