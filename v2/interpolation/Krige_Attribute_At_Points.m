function [prediction,variance,nearest_sample_m,in_hull] = ...
    Krige_Attribute_At_Points(model,x_m,y_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Krige_Attribute_At_Points
%%
%% Evaluate a fitted ordinary kriging model (Fit_Ordinary_Kriging) at
%% arbitrary scattered target points, e.g. a field's simulation cell
%% centres, rather than the regular display mesh Krige_Attribute_Grid
%% draws. Identical block-solve math to Krige_Attribute_Grid (same
%% exactness at a sample location, same variance definition); this is its
%% scattered-point sibling, not a replacement -- Krige_Attribute_Grid
%% keeps doing the figures.
%%
%% Also reports, per target, the distance to its nearest sample and
%% whether it falls inside the convex hull of the samples used to fit the
%% model, since a prediction outside that hull is an extrapolation the
%% variogram was never fitted to support.
%%
%% INPUTS
%%   model   struct from Fit_Ordinary_Kriging
%%   x_m,y_m target point coordinates, same projected/local frame the
%%           model was fit in
%%
%% OUTPUTS
%%   prediction, variance   [N x 1], NaN nowhere (unlike Krige_Attribute_Grid,
%%                          there is no boundary polygon to clip to here --
%%                          the caller's target points are already the
%%                          cells that matter)
%%   nearest_sample_m       [N x 1] distance to the closest fitting sample
%%   in_hull                [N x 1] logical, true if the target lies inside
%%                          (or on) the convex hull of the fitting samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(model) || ~isfield(model,'system_matrix') || ...
        ~isfield(model,'values')
    error('Krige_Attribute_At_Points:BadModel', ...
        'model must come from Fit_Ordinary_Kriging.');
end
x_m = reshape(double(x_m),[],1);
y_m = reshape(double(y_m),[],1);
if numel(x_m) ~= numel(y_m)
    error('Krige_Attribute_At_Points:SizeMismatch', ...
        'x_m and y_m must have the same number of elements.');
end
if ~all(isfinite(x_m) & isfinite(y_m))
    error('Krige_Attribute_At_Points:NonFiniteTarget', ...
        'Every target point must have finite coordinates.');
end

target_count = numel(x_m);
n = model.sample_count;
prediction = zeros(target_count,1);
variance = zeros(target_count,1);

block_size = max(1,floor(4e6/(n+1)));
for first = 1:block_size:target_count
    last = min(target_count,first+block_size-1);
    block = (first:last).';
    separation_m = hypot(model.x_m-reshape(x_m(block),1,[]), ...
        model.y_m-reshape(y_m(block),1,[]));
    right_hand_side = ones(n+1,numel(block));
    right_hand_side(1:n,:) = Variogram_Semivariance( ...
        model.variogram_model,model.model_parameters,separation_m);
    right_hand_side([separation_m <= model.exact_tolerance_m; ...
        false(1,numel(block))]) = 0;
    weights = model.system_matrix\right_hand_side;
    prediction(block) = weights(1:n,:).'*model.values;
    variance(block) = sum(weights.*right_hand_side,1).';
end

if nargout >= 3
    d = hypot(model.x_m-x_m.',model.y_m-y_m.');   % [n x target_count]
    nearest_sample_m = min(d,[],1).';
end
if nargout >= 4
    if n >= 3
        try
            k = convhull(model.x_m,model.y_m);
            hull = polyshape(model.x_m(k),model.y_m(k));
            in_hull = isinterior(hull,x_m,y_m) | ...
                on_hull_boundary(hull,x_m,y_m);
        catch
            %% Degenerate (collinear) sample layout: no interior exists.
            in_hull = false(target_count,1);
        end
    else
        in_hull = false(target_count,1);
    end
end
end

%% =======================================================================
function tf = on_hull_boundary(hull,x_m,y_m)
%% isinterior(polyshape,...) can report a point exactly on an edge as
%% outside; catch that sliver case by a direct distance-to-boundary test.
tol = 1e-9*max(1,max(abs(hull.Vertices(:))));
tf = false(size(x_m));
[xb,yb] = boundary(hull);
for i = 1:numel(x_m)
    tf(i) = min(hypot(xb-x_m(i),yb-y_m(i))) <= tol;
end
end
