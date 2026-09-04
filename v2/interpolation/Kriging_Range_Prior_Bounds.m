function [bounds_m,diag] = Kriging_Range_Prior_Bounds(x_m,y_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kriging_Range_Prior_Bounds
%%
%% The interval of variogram ranges a sample cloud can actually speak to,
%% used as the support of the range prior in Krige_Attribute_Range_Ensemble
%% when the sampling is too sparse to FIT a range (the usual case for a
%% ten-sample control field; see Fit_Empirical_Variogram_Range .resolved).
%% Both bounds are derived from the sample geometry alone, never assumed:
%%
%%   lower = median nearest-neighbour spacing.  An exponential model whose
%%           effective range equals that spacing already correlates
%%           neighbours at only e^-3 ~ 5%; any shorter range gives the same
%%           (pure-nugget) predictions, so extending the prior below it adds
%%           no distinct hypothesis, it only over-weights that limit.
%%   upper = sample-cloud diameter (longest pairwise separation).  A range
%%           longer than the cloud is indistinguishable, over the sampled
%%           domain, from "no decay across the field".
%%
%% OUTPUTS
%%   bounds_m  [1 x 2]  [lower upper] in the samples' length unit
%%   diag      struct: nearest_neighbour_m (per sample), diameter_m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(double(x_m),[],1);
y = reshape(double(y_m),[],1);
if numel(x) ~= numel(y) || numel(x) < 2
    error('Kriging_Range_Prior_Bounds:BadInput', ...
        'x_m and y_m must be equal-length vectors of at least 2 samples.');
end
if ~all(isfinite(x) & isfinite(y))
    error('Kriging_Range_Prior_Bounds:NonFinite', ...
        'Every sample must have finite coordinates.');
end

d = hypot(x-x.',y-y.');
diameter = max(d(:));
d(1:numel(x)+1:end) = inf;
nn = min(d,[],2);

positive = nn(nn > 0);
if isempty(positive) || ~(diameter > 0)
    error('Kriging_Range_Prior_Bounds:CoincidentSamples', ...
        'All samples share one location; no spatial scale exists.');
end
lower = median(positive);   % duplicates at zero spacing carry no scale
upper = max(diameter,lower);

bounds_m = [lower upper];
diag = struct('nearest_neighbour_m',nn,'diameter_m',diameter);
end
