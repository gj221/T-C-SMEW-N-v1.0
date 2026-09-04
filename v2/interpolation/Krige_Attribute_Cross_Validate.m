function cv = Krige_Attribute_Cross_Validate(x_m,y_m,values,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Krige_Attribute_Cross_Validate
%%
%% Leave-one-out cross-validation of the kriging estimator exactly as it
%% is deployed by Krige_Attribute_Range_Ensemble: each sample in turn is
%% withheld, predicted from the others (same variogram model, same fixed
%% sill/nugget, same range prior and the same random draws), and its
%% error is standardised by the estimator's own reported std.  This is the
%% empirical check that the model-internal uncertainty is honest -- the
%% kriging variance is a statement of the variogram, not of the data, so
%% with an assumed range it can be arbitrarily overconfident.
%%
%%   rms_z  = sqrt(mean(z^2)),  z = (predicted - observed)/std
%%          ~ 1   calibrated
%%          > 1   overconfident (errors larger than the std claims)
%%          < 1   conservative
%%   mean_z ~ 0   unbiased
%%
%% rms_z_interpolation is the same statistic against the interpolation-
%% only std (range treated as known): comparing the two shows how much
%% of the calibration the range term restores.
%%
%% Name/value: as Krige_Attribute_Range_Ensemble.  Each fold is the
%% estimator AS DEPLOYED ON THE n-1 RETAINED SAMPLES: sill, nugget and
%% the range prior are all re-derived from those samples, never from the
%% full set.  Sharing the full-set sill would leak the held-out value
%% into every fold's variance and make rms_z one-sidedly optimistic.
%% A .sill/.nugget given explicitly in VariogramParameters is honoured
%% (it is then, by the caller's choice, external to the samples).
%%
%% OUTPUT struct cv (n = number of samples, needs n >= 5 so that every
%% fold keeps >= 4 samples and the sill inflation stays defined)
%%   .observed .predicted .error .std .std_interpolation .z
%%   .z_interpolation                                          [n x 1]
%%   .rmse .mean_error .rms_z .mean_z .rms_z_interpolation
%%   .n .range_bounds_m (of the full set, for reference) .n_draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(here);

opt = struct('VariogramModel','exponential','VariogramParameters',struct(), ...
    'RangeDraws',200,'RangeBoundsM',[],'Seed',0,'SillUncertainty',true);
known = fieldnames(opt);
for k = 1:2:numel(varargin)
    if ~ismember(varargin{k},known)
        error('Krige_Attribute_Cross_Validate:UnknownOption', ...
            'Unknown option "%s".',char(string(varargin{k})));
    end
    opt.(varargin{k}) = varargin{k+1};
end

x = reshape(double(x_m),[],1); y = reshape(double(y_m),[],1);
v = reshape(double(values),[],1);
if ~isequal(numel(x),numel(y),numel(v))
    error('Krige_Attribute_Cross_Validate:SizeMismatch', ...
        'x_m, y_m and values must hold the same number of samples.');
end
if ~all(isfinite(x) & isfinite(y) & isfinite(v))
    error('Krige_Attribute_Cross_Validate:NonFiniteSample', ...
        'Every sample must have a finite coordinate pair and value.');
end
n = numel(v);
if n < 5
    error('Krige_Attribute_Cross_Validate:TooFewSamples', ...
        ['Leave-one-out kriging needs at least 5 samples (4 in every ' ...
         'fold), got %d.'],n);
end

%% Only an explicitly supplied sill/nugget/range is shared across folds;
%% everything else is derived inside each fold from its own samples.
vp = opt.VariogramParameters;
bounds_full = opt.RangeBoundsM;
if isempty(bounds_full) && opt.RangeDraws > 0
    bounds_full = Kriging_Range_Prior_Bounds(x,y);   % reported, not used
end

pred = zeros(n,1); sd = zeros(n,1); sd_i = zeros(n,1); n_draws = 0;
for i = 1:n
    keep = true(n,1); keep(i) = false;
    e = Krige_Attribute_Range_Ensemble(x(keep),y(keep),v(keep),x(i),y(i), ...
        'VariogramModel',opt.VariogramModel,'VariogramParameters',vp, ...
        'RangeDraws',opt.RangeDraws,'RangeBoundsM',opt.RangeBoundsM, ...
        'Seed',opt.Seed,'SillUncertainty',opt.SillUncertainty);
    pred(i) = e.prediction; sd(i) = e.std; sd_i(i) = e.std_interpolation;
    n_draws = e.n_draws;
end
bounds = bounds_full;

err = pred - v;
z   = err./max(sd,eps);
z_i = err./max(sd_i,eps);

cv = struct();
cv.observed = v;   cv.predicted = pred;   cv.error = err;
cv.std = sd;       cv.std_interpolation = sd_i;
cv.z = z;          cv.z_interpolation = z_i;
cv.rmse       = sqrt(mean(err.^2));
cv.mean_error = mean(err);
cv.rms_z      = sqrt(mean(z.^2));
cv.mean_z     = mean(z);
cv.rms_z_interpolation = sqrt(mean(z_i.^2));
cv.n = n;
cv.range_bounds_m = bounds;
cv.n_draws = n_draws;
end
