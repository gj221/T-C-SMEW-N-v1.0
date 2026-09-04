function ens = Krige_Attribute_Range_Ensemble(x_m,y_m,values,xt_m,yt_m,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Krige_Attribute_Range_Ensemble
%%
%% Ordinary kriging at scattered targets under an UNCERTAIN variogram
%% range, by Monte Carlo.  A single Fit_Ordinary_Kriging call treats its
%% range as known -- and when the sampling is too sparse to fit one (a
%% ten-sample field cannot resolve a semivariogram, Fit_Empirical_
%% Variogram_Range .resolved = false) the range is an assumption, and the
%% kriging variance it reports carries none of that assumption's
%% uncertainty.  Here the range is instead drawn repeatedly from a prior
%% whose support the sample geometry itself sets (Kriging_Range_Prior_
%% Bounds: median nearest-neighbour spacing .. cloud diameter; log-uniform,
%% the scale-invariant choice for a length parameter), the field is kriged
%% at every draw, and the draws are combined by the law of total variance
%%
%%     Var[Z]  =  E[ Var(Z | range) ]  +  Var( E[Z | range] )
%%              = variance_interpolation +  variance_range
%%
%% so the reported uncertainty prices in both the ordinary interpolation
%% error at a given range and how far the prediction itself moves as the
%% assumed range moves.  The sill and nugget are held fixed (they scale the
%% variance but, at a fixed nugget/sill ratio, not the kriging weights);
%% only the range -- the parameter with no empirical support -- is
%% randomised.  Every draw is an exact interpolator, so the ensemble stays
%% exact (zero total variance) at a sample location.
%%
%% INPUTS
%%   x_m,y_m,values   fitting samples (same frame as the targets)
%%   xt_m,yt_m        target points
%%
%% Name/value
%%   'VariogramModel'       default 'exponential'
%%   'VariogramParameters'  struct passed to Fit_Ordinary_Kriging: .sill,
%%                          .nugget (defaults: population variance, 5% of
%%                          it).  A .range is honoured only when RangeDraws
%%                          is 0 (single fit).
%%   'RangeDraws'           Monte Carlo draws (default 200).  0 = one fit
%%                          with the given/default range: the classical
%%                          result, with ens.range_assumed = true.
%%   'RangeBoundsM'         [lower upper] prior support (default: from
%%                          Kriging_Range_Prior_Bounds on the samples)
%%   'Seed'                 RandStream seed for reproducible draws (default 0)
%%   'SillUncertainty'      inflate E[Var|range] by nu/(nu-2), nu = n-1,
%%                          the mean of the sill's scaled inverse-chi-
%%                          square sampling distribution (default true;
%%                          needs n >= 4).  false keeps the sample
%%                          variance as the sill.
%%
%% OUTPUT struct ens
%%   .prediction               [Nt x 1] ensemble-mean prediction
%%   .variance .std            total (law of total variance)
%%   .variance_interpolation   E[kriging variance | range]
%%   .variance_range           Var(prediction | range) across draws
%%   .std_interpolation .std_range
%%   .nearest_sample_m .in_hull   as Krige_Attribute_At_Points
%%   .range_draws_m            [D x 1] the ranges used
%%   .range_bounds_m  .n_draws  .range_assumed
%%   .variogram_model .sill .nugget .n_samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(here);

opt = struct('VariogramModel','exponential','VariogramParameters',struct(), ...
    'RangeDraws',200,'RangeBoundsM',[],'Seed',0,'SillUncertainty',true);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

x = reshape(double(x_m),[],1); y = reshape(double(y_m),[],1);
v = reshape(double(values),[],1);
xt = reshape(double(xt_m),[],1); yt = reshape(double(yt_m),[],1);
if ~isequal(numel(x),numel(y),numel(v))
    error('Krige_Attribute_Range_Ensemble:SizeMismatch', ...
        'x_m, y_m and values must hold the same number of samples.');
end
if ~all(isfinite(x) & isfinite(y) & isfinite(v))
    error('Krige_Attribute_Range_Ensemble:NonFiniteSample', ...
        'Every sample must have a finite coordinate pair and value.');
end
if numel(xt) ~= numel(yt)
    error('Krige_Attribute_Range_Ensemble:TargetMismatch', ...
        'xt_m and yt_m must have the same number of elements.');
end
n_draws = round(double(opt.RangeDraws));
if ~(isscalar(n_draws) && n_draws >= 0)
    error('Krige_Attribute_Range_Ensemble:BadDraws', ...
        'RangeDraws must be a non-negative integer.');
end
if ~isstruct(opt.VariogramParameters) || ~isscalar(opt.VariogramParameters)
    error('Krige_Attribute_Range_Ensemble:BadParameters', ...
        'VariogramParameters must be a scalar struct.');
end

%% Fix the sill and nugget once, from the fitting set -- Fit_Ordinary_
%% Kriging's own defaults -- so every draw shares one variance scale.
vp = opt.VariogramParameters;
if ~isfield(vp,'sill')   || isempty(vp.sill);   vp.sill   = var(v,1);       end
if ~isfield(vp,'nugget') || isempty(vp.nugget); vp.nugget = 0.05*vp.sill;   end

Nt = numel(xt);
model_name = lower(char(opt.VariogramModel));
if n_draws > 0 && any(strcmp(model_name,{'linear','power'}))
    %% These models have no range to integrate out; Fit_Ordinary_Kriging
    %% would ignore .range and .sill, and perturbing only the nugget would
    %% move the weights -- the one thing the decomposition assumes fixed.
    error('Krige_Attribute_Range_Ensemble:NoRange', ...
        ['The %s variogram has no range parameter; use RangeDraws = 0 ' ...
         'or an exponential / gaussian / spherical model.'],model_name);
end
if ~(vp.sill > 0)
    %% Every sample carries the same value: the field is that constant,
    %% with no uncertainty.  A zero sill makes the kriging system singular,
    %% so answer directly rather than solve it.
    geom = Fit_Ordinary_Kriging(x,y,v,'exponential', ...
        struct('sill',1,'nugget',0.05,'range',max(1,max(hypot(x-x.',y-y.'),[],'all'))));
    [~,~,near_m,in_hull] = Krige_Attribute_At_Points(geom,xt,yt);
    ens = struct('prediction',repmat(v(1),Nt,1), ...
        'variance_interpolation',zeros(Nt,1),'variance_range',zeros(Nt,1), ...
        'variance',zeros(Nt,1),'std',zeros(Nt,1),'std_interpolation',zeros(Nt,1), ...
        'std_range',zeros(Nt,1),'nearest_sample_m',near_m,'in_hull',in_hull, ...
        'range_draws_m',NaN,'range_bounds_m',[NaN NaN],'n_draws',0, ...
        'range_assumed',true,'variogram_model',model_name,'sill',0,'nugget',0, ...
        'n_samples',numel(v),'sill_inflation',1);
    return
end
if n_draws == 0
    %% Classical single fit: the range is whatever the caller (or the
    %% Fit_Ordinary_Kriging default) says -- an assumption, and flagged so.
    model = Fit_Ordinary_Kriging(x,y,v,opt.VariogramModel,vp);
    [P,V,near_m,in_hull] = Krige_Attribute_At_Points(model,xt,yt);
    range_draws = model.model_parameters(2);
    bounds = [range_draws range_draws];
    var_interp = V;
    var_range  = zeros(Nt,1);
    range_assumed = true;
    sill_inflation = 1;
else
    bounds = opt.RangeBoundsM;
    if isempty(bounds)
        bounds = Kriging_Range_Prior_Bounds(x,y);
    end
    bounds = double(bounds(:)).';
    if numel(bounds) ~= 2 || ~all(isfinite(bounds)) || bounds(1) <= 0 || ...
            bounds(2) < bounds(1)
        error('Krige_Attribute_Range_Ensemble:BadBounds', ...
            'RangeBoundsM must be [lower upper] with 0 < lower <= upper.');
    end
    %% Log-uniform draws on [lower, upper]; a private stream so the
    %% caller's global generator is untouched.
    rs = RandStream('mt19937ar','Seed',opt.Seed);
    u = rand(rs,n_draws,1);
    range_draws = bounds(1)*(bounds(2)/bounds(1)).^u;

    %% The SILL is estimated from the same handful of samples, and the
    %% kriging variance is proportional to it, so treating the sample
    %% variance as the truth is a second source of overconfidence.  Its
    %% sampling distribution is the scaled inverse-chi-square,
    %% sill ~ s2 * nu/chi2_nu with nu = n-1, whose mean is nu/(nu-2).
    %% Because scaling sill and nugget together leaves the kriging
    %% WEIGHTS -- and so the prediction -- untouched and scales the
    %% kriging variance linearly, the expectation is available in closed
    %% form: no draw is needed, and a draw would only add sampling noise
    %% (with infinite variance for nu <= 4).  nu = n-1 is optimistic when
    %% the samples are spatially correlated (the effective sample size is
    %% then smaller), so this narrows the calibration gap without
    %% claiming to close it.
    if opt.SillUncertainty && numel(v) >= 4          % nu >= 3 > 2
        nu = numel(v)-1;
        sill_inflation = nu/(nu-2);
    else
        sill_inflation = 1;
    end

    P = zeros(Nt,n_draws); V = zeros(Nt,n_draws);
    near_m = []; in_hull = [];
    for d = 1:n_draws
        vpd = vp; vpd.range = range_draws(d);
        model = Fit_Ordinary_Kriging(x,y,v,opt.VariogramModel,vpd);
        if d == 1
            [P(:,d),V(:,d),near_m,in_hull] = Krige_Attribute_At_Points(model,xt,yt);
        else
            [P(:,d),V(:,d)] = Krige_Attribute_At_Points(model,xt,yt);
        end
    end
    var_interp = sill_inflation*mean(V,2);     % E[ Var(Z | range,sill) ]
    var_range  = var(P,0,2);                   % Var( E[Z | range] )
    P = mean(P,2);
    range_assumed = false;
end

var_interp = max(var_interp,0);
ens = struct();
ens.prediction             = P;
ens.variance_interpolation = var_interp;
ens.variance_range         = var_range;
ens.variance               = var_interp + var_range;
ens.std                    = sqrt(ens.variance);
ens.std_interpolation      = sqrt(var_interp);
ens.std_range              = sqrt(var_range);
ens.nearest_sample_m       = near_m;
ens.in_hull                = in_hull;
ens.range_draws_m          = range_draws(:);
ens.range_bounds_m         = bounds;
ens.n_draws                = n_draws;
ens.range_assumed          = range_assumed;
ens.variogram_model        = lower(char(opt.VariogramModel));
ens.sill                   = vp.sill;
ens.nugget                 = vp.nugget;
ens.n_samples              = numel(v);
ens.sill_inflation         = sill_inflation;   % nu/(nu-2) applied to E[Var|range]
end
