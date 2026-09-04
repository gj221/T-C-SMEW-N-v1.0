function test_kriging_uncertainty()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_kriging_uncertainty
%%
%% Checks on the range-uncertainty Monte Carlo (Krige_Attribute_Range_
%% Ensemble), its prior support (Kriging_Range_Prior_Bounds) and the
%% leave-one-out calibration check (Krige_Attribute_Cross_Validate).
%% Synthetic Gaussian random fields throughout (Cholesky of the true
%% exponential covariance), so the truth is known.
%%
%%   1  RangeDraws = 0 reproduces the classical single fit bit-for-bit
%%   2  a degenerate prior [r r] reproduces the single fit at r
%%   3  the law-of-total-variance identity holds and the range term is live
%%   4  the sill draw inflates the variance by nu/(nu-2); draws are
%%      reproducible by seed, within bounds, log-uniform
%%   5  the ensemble stays exact at a sample location
%%   6  prior bounds from a regular grid are the spacing and the diagonal
%%   7  LOOCV is calibrated at the true range and flags an over-long
%%      (over-confident) assumed range
%%   8  LOOCV of the ensemble estimator runs and reports its prior
%%   9  input guards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_kriging_uncertainty ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
v2   = fileparts(code);
addpath(code, fullfile(v2,'interpolation'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% --- synthetic field: exponential, true range 30 m, sill 1, nugget 0.05
rng(11);
n = 40; L = 100;
xs = L*rand(n,1); ys = L*rand(n,1);
true_range = 30; sill = 1; nugget = 0.05; mu = 6;
vs = mu + gaussian_field(xs,ys,true_range,sill,nugget);
vp = struct('sill',sill,'nugget',nugget);
[XT,YT] = meshgrid(10:20:90,10:20:90); xt = XT(:); yt = YT(:);

%% ---------------------------------------------------------------------
%% 1. RangeDraws = 0 == classical single fit
%% ---------------------------------------------------------------------
vp80 = vp; vp80.range = 80;
model = Fit_Ordinary_Kriging(xs,ys,vs,'exponential',vp80);
[p0,v0] = Krige_Attribute_At_Points(model,xt,yt);
e0 = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp80,'RangeDraws',0);
check('RangeDraws=0 reproduces the single fit exactly', ...
    isequal(e0.prediction,p0) && isequal(e0.variance,v0) && ...
    e0.range_assumed && e0.n_draws == 0 && all(e0.std_range == 0));

%% ---------------------------------------------------------------------
%% 2. Degenerate prior [r r] == single fit at r
%% ---------------------------------------------------------------------
vpr = vp; vpr.range = 25;
[pr,vr] = Krige_Attribute_At_Points( ...
    Fit_Ordinary_Kriging(xs,ys,vs,'exponential',vpr),xt,yt);
er = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',50,'RangeBoundsM',[25 25], ...
    'SillUncertainty',false);
check('degenerate prior [r r], sill fixed, reproduces the single fit at r', ...
    max(abs(er.prediction-pr)) < 1e-10 && ...
    max(abs(er.variance_interpolation-vr)) < 1e-10 && ...
    max(er.std_range) < 1e-10 && ~er.range_assumed);

%% ---------------------------------------------------------------------
%% 3. Law of total variance, range term live
%% ---------------------------------------------------------------------
ew = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',300,'RangeBoundsM',[10 150]);
check('variance = variance_interpolation + variance_range', ...
    max(abs(ew.variance-(ew.variance_interpolation+ew.variance_range))) < 1e-12);
check('the range term is non-zero away from the samples', ...
    max(ew.std_range) > 1e-3, sprintf('max std_range %.3g',max(ew.std_range)));
check('total std >= interpolation-only std everywhere', ...
    all(ew.std >= ew.std_interpolation - 1e-12) && all(isfinite(ew.std)));

es_off = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',400,'RangeBoundsM',[10 150], ...
    'SillUncertainty',false,'Seed',5);
es_on = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',400,'RangeBoundsM',[10 150], ...
    'SillUncertainty',true,'Seed',5);
check('drawing the sill inflates the variance (n-1 dof, factor ~ nu/(nu-2))', ...
    mean(es_on.variance) > mean(es_off.variance), ...
    sprintf('on %.4g vs off %.4g',mean(es_on.variance),mean(es_off.variance)));
expected = (n-1)/((n-1)-2);
ratio = mean(es_on.variance_interpolation)/mean(es_off.variance_interpolation);
check('the inflation matches the scaled inverse-chi-square mean', ...
    abs(ratio-expected) < 0.15*expected, ...
    sprintf('ratio %.3f, expected %.3f',ratio,expected));

%% ---------------------------------------------------------------------
%% 4. Reproducible, bounded, log-uniform draws
%% ---------------------------------------------------------------------
ea = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',400,'RangeBoundsM',[20 80],'Seed',3);
eb = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',400,'RangeBoundsM',[20 80],'Seed',3);
ec = Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramParameters',vp,'RangeDraws',400,'RangeBoundsM',[20 80],'Seed',4);
check('same seed -> identical ensemble', ...
    isequal(ea.prediction,eb.prediction) && isequal(ea.range_draws_m,eb.range_draws_m));
check('different seed -> different draws', ~isequal(ea.range_draws_m,ec.range_draws_m));
check('every draw lies inside the prior bounds', ...
    all(ea.range_draws_m >= 20) && all(ea.range_draws_m <= 80));
lmid = mean(log([20 80]));
check('draws are log-uniform (median of log-range at the log-midpoint)', ...
    abs(median(log(ea.range_draws_m)) - lmid) < 0.15, ...
    sprintf('median log r = %.3f, expected %.3f',median(log(ea.range_draws_m)),lmid));
check('the global rng is not disturbed by the private stream', ...
    isequal(rng_state_after(@() Krige_Attribute_Range_Ensemble( ...
        xs,ys,vs,xt,yt,'VariogramParameters',vp,'RangeDraws',5)), ...
        rng_state_after(@() [])));

%% ---------------------------------------------------------------------
%% 5. Exact at a sample location
%% ---------------------------------------------------------------------
es = Krige_Attribute_Range_Ensemble(xs,ys,vs,xs(5),ys(5), ...
    'VariogramParameters',vp,'RangeDraws',100,'RangeBoundsM',[10 150]);
check('ensemble prediction is exact at a sample location', ...
    abs(es.prediction-vs(5)) < 1e-8, sprintf('|d| = %.3g',abs(es.prediction-vs(5))));
check('ensemble total std is ~0 at a sample location', es.std < 1e-6);

%% ---------------------------------------------------------------------
%% 6. Prior bounds from a regular grid
%% ---------------------------------------------------------------------
[GX,GY] = meshgrid(0:12:48,0:12:48);
[b,bd] = Kriging_Range_Prior_Bounds(GX(:),GY(:));
check('regular-grid lower bound is the grid spacing', abs(b(1)-12) < 1e-9);
check('regular-grid upper bound is the cloud diagonal', ...
    abs(b(2)-hypot(48,48)) < 1e-9 && abs(bd.diameter_m-b(2)) < 1e-9);
[bdup,~] = Kriging_Range_Prior_Bounds([0 0 12 24],[0 0 0 0]);
check('a duplicated sample location does not collapse the lower bound', ...
    bdup(1) == 12);

%% ---------------------------------------------------------------------
%% 7. LOOCV calibration
%% ---------------------------------------------------------------------
vpt = vp; vpt.range = true_range;
cv_true = Krige_Attribute_Cross_Validate(xs,ys,vs, ...
    'VariogramParameters',vpt,'RangeDraws',0);
vpl = vp; vpl.range = 10*true_range;
cv_long = Krige_Attribute_Cross_Validate(xs,ys,vs, ...
    'VariogramParameters',vpl,'RangeDraws',0);
check('LOOCV at the true range is calibrated (rms_z in [0.6, 1.5])', ...
    cv_true.rms_z > 0.6 && cv_true.rms_z < 1.5, sprintf('rms_z = %.2f',cv_true.rms_z));
check('LOOCV is roughly unbiased at the true range (|mean_z| < 0.5)', ...
    abs(cv_true.mean_z) < 0.5, sprintf('mean_z = %.2f',cv_true.mean_z));
check('an over-long assumed range is flagged as overconfident (rms_z rises)', ...
    cv_long.rms_z > cv_true.rms_z, ...
    sprintf('rms_z true %.2f vs 10x range %.2f',cv_true.rms_z,cv_long.rms_z));
check('LOOCV reports one fold per sample', ...
    cv_true.n == n && numel(cv_true.z) == n && all(isfinite(cv_true.z)));

%% ---------------------------------------------------------------------
%% 8. LOOCV of the ensemble estimator
%% ---------------------------------------------------------------------
cv_ens = Krige_Attribute_Cross_Validate(xs,ys,vs, ...
    'VariogramParameters',vp,'RangeDraws',60,'Seed',1);
bfull = Kriging_Range_Prior_Bounds(xs,ys);
check('ensemble LOOCV uses the full-set prior in every fold', ...
    isequal(cv_ens.range_bounds_m,bfull) && cv_ens.n_draws == 60);
check('ensemble LOOCV statistics are finite', ...
    all(isfinite([cv_ens.rmse cv_ens.rms_z cv_ens.mean_z cv_ens.rms_z_interpolation])));
check('the range term never makes the ensemble MORE confident than interpolation alone', ...
    cv_ens.rms_z <= cv_ens.rms_z_interpolation + 1e-12);

%% ---------------------------------------------------------------------
%% 9. Guards
%% ---------------------------------------------------------------------
assert_error(@() Krige_Attribute_Cross_Validate([0 1 2 3],[0 1 2 3],[1 2 3 4]), ...
    'Krige_Attribute_Cross_Validate:TooFewSamples');
assert_error(@() Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'VariogramModel','linear','RangeDraws',10), ...
    'Krige_Attribute_Range_Ensemble:NoRange');
assert_error(@() Krige_Attribute_Cross_Validate(xs,ys,vs,'NotAnOption',1), ...
    'Krige_Attribute_Cross_Validate:UnknownOption');
ek = Krige_Attribute_Range_Ensemble(xs,ys,7*ones(size(vs)),xt,yt,'RangeDraws',20);
check('a constant field short-circuits to the constant with zero variance', ...
    all(ek.prediction == 7) && all(ek.std == 0) && ek.range_assumed && ...
    numel(ek.in_hull) == numel(xt));
assert_error(@() Krige_Attribute_Range_Ensemble(xs,ys,vs,xt,yt, ...
    'RangeBoundsM',[50 10]),'Krige_Attribute_Range_Ensemble:BadBounds');
assert_error(@() Kriging_Range_Prior_Bounds([1 1 1],[2 2 2]), ...
    'Kriging_Range_Prior_Bounds:CoincidentSamples');
check('input guards raise their documented errors', true);

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_kriging_uncertainty: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function f = gaussian_field(x,y,range_m,sill,nugget)
%% Zero-mean Gaussian field with the pykrige exponential covariance
%% (decay length range/3) plus a nugget, via Cholesky.
d = hypot(x-x.',y-y.');
C = (sill-nugget)*exp(-d/(range_m/3)) + nugget*eye(numel(x));
f = chol(C,'lower')*randn(numel(x),1);
end

%% =======================================================================
function s = rng_state_after(fcn)
rng(99);
fcn();
s = rng;
end

%% =======================================================================
function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end
