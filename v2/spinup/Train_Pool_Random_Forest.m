function model = Train_Pool_Random_Forest(X, Y, feature_names, target_names, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train_Pool_Random_Forest
%%
%% Fit one random forest per conserved-pool target (Lian et al. 2025,
%% Fig. 3e) mapping topographic / soil / vegetation covariates onto the
%% steady-state pool totals of the spun-up subsample cells, and report
%% k-fold cross-validation error (Lian et al. target: < 10 % relative).
%%
%% Requires the Statistics and Machine Learning Toolbox (TreeBagger).
%%
%% INPUTS
%%   X             [n_train x n_feat]   covariates of the subsample cells
%%   Y             [n_train x n_target] steady-state conserved pool totals
%%   feature_names {1 x n_feat}   covariate names (provenance only)
%%   target_names  {1 x n_target} pool names, e.g. {'POC_gCm2',...}
%%
%% Name/value
%%   'NumTrees'   trees per forest (default 200)
%%   'KFold'      cross-validation folds (default 5)
%%   'MinLeaf'    minimum leaf size (default 3)
%%   'Seed'       RNG seed (default 42)
%%
%% OUTPUT struct
%%   .forests{t}          the TreeBagger for target t
%%   .predict(Xq)         function handle -> [n_q x n_target] predictions
%%                        (non-negative, as pool totals must be)
%%   .cv_rel_error(t)     k-fold mean |pred-obs| / mean(obs)
%%   .cv_r2(t)            k-fold R^2
%%   .feature_names .target_names .n_train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('TreeBagger','file') ~= 2
    error('Train_Pool_Random_Forest:NoToolbox', ...
        ['TreeBagger (Statistics and Machine Learning Toolbox) is required ' ...
         'for the RF pool reconstruction. Use Run_Field_Cell_Spinup with ' ...
         'Method=''per_cell'' instead.']);
end
opt = struct('NumTrees',200,'KFold',5,'MinLeaf',3,'Seed',42);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

X = double(X); Y = double(Y);
[n, n_feat] = size(X);
n_target = size(Y,2);
assert(size(Y,1) == n, 'X and Y must have the same row count.');
if n < 2*opt.KFold
    error('Train_Pool_Random_Forest:TooFewTrain', ...
        'Need at least %d training cells for %d-fold CV, got %d.', ...
        2*opt.KFold, opt.KFold, n);
end

%% Seed the global stream for TreeBagger's bootstrap draws, but hand the
%% caller's stream back on exit -- a seed left behind would make every
%% later stochastic step in the session restart from the same point.
prev_rng = rng; restore_rng = onCleanup(@() rng(prev_rng)); %#ok<NASGU>
rng(opt.Seed);
forests = cell(1,n_target);
cv_rel  = nan(1,n_target);
cv_r2   = nan(1,n_target);
folds = mod(randperm(n), opt.KFold) + 1;

for t = 1:n_target
    y = Y(:,t);
    %% k-fold CV
    yhat = nan(n,1);
    for f = 1:opt.KFold
        tr = folds ~= f; te = ~tr;
        if nnz(tr) < 2 || nnz(te) < 1; continue; end
        b = TreeBagger(opt.NumTrees, X(tr,:), y(tr), 'Method','regression', ...
            'MinLeafSize',opt.MinLeaf, 'OOBPrediction','off');
        yhat(te) = max(0, predict(b, X(te,:)));
    end
    ok = isfinite(yhat);
    if any(ok) && mean(y(ok)) ~= 0
        cv_rel(t) = mean(abs(yhat(ok)-y(ok))) / mean(abs(y(ok)));
        ss_res = sum((y(ok)-yhat(ok)).^2);
        ss_tot = sum((y(ok)-mean(y(ok))).^2);
        cv_r2(t) = 1 - ss_res/max(ss_tot,eps);
    end
    %% Final forest on all training data
    forests{t} = TreeBagger(opt.NumTrees, X, y, 'Method','regression', ...
        'MinLeafSize',opt.MinLeaf, 'OOBPrediction','on');
end

model = struct();
model.forests = forests;
model.feature_names = feature_names;
model.target_names = target_names;
model.n_train = n;
model.n_features = n_feat;
model.cv_rel_error = cv_rel;
model.cv_r2 = cv_r2;
model.predict = @(Xq) predict_all(forests, double(Xq));

fprintf('Train_Pool_Random_Forest: %d targets, %d train cells, %d trees\n', ...
    n_target, n, opt.NumTrees);
for t = 1:n_target
    fprintf('  %-16s  CV rel err %5.1f%%   R2 %5.2f\n', ...
        target_names{t}, 100*cv_rel(t), cv_r2(t));
end
worst = max(cv_rel);
if worst > 0.10
    warning('Train_Pool_Random_Forest:HighCVError', ...
        ['Worst CV relative error is %.1f%% (> 10%%). More subsample ' ...
         'cells or fewer targets may be needed.'], 100*worst);
end
end

%% =======================================================================
function P = predict_all(forests, Xq)
P = zeros(size(Xq,1), numel(forests));
for t = 1:numel(forests)
    P(:,t) = max(0, predict(forests{t}, Xq));
end
end
