function [idx, strata] = Select_Topographic_Subsample(topo, resolved, fraction, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select_Topographic_Subsample
%%
%% Choose the cells to spin up in full for the Lian et al. (2025)
%% random-forest pool reconstruction: a stratified sample that spans the
%% joint distribution of the topographic features (elevation, slope,
%% aspect, curvature, TWI) and the soil / vegetation classes, so the
%% forest is trained across the whole feature space rather than a random
%% cluster of it.
%%
%% INPUTS
%%   topo      topography-attached topology (Attach_Topography_To_Topology)
%%   resolved  per-cell agronomic inputs (Resolve_Cell_Agronomic_Inputs);
%%             its soil / vegetation columns become extra strata.  [] to
%%             stratify on topography alone.
%%   fraction  target share of cells to spin up in full (default 0.12,
%%             Lian et al. used ~0.10)
%%
%% Name/value
%%   'MinPerStratum'  at least this many cells per non-empty stratum
%%                    (default 1)
%%   'Bins'           quantile bins per continuous feature (default 3)
%%   'Seed'           RNG seed for reproducibility (default 42)
%%
%% OUTPUTS
%%   idx      column vector of selected cell indices (into 1..topo.N_cells)
%%   strata   [N_cells x 1] stratum id of every cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(topo,'topography_attached') || ~topo.topography_attached
    error('Select_Topographic_Subsample:NoTopography', ...
        'Attach a DEM to the topology first (Attach_Topography_To_Topology).');
end
if nargin < 3 || isempty(fraction); fraction = 0.12; end
opt = struct('MinPerStratum',1,'Bins',3,'Seed',42);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

N = topo.N_cells;
feat = [topo.cell_z_m(:), topo.cell_slope(:), ...
    cosd(topo.cell_aspect_deg(:)), topo.cell_curv_prof(:), topo.cell_twi(:)];
%% Quantile-bin each continuous feature.
codes = zeros(N, size(feat,2));
for f = 1:size(feat,2)
    v = feat(:,f);
    if range(v) < eps
        codes(:,f) = 1;
    else
        edges = quantile(v, linspace(0,1,opt.Bins+1));
        edges(1) = -inf; edges(end) = inf;
        [~,~,codes(:,f)] = histcounts(v, unique(edges));
    end
end
%% Soil / vegetation classes, if supplied.
if ~isempty(resolved) && istable(resolved)
    class_cols = intersect(lower({'soil_class','texture_class','pft', ...
        'vegetation','land_cover'}), lower(resolved.Properties.VariableNames));
    for c = 1:numel(class_cols)
        hit = find(strcmpi(resolved.Properties.VariableNames,class_cols{c}),1);
        col = resolved.(resolved.Properties.VariableNames{hit});
        [~,~,cc] = unique(string(col));
        codes = [codes, cc(:)]; %#ok<AGROW>
    end
end

[strata_u,~,strata] = unique(codes,'rows');
n_target = max(1, round(fraction*N));

%% Seed the local draw but hand the caller's stream back on exit: a seed
%% left behind would restart every later stochastic step in the session
%% from the same point.
prev_rng = rng; restore_rng = onCleanup(@() rng(prev_rng)); %#ok<NASGU>
rng(opt.Seed);

%% Two-stage allocation.  Stage 1 gives every non-empty stratum its
%% guaranteed MinPerStratum cells; those are never trimmed, so the
%% forest is always trained across the whole feature space.  Stage 2
%% distributes what is left of the budget in proportion to stratum size
%% (largest remainder), so the sample also tracks the field's actual
%% composition.  Trimming at random -- the previous behaviour -- could
%% drop whole strata, leaving their cells both extrapolated by the
%% forest and seeded from a neighbour in a different stratum.
n_strata = size(strata_u,1);
members = arrayfun(@(s) find(strata == s), 1:n_strata, 'UniformOutput',false);
sizes = cellfun(@numel, members).';
guaranteed = min(sizes, opt.MinPerStratum);
if sum(guaranteed) > n_target
    n_target = sum(guaranteed);
    fprintf(['Select_Topographic_Subsample: %d strata x %d minimum ' ...
        'exceeds the %.0f%% budget; raising the sample to %d cells so no ' ...
        'stratum is dropped\n'], n_strata, opt.MinPerStratum, ...
        100*fraction, n_target);
end
spare = sizes - guaranteed;
extra = zeros(n_strata,1);
budget = n_target - sum(guaranteed);
if budget > 0 && sum(spare) > 0
    share = budget*spare/sum(spare);
    extra = min(spare, floor(share));
    [~,order] = sort(share-floor(share),'descend');
    left = budget - sum(extra);
    for k = order(:).'
        if left <= 0; break; end
        if extra(k) < spare(k); extra(k) = extra(k)+1; left = left-1; end
    end
end
idx = cell(n_strata,1);
for s = 1:n_strata
    take = guaranteed(s) + extra(s);
    idx{s} = members{s}(randperm(sizes(s), take));
end
idx = sort(unique(vertcat(idx{:})));

fprintf(['Select_Topographic_Subsample: %d of %d cells (%.1f%%) across ' ...
    '%d strata\n'], numel(idx), N, 100*numel(idx)/N, size(strata_u,1));
end
