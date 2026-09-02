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

rng(opt.Seed);
idx = [];
per_stratum = max(opt.MinPerStratum, ceil(n_target/size(strata_u,1)));
for s = 1:size(strata_u,1)
    members = find(strata == s);
    take = min(numel(members), per_stratum);
    idx = [idx; members(randperm(numel(members),take))]; %#ok<AGROW>
end
%% Top up (or trim) to the target count.
if numel(idx) < n_target
    rest = setdiff((1:N).', idx);
    add  = rest(randperm(numel(rest), min(numel(rest), n_target-numel(idx))));
    idx  = [idx; add];
elseif numel(idx) > n_target
    idx = idx(randperm(numel(idx), n_target));
end
idx = sort(unique(idx));

fprintf(['Select_Topographic_Subsample: %d of %d cells (%.1f%%) across ' ...
    '%d strata\n'], numel(idx), N, 100*numel(idx)/N, size(strata_u,1));
end
