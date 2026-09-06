function out = Plot_Topographic_Redistribution(topo, pools_topo, pools_flat, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Topographic_Redistribution
%%
%% Lian et al. (2025) Fig. 7e-f / Fig. 8-9 analogue for a T&C-SMEW-N
%% field run: how topography-driven lateral transport reshapes the soil
%% carbon, alkalinity and nutrient pools, and which topographic feature
%% each pool is most sensitive to.
%%
%% INPUTS
%%   topo         topography-attached topology (Attach_Topography_To_Topology)
%%   pools_topo   [N_cells x N_var] long-term-mean pool values from the
%%                topographic 2-D run (elevation-head Darcy + overland
%%                routing on)
%%   pools_flat   [N_cells x N_var] the same variables from the flat run
%%                (no DEM attached) OR from per-cell 1-D runs.  [] to plot
%%                only the topographic field and its feature sensitivity.
%%
%% Name/value
%%   'VarNames'        {1 x N_var} labels (default var1..varN)
%%   'FigureDirectory' output folder (default v2/figures)
%%   'Tag'             filename tag (default 'run')
%%
%% OUTPUT struct
%%   .sensitivity   table: variable x {elevation TWI aspect curvature}
%%                  relative importance from the standardized gradient
%%                  method (Lian et al. Eq. 1-2)
%%   .twi_slope     [N_var x 1] OLS slope of (topo - flat) on TWI
%%   .figure_file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'));

opt = struct('VarNames',{{}},'FigureDirectory',fullfile(v2,'figures'), ...
    'Tag','run');
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if ~exist(opt.FigureDirectory,'dir'); mkdir(opt.FigureDirectory); end

N = topo.N_cells;
V = size(pools_topo,2);
names = opt.VarNames;
if isempty(names); names = arrayfun(@(k) sprintf('var%d',k), 1:V, ...
        'UniformOutput',false); end
have_flat = ~isempty(pools_flat) && isequal(size(pools_flat),size(pools_topo));

x = topo.cell_xc_projected_m(:); y = topo.cell_yc_projected_m(:);
twi = topo.cell_twi(:); curv = topo.cell_curv_prof(:);
elev = topo.cell_z_m(:);
%% aspect is degrees clockwise from north (downslope direction), so a
%% south-facing cell has aspect 180: negate the cosine to get the Lian et
%% al. encoding (+1 = south-facing, -1 = north-facing).
asp  = -cosd(topo.cell_aspect_deg(:));

%% --- standardized gradient-method sensitivity (Lian et al. Eq. 1) ----
feat = [elev, twi, asp, curv];
feat_sd = std(feat,0,1);
%% A feature with no variation across the domain carries no information
%% about the pools, and leaving it in makes the design matrix rank
%% deficient.  Drop it and report its importance as NaN -- "undetermined"
%% -- rather than 0, which would read as "this driver does not matter".
%% The elevation-injection topology path is the usual case: it reports a
%% single slope/aspect and zero curvature for every cell.
%% A feature that is NaN somewhere (TWI on the elevation-injection path)
%% is undefined, not constant; it is dropped the same way but said so.
fn = {'elevation','TWI','aspect','curvature'};
undefined = ~isfinite(feat_sd);
varying = isfinite(feat_sd) & feat_sd > 1e-9*max(1,max(abs(feat),[],1));
if any(undefined)
    fprintf('Plot_Topographic_Redistribution: %s undefined (NaN) for some cells; importance NaN\n', ...
        strjoin(fn(undefined),', '));
end
featz = (feat - mean(feat,1)) ./ max(feat_sd,eps);
featz(:,~varying) = 0;
feature_names = {'elevation','TWI','aspect','curvature'};
if ~all(varying)
    fprintf(['Plot_Topographic_Redistribution: %s constant across the ' ...
        'domain; importance reported as NaN\n'], ...
        strjoin(feature_names(~varying),', '));
end
importance = nan(V,4);
for v = 1:V
    yv = pools_topo(:,v);
    keep = yv <= quantile(yv,0.99) & isfinite(yv);      % drop top-1% outliers
    yz = (yv(keep) - mean(yv(keep))) / max(std(yv(keep)),eps);
    Xz = featz(keep,:);
    %% Solve on the varying columns only (a constant column would make the
    %% least-squares system rank deficient); the dropped terms contribute
    %% zero to every derivative below because their Xz column is zero.
    %% full quadratic + interaction design (Lian et al. Eq. 1):
    %% columns 1-4 linear, 5-10 the six pairwise interactions
    %% (12,13,14,23,24,34), 11-14 the squares, 15 the intercept.
    D = [Xz, ...
         Xz(:,1).*Xz(:,2), Xz(:,1).*Xz(:,3), Xz(:,1).*Xz(:,4), ...
         Xz(:,2).*Xz(:,3), Xz(:,2).*Xz(:,4), Xz(:,3).*Xz(:,4), ...
         Xz.^2, ones(nnz(keep),1)];
    %% Minimum-norm solution: a dropped (constant) feature leaves zero
    %% columns in D, which makes a backslash solve rank deficient and warn.
    %% lsqminnorm gives those terms exactly zero coefficients, which is
    %% what a predictor with no variation should contribute.
    b = lsqminnorm(D,yz);
    %% mean |partial derivative| w.r.t. each standardized feature
    dfdx = zeros(nnz(keep),4);
    dfdx(:,1) = b(1) + b(5)*Xz(:,2) + b(6)*Xz(:,3) + b(7)*Xz(:,4)  + 2*b(11)*Xz(:,1);
    dfdx(:,2) = b(2) + b(5)*Xz(:,1) + b(8)*Xz(:,3) + b(9)*Xz(:,4)  + 2*b(12)*Xz(:,2);
    dfdx(:,3) = b(3) + b(6)*Xz(:,1) + b(8)*Xz(:,2) + b(10)*Xz(:,4) + 2*b(13)*Xz(:,3);
    dfdx(:,4) = b(4) + b(7)*Xz(:,1) + b(9)*Xz(:,2) + b(10)*Xz(:,3) + 2*b(14)*Xz(:,4);
    %% Normalise over the features that actually vary, so the reported
    %% shares still sum to 1 across them.
    m = mean(abs(dfdx),1);
    m(~varying) = NaN;
    importance(v,:) = m / max(sum(m(varying)),eps);
end
sensitivity = array2table(importance, ...
    'VariableNames',{'elevation','TWI','aspect','curvature'}, ...
    'RowNames',names);

%% --- TWI relationship of the 2-D minus flat difference --------------
twi_slope = nan(V,1);
if have_flat
    for v = 1:V
        d = pools_topo(:,v) - pools_flat(:,v);
        ok = isfinite(d) & isfinite(twi);
        if nnz(ok) >= 2
            p = polyfit(twi(ok), d(ok), 1);
            twi_slope(v) = p(1);
        end
    end
end

%% --- figure --------------------------------------------------------
rows = V; cols = ternary(have_flat, 3, 2);
fig = figure('Color','w','Position',[60 60 360*cols 300*rows]);
try; fig.Theme = 'light'; catch; end       % batch sessions default to dark axes
tl = tiledlayout(fig, rows, cols, 'TileSpacing','compact','Padding','compact');
%% Local metres: projected coordinates (7.26e5, 4.07e6) would push an
%% axis exponent into the panel titles.
xl = x - min(x); yl = y - min(y);
for v = 1:V
    nexttile(tl);
    scatter(xl, yl, 18, pools_topo(:,v), 'filled'); axis image; colorbar;
    title(sprintf('%s  (2-D topographic)', names{v}), 'FontWeight','normal');
    xlabel('x [m]'); ylabel('y [m]');
    if v == 1; ylabel(sprintf('%s   y [m]',opt.Tag),'FontWeight','bold','Interpreter','none'); end

    if have_flat
        nexttile(tl);
        d = pools_topo(:,v) - pools_flat(:,v);
        scatter(xl, yl, 18, d, 'filled'); axis image; colorbar;
        try; colormap(gca, redbluecmap); catch; colormap(gca, parula); end
        title(sprintf('%s  (2-D minus flat)', names{v}), 'FontWeight','normal');
        xlabel('x [m]'); ylabel('y [m]');
    end

    nexttile(tl);
    if have_flat
        d = pools_topo(:,v) - pools_flat(:,v);
        scatter(twi, d, 14, curv, 'filled'); ylabel('2-D minus flat');
        ok = isfinite(d) & isfinite(twi);
        if nnz(ok) >= 2
            hold on; xx = linspace(min(twi(ok)),max(twi(ok)),20);
            plot(xx, polyval(polyfit(twi(ok),d(ok),1), xx), 'k-');
        end
        title(sprintf('slope %.3g per TWI unit', twi_slope(v)), 'FontWeight','normal');
        ax = gca; ax.YAxis.Exponent = 0;    % keep the exponent out of the title
    else
        scatter(twi, pools_topo(:,v), 14, curv, 'filled');
        ylabel(names{v}); title('vs TWI', 'FontWeight','normal');
    end
    xlabel('topographic wetness index'); cb = colorbar; cb.Label.String = 'profile curvature';
end
title(tl, sprintf(['Topographic redistribution - %s   ' ...
    '(mean feature importance: elev %.0f%%, TWI %.0f%%, aspect %.0f%%, curv %.0f%%)'], ...
    opt.Tag, 100*mean(importance,1)), 'Interpreter','none');

figure_file = fullfile(opt.FigureDirectory, ...
    sprintf('Topographic_Redistribution_%s.png', opt.Tag));
exportgraphics(fig, figure_file, 'Resolution',150);
savefig(fig, strrep(figure_file,'.png','.fig'));

out = struct('sensitivity',sensitivity,'twi_slope',twi_slope, ...
    'figure_file',figure_file);
fprintf('Plot_Topographic_Redistribution: wrote %s\n', figure_file);
disp(sensitivity);
end

function o = ternary(c,a,b); if c; o = a; else; o = b; end; end
