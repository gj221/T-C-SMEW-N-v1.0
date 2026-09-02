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
asp  = cosd(topo.cell_aspect_deg(:));         % 1 = south-facing (Lian et al.)

%% --- standardized gradient-method sensitivity (Lian et al. Eq. 1) ----
feat = [elev, twi, asp, curv];
featz = (feat - mean(feat,1)) ./ max(std(feat,0,1),eps);
importance = zeros(V,4);
for v = 1:V
    yv = pools_topo(:,v);
    keep = yv <= quantile(yv,0.99) & isfinite(yv);      % drop top-1% outliers
    yz = (yv(keep) - mean(yv(keep))) / max(std(yv(keep)),eps);
    Xz = featz(keep,:);
    %% full quadratic + interaction design (Lian et al. Eq. 1)
    D = [Xz, ...
         Xz(:,1).*Xz(:,2), Xz(:,1).*Xz(:,3), Xz(:,1).*Xz(:,4), ...
         Xz(:,2).*Xz(:,3), Xz(:,2).*Xz(:,4), Xz(:,3).*Xz(:,4), ...
         Xz.^2, ones(nnz(keep),1)];
    b = D \ yz;
    %% mean |partial derivative| w.r.t. each standardized feature
    dfdx = zeros(nnz(keep),4);
    dfdx(:,1) = b(1) + b(5)*Xz(:,2) + b(6)*Xz(:,3) + b(7)*Xz(:,4) + 2*b(14)*Xz(:,1);
    dfdx(:,2) = b(2) + b(5)*Xz(:,1) + b(8)*Xz(:,3) + b(9)*Xz(:,4) + 2*b(15)*Xz(:,2);
    dfdx(:,3) = b(3) + b(6)*Xz(:,1) + b(8)*Xz(:,2) + b(10)*Xz(:,4) + 2*b(16)*Xz(:,3);
    dfdx(:,4) = b(4) + b(7)*Xz(:,1) + b(9)*Xz(:,2) + b(10)*Xz(:,3) + 2*b(17)*Xz(:,4);
    m = mean(abs(dfdx),1);
    importance(v,:) = m / max(sum(m),eps);
end
sensitivity = array2table(importance, ...
    'VariableNames',{'elevation','TWI','aspect','curvature'}, ...
    'RowNames',names);

%% --- TWI relationship of the 2-D minus flat difference --------------
twi_slope = nan(V,1);
if have_flat
    for v = 1:V
        d = pools_topo(:,v) - pools_flat(:,v);
        p = polyfit(twi(isfinite(d)), d(isfinite(d)), 1);
        twi_slope(v) = p(1);
    end
end

%% --- figure --------------------------------------------------------
rows = V; cols = ternary(have_flat, 3, 2);
fig = figure('Color','w','Position',[60 60 360*cols 300*rows]);
tl = tiledlayout(fig, rows, cols, 'TileSpacing','compact','Padding','compact');
for v = 1:V
    nexttile(tl);
    scatter(x, y, 18, pools_topo(:,v), 'filled'); axis image; colorbar;
    title(sprintf('%s  (2-D topographic)', names{v}), 'FontWeight','normal');
    if v == 1; ylabel(sprintf('%s',opt.Tag),'FontWeight','bold'); end

    if have_flat
        nexttile(tl);
        d = pools_topo(:,v) - pools_flat(:,v);
        scatter(x, y, 18, d, 'filled'); axis image; colorbar;
        try; colormap(gca, redbluecmap); catch; colormap(gca, parula); end
        title(sprintf('%s  (2-D minus flat)', names{v}), 'FontWeight','normal');
    end

    nexttile(tl);
    if have_flat
        d = pools_topo(:,v) - pools_flat(:,v);
        scatter(twi, d, 14, curv, 'filled'); ylabel('2-D minus flat');
        hold on; xx = linspace(min(twi),max(twi),20);
        plot(xx, polyval(polyfit(twi(isfinite(d)),d(isfinite(d)),1), xx), 'k-');
        title(sprintf('slope %.3g per TWI unit', twi_slope(v)), 'FontWeight','normal');
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
