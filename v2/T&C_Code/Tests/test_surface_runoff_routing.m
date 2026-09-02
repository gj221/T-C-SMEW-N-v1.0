function test_surface_runoff_routing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_surface_runoff_routing
%%
%%   1  no topography attached -> every cell keeps its own runoff
%%   2  MFD overland weights out of a cell sum to 1 where it has a
%%      downslope neighbour
%%   3  water conservation: delivered + exited == input
%%   4  on a tilted chain, runoff hops one cell downslope; only the
%%      outlet-adjacent share exits, the rest is delivered for re-infiltration
%%   5  solute rides the water, upwind from the contributing cell
%%   6  a convergent hollow concentrates runoff along the valley line
%%   7  an interior pit keeps its water (offered back, nothing exits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_surface_runoff_routing ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,here);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% ---------------------------------------------------------------------
%% 1. No topography -> pass-through
%% ---------------------------------------------------------------------
row = row_topology(5,10);
runoff = [3;1;4;1;5];
[q,~,exit_w] = Surface_Runoff_Routing(row,runoff,[],1);
check('no topography: q_deliver == runoff, nothing exits', ...
    isequal(q,runoff) && exit_w == 0);

%% ---------------------------------------------------------------------
%% 2-4. Tilted chain: 6% slope toward cell n
%% ---------------------------------------------------------------------
n = 8; h = 5;
chain = row_topology(n,h);
z = ((n:-1:1)'-1)*h*0.06;
chain = Attach_Topography_To_Topology(chain,struct('cell_z_m',z));

% weights out of an interior cell sum to 1
c = 4;
w_out = 0;
for e = 1:chain.N_edges
    if chain.edge_a(e) == c; w_out = w_out + chain.edge_overland_weight_ab(e); end
    if chain.edge_b(e) == c; w_out = w_out + chain.edge_overland_weight_ba(e); end
end
check('MFD weights out of an interior cell sum to 1', abs(w_out-1) < 1e-9, ...
    sprintf('sum = %.6f',w_out));

runoff = 2*ones(n,1);
[q,~,exit_w,~,diag] = Surface_Runoff_Routing(chain,runoff,[],1);
A = chain.cell_area_m2(:);
in_w  = sum(runoff/1000.*A);
out_w = sum(q/1000.*A) + exit_w;
check('water conservation (delivered + exited == input)', ...
    abs(in_w-out_w) < 1e-9*in_w, sprintf('in %.6g out %.6g',in_w,out_w));
% one hop: cell n (outlet) own runoff + cell n-1 -> cell n both exit
exit_expect = (runoff(n)+runoff(n-1))/1000*A(n);
check('only the outlet-adjacent share exits in one step', ...
    abs(exit_w - exit_expect) < 1e-9*max(1,exit_expect), ...
    sprintf('exit %.6g expected %.6g',exit_w,exit_expect));
check('interior cells receive their upslope neighbour''s runoff', ...
    all(q(2:n-1) > 0) && q(1) == 0 && q(n) == 0);
check('the outlet is the downslope end cell', chain.cell_is_outlet(n) && ...
    ~any(chain.cell_is_outlet(1:n-1)));
check('every interior edge carried water downslope (a->b positive)', ...
    all(diag.edge_flux_m3 > 0));

%% ---------------------------------------------------------------------
%% 5. Solute rides the water, upwind from the contributing cell
%% ---------------------------------------------------------------------
C = zeros(n,2);
C(1,1) = 40;                        % only the ridge cell has species 1
[~,C5,~,exit_m5] = Surface_Runoff_Routing(chain,runoff,C,1);
check('ridge cell 1 delivers its species-1 load to cell 2 at conc 40', ...
    abs(C5(2,1) - 40) < 1e-9 && all(C5([1 3:n],1) == 0), ...
    sprintf('C5(2,1) = %.6g',C5(2,1)));
check('species-1 does not reach the outlet in one hop from the ridge', ...
    exit_m5(1) == 0);
check('species 2 (absent everywhere) moves nothing', abs(exit_m5(2)) < 1e-12);

%% ---------------------------------------------------------------------
%% 6. Convergent hollow concentrates runoff along the valley
%% ---------------------------------------------------------------------
dem = synthetic_hillslope_dem('hollow','Extent',[40 40],'CellSize',2, ...
    'Slope',0.05,'Relief',2);
grid_topo = grid_topology(dem);
grid_topo = Attach_Topography_To_Topology(grid_topo,dem);
runoff_g = ones(grid_topo.N_cells,1);
[qg,~,exit_wg] = Surface_Runoff_Routing(grid_topo,runoff_g,[],1);
% valley cells (nearest the centre column) should deliver/collect more than
% ridge cells before exiting
xmid = 0.5*(min(grid_topo.cell_xc)+max(grid_topo.cell_xc));
near_valley = abs(grid_topo.cell_xc - xmid) < 3;
far_ridge   = abs(grid_topo.cell_xc - xmid) > 15;
check('hollow: valley cells collect more runoff than ridge cells', ...
    mean(qg(near_valley)) > mean(qg(far_ridge)) + 1e-9, ...
    sprintf('valley %.3f vs ridge %.3f mm',mean(qg(near_valley)),mean(qg(far_ridge))));
check('hollow: some runoff still leaves the domain', exit_wg > 0);

%% ---------------------------------------------------------------------
%% 7. Interior pit keeps its water
%% ---------------------------------------------------------------------
pit = row_topology(3,10);
pit = Attach_Topography_To_Topology(pit,struct('cell_z_m',[2;0;2]), ...
    'OutletCells',[]);
% cell 2 is a pit: lower than both neighbours, and (as a strip) it is on
% the perimeter, so derive_outlets will actually call it an outlet. Force
% the interior-pit case with an explicit non-outlet override:
pit.cell_is_outlet = false(3,1);
[qp,~,exit_wp,~,dp] = Surface_Runoff_Routing(pit,[1;1;1],[],1);
check('pit: nothing exits when no cell is an outlet', exit_wp == 0);
check('pit: all water ends up delivered (offered back)', ...
    abs(sum(qp/1000.*pit.cell_area_m2(:)) - sum([1;1;1]/1000.*pit.cell_area_m2(:))) < 1e-9);
check('pit: the low centre cell holds the most water', qp(2) == max(qp));

%% ---------------------------------------------------------------------
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_surface_runoff_routing: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function topo = row_topology(n,h)
topo.N_cells = n;
topo.N_edges = n-1;
topo.edge_a = (1:n-1)';
topo.edge_b = (2:n)';
topo.edge_face_width_m      = repmat(h,n-1,1);
topo.edge_face_distance_m   = repmat(h,n-1,1);
topo.edge_face_distance_a_m = repmat(h/2,n-1,1);
topo.edge_face_distance_b_m = repmat(h/2,n-1,1);
topo.edge_active = ones(n-1,1);
topo.cell_area_m2 = repmat(h^2,n,1);
topo.cell_xc = ((1:n)'-0.5)*h;
topo.cell_yc = repmat(h/2,n,1);
topo.outer_face_width_m = repmat(2*h,n,1);
topo.outer_face_width_m([1 n]) = 3*h;
topo.adj = sparse([topo.edge_a;topo.edge_b],[topo.edge_b;topo.edge_a],1,n,n) > 0;
topo.cell_region = repmat({'main'},n,1);
end

%% =======================================================================
function topo = grid_topology(dem)
%% Regular ny x nx cell graph over a synthetic DEM (test helper).
h = dem.options.CellSize;
nx = numel(dem.x); ny = numel(dem.y);
[C,Rr] = meshgrid(1:nx,1:ny);
id = reshape(1:nx*ny,ny,nx);
topo.N_cells = nx*ny;
topo.cell_xc = dem.x(C(:)).';
topo.cell_yc = dem.y(Rr(:)).';
topo.cell_area_m2 = repmat(h^2,nx*ny,1);
ea = []; eb = [];
for j = 1:nx
    for i = 1:ny
        if j < nx; ea(end+1,1) = id(i,j); eb(end+1,1) = id(i,j+1); end %#ok<AGROW>
        if i < ny; ea(end+1,1) = id(i,j); eb(end+1,1) = id(i+1,j); end %#ok<AGROW>
    end
end
topo.N_edges = numel(ea);
topo.edge_a = ea; topo.edge_b = eb;
topo.edge_face_width_m      = repmat(h,topo.N_edges,1);
topo.edge_face_distance_m   = repmat(h,topo.N_edges,1);
topo.edge_face_distance_a_m = repmat(h/2,topo.N_edges,1);
topo.edge_face_distance_b_m = repmat(h/2,topo.N_edges,1);
topo.edge_active = ones(topo.N_edges,1);
topo.outer_face_width_m = zeros(topo.N_cells,1);
for i = 1:ny
    for j = 1:nx
        e = (i==1)+(i==ny)+(j==1)+(j==nx);
        topo.outer_face_width_m(id(i,j)) = e*h;
    end
end
topo.adj = sparse([ea;eb],[eb;ea],1,topo.N_cells,topo.N_cells) > 0;
topo.cell_region = repmat({'main'},topo.N_cells,1);
end
