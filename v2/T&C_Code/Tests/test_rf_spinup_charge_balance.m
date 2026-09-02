function test_rf_spinup_charge_balance()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_rf_spinup_charge_balance
%%
%% Unit checks on the Lian-style RF pool-reconstruction machinery
%% (Select_Topographic_Subsample, Train_Pool_Random_Forest) and the
%% charge-closure primitive ported from the lian-spinup-convergence branch
%% (ERW_Initial_H_From_Anion_Budget).  The full field reconstruction
%% (Reconstruct_Cell_Pools_RF) is validated by a real field spin-up run.
%%
%%   1  the subsample spans every topographic stratum
%%   2  the forest recovers a known topo->pool relationship within its CV
%%      error, and CV error is reported
%%   3  Initialize_BG_State_From_Observations reproduces the totals it was
%%      given (conserved-total roundtrip)
%%   4  ERW_Initial_H_From_Anion_Budget returns a pore-water H whose
%%      residual anion equals the sourced budget (charge balance closes)
%%   5  a larger sourced anion budget drives a lower pore-water pH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_rf_spinup_charge_balance ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
v2   = fileparts(code);
addpath(code, fullfile(code,'Parameters'), fullfile(v2,'spinup'), here);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% --- synthetic topography-attached grid -----------------------------
dem = synthetic_hillslope_dem('hollow','Extent',[80 80],'CellSize',4, ...
    'Slope',0.06,'Relief',3);
topo = grid_topology(dem);
topo = Attach_Topography_To_Topology(topo, dem);

%% 1. subsample spans the strata
[idx,strata] = Select_Topographic_Subsample(topo,[],0.15,'Bins',3);
covered = unique(strata(idx));
check('subsample covers every non-empty stratum', ...
    isequal(sort(covered(:)), sort(unique(strata(:)))), ...
    sprintf('%d of %d strata',numel(covered),numel(unique(strata))));
check('subsample is ~15%% of cells', ...
    numel(idx) >= 0.10*topo.N_cells && numel(idx) <= 0.30*topo.N_cells, ...
    sprintf('%d of %d',numel(idx),topo.N_cells));

%% 2. RF recovers a known relationship
if exist('TreeBagger','file') == 2
    z = topo.cell_z_m(:); tw = topo.cell_twi(:);
    %% "SOC" rises with elevation (temperature proxy) and TWI, plus noise
    rng(1);
    SOC = 4000 + 120*(z-min(z)) + 300*(tw-mean(tw)) + 200*randn(topo.N_cells,1);
    SON = SOC/12 + 10*randn(topo.N_cells,1);
    X = [z, topo.cell_slope(:), cosd(topo.cell_aspect_deg(:)), tw];
    tr = idx;
    model = Train_Pool_Random_Forest(X(tr,:), [SOC(tr) SON(tr)], ...
        {'elev','slope','aspcos','twi'}, {'SOC','SON'}, 'NumTrees',150);
    Yhat = model.predict(X);
    rel = mean(abs(Yhat(:,1)-SOC))/mean(SOC);
    check('RF recovers the synthetic SOC field (< 15%% mean error)', ...
        rel < 0.15, sprintf('mean rel err %.3f',rel));
    check('RF reports a k-fold CV error per target', ...
        all(isfinite(model.cv_rel_error)) && numel(model.cv_rel_error) == 2);
else
    fprintf('  SKIP  RF checks (Statistics and Machine Learning Toolbox absent)\n');
end

%% 3. Initialize_BG_State_From_Observations roundtrip
B0 = seed_79_pool();
targets = struct('POC_gCm2',5200,'MAOC_gCm2',3100,'DOC_gCm2',9, ...
    'total_N_gNm2',480,'NO3_gNm2',1.4,'NH4_gNm2',0.7);
B1 = Initialize_BG_State_From_Observations(B0, targets);
got = struct('POC_gCm2',sum(B1(9:10)),'MAOC_gCm2',B1(11), ...
    'DOC_gCm2',sum(B1(12:13)),'total_N_gNm2',sum(B1(23:34)), ...
    'NO3_gNm2',B1(32),'NH4_gNm2',B1(31));
f = fieldnames(targets); ok = true; worst = '';
for k = 1:numel(f)
    if abs(got.(f{k})-targets.(f{k})) >= 1e-6*max(1,targets.(f{k}))
        ok = false; worst = f{k};
    end
end
check('observed totals are reproduced on the reconstructed 79-pool state', ...
    ok, worst);

%% 4-5. charge-balance closure from the sourced anion budget
conv_mol = 1e6; Al_CEC_scale = 1e3;
k1 = 4.45e-7*conv_mol; k2 = 4.69e-11*conv_mol; k_w = 1e-14*conv_mol^2;
k_nh3 = 5.6e-10*conv_mol;
K1 = 10^-5*conv_mol; K2 = 10^-5.1*conv_mol; K3 = 10^-6.7*conv_mol; K4 = 10^-6.2*conv_mol;
K_CEC = [0.5 0.6 3.0 0.8 1.2];              % Gaines-Thomas selectivity
K_AlOH3 = 8.5;
Ca = 900; Mg = 250; Na = 120; K_aq = 90; NO3 = 60; PO4 = 1; NH4 = 5;
CO2_w = 40;
An_target = 300;                            % [umolc/L] sourced strong anion
[H1,~,Al_w1] = ERW_Initial_H_From_Anion_Budget(An_target, [], Ca, Mg, Na, ...
    K_aq, NO3, PO4, NH4, CO2_w, k1, k2, k_w, k_nh3, K1, K2, K3, K4, ...
    K_CEC, K_AlOH3, Al_CEC_scale, conv_mol);
An_check = ERW_Residual_An_From_Charge(H1, CO2_w, Ca, Mg, Na, K_aq, NH4, ...
    NO3, PO4, Al_w1, k1, k2, k_w, k_nh3, K1, K2, K3, K4, Al_CEC_scale);
check('closed pore-water H reproduces the sourced anion budget', ...
    abs(An_check - An_target) < 1e-4*max(1,An_target), ...
    sprintf('An(H) = %.4f vs target %.4f',An_check,An_target));
pH1 = -log10(H1/conv_mol);
check('closed pore-water pH is in the physical range', pH1 > 3 && pH1 < 9, ...
    sprintf('pH = %.2f',pH1));

H2 = ERW_Initial_H_From_Anion_Budget(An_target+400, [], Ca, Mg, Na, K_aq, ...
    NO3, PO4, NH4, CO2_w, k1, k2, k_w, k_nh3, K1, K2, K3, K4, ...
    K_CEC, K_AlOH3, Al_CEC_scale, conv_mol);
check('a larger sourced anion budget lowers the pore-water pH', ...
    -log10(H2/conv_mol) < pH1, ...
    sprintf('pH %.2f -> %.2f',pH1,-log10(H2/conv_mol)));

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_rf_spinup_charge_balance: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function B = seed_79_pool()
B = zeros(1,79);
B(9:10) = [3000 1200];       % POC (litter/particulate)
B(11)   = 2000;              % MAOC
B(12:13)= [4 3];             % DOC bacteria/fungi
B(18:19)= [60 40];           % MBC
B(23:34)= [10 8 6 40 30 20 15 12 1 0.8 5 3];   % organic + inorganic N pools
B(35:47)= [30 25 20 200 150 100 80 50 5 4 3 2 1];% P pools
B(60)=1200; B(61)=800; B(68)=300; B(69)=200; B(52)=90; B(53)=60; B(76)=400;
end

%% =======================================================================
function topo = grid_topology(dem)
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
        topo.outer_face_width_m(id(i,j)) = ((i==1)+(i==ny)+(j==1)+(j==nx))*h;
    end
end
topo.adj = sparse([ea;eb],[eb;ea],1,topo.N_cells,topo.N_cells) > 0;
topo.cell_region = repmat({'main'},topo.N_cells,1);
end
