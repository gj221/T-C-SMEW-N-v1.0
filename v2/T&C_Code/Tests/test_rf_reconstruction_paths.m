function test_rf_reconstruction_paths()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_rf_reconstruction_paths
%%
%% Reconstruct_Cell_Pools_RF (the forest stage of the flux-tracking field
%% spin-up) end to end on a synthetic hillslope, with the expensive
%% Run_Automatic_Site_Spinup replaced by a stub (a function of the same
%% name placed first on the path) that writes a validated state from the
%% seed it is given (the handoff-only run of a reconstructed cell), or from
%% a converged template on a cold start, and fails on demand.  The replayed
%% states are the archived, converged Kantola state (Data_Kantola_BG.mat)
%% scaled smoothly with elevation, so the seed mapping is exercised on a
%% real converged stoichiometry rather than on round numbers:
%%   1  every untracked cell is reconstructed (no silent fallback), and
%%      only through its handoff gate (SPINUP_HANDOFF_ONLY): the seed
%%      mapping accepts the template's own converged partition, takes the
%%      cell's own soil moisture from the tracking run, and records its
%%      template
%%   2  every returned state carries reconstruction_method and
%%      reconstruction_failed, and the array is one uniform struct
%%   3  a failed reconstruction falls back to a cold start, which runs
%%   4  a cell whose reconstruction and cold run both fail is returned
%%      flagged 'failed' with its OWN cell_id and unconverged / not
%%      handoff-ready flags, and the reconstruction of the other cells
%%      survives
%%   5  the reconstruction and the cold run of a cell carry different
%%      spin-up provenance keys (so neither resumes the other's
%%      checkpoints), a re-run of the same reconstruction carries the same
%%      key, and nothing on disk is deleted to tell them apart
%%   6  CellIndex reconstructs one cell, identical to the full run's
%%   7  an unknown option is rejected
%%   8  Select_Topographic_Subsample coarsens its bins when three bins
%%      leave more strata than the budget, and still covers every stratum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_rf_reconstruction_paths ===\n');
here = fileparts(mfilename('fullpath'));
helper = which('synthetic_hillslope_dem');          % the Tests folder, wherever this file is run from
if ~isempty(helper); here = fileparts(helper); end
code = fileparts(here); v2 = fileparts(code);
addpath(code, fullfile(code,'Parameters'), fullfile(v2,'spinup'), here);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
if exist('TreeBagger','file') ~= 2
    fprintf('  SKIP  Statistics and Machine Learning Toolbox not available\n');
    return
end

%% --- synthetic grid with topography --------------------------------------
dem = synthetic_hillslope_dem('hollow','Extent',[48 48],'CellSize',4, ...
    'Slope',0.06,'Relief',3);
topo = grid_topology(dem);
topo = Attach_Topography_To_Topology(topo, dem);
N = topo.N_cells;
topo.cell_id = (1:N).';
resolved = table((1:N).', 'VariableNames', {'cell_id'});
ident = arrayfun(@(c) sprintf('cell=%d;test', c), (1:N).', 'UniformOutput', false);
CellInputs = repmat(struct('cell_id',0,'MOD_PARAM',table(1,'VariableNames',{'Number'}), ...
    'ERW_PARAM',table(1,'VariableNames',{'Number'}),'agronomic_identity',''), N, 1);
for c = 1:N
    CellInputs(c).cell_id = c; CellInputs(c).agronomic_identity = ident{c};
end

%% --- converged template and subsample states ------------------------------
K = load(fullfile(code,'Inputs','Data_Kantola_BG.mat'));
B0 = K.B_out(end,:);
z = topo.cell_z_m(:); zn = (z - min(z))/max(range(z),eps);
[sub_idx,~] = Select_Topographic_Subsample(topo, [], 0.25);
sub_idx = sub_idx(:);
sub_states = repmat(mock_state(B0, sub_idx(1), ident{sub_idx(1)}, K), numel(sub_idx), 1);
for s = 1:numel(sub_idx)
    c = sub_idx(s);
    B = B0; f = 0.9 + 0.2*zn(c);               % organic stocks scale with elevation,
    B(9:47) = B0(9:47)*f;                      % stoichiometry preserved
    sub_states(s) = mock_state(B, c, ident{c}, K);
end

%% the tracking run's handoff states: every cell's own soil moisture
handoff = repmat(struct('O_out',K.O_out(end,:)),N,1);
for c = 1:N; handoff(c).O_out = K.O_out(end,:)*(0.95 + 0.1*zn(c)); end

%% --- the stub spin-up, first on the path -----------------------------------
stub_dir = tempname; mkdir(stub_dir);
out_dir = tempname; mkdir(out_dir);
template_file = fullfile(stub_dir,'template_state.mat');
save(template_file,'-struct','K');
write_stub(stub_dir);
old_dir = pwd; old_path = path;
cleanup = onCleanup(@() restore(old_dir, old_path, stub_dir, out_dir));
addpath(stub_dir,'-begin'); cd(stub_dir); rehash; clear Run_Automatic_Site_Spinup
log_file = fullfile(stub_dir,'calls.log');
setenv('TC_STUB_LOG',log_file); setenv('TC_STUB_TEMPLATE',template_file);
setenv('TC_STUB_FAIL_SEEDED',''); setenv('TC_STUB_FAIL_COLD','');
rf_opts = {'NumTrees',30,'KFold',3};

%% 1, 2, 5: nothing fails
rf_cells = setdiff((1:N).', sub_idx);
other_run = fullfile(out_dir, sprintf('cell_%04d_BG_checkpoint_p3_loop.mat', rf_cells(1)));
dummy = 1; save(other_run,'dummy'); %#ok<NASGU>
[S1,rep1] = Reconstruct_Cell_Pools_RF(topo,resolved,CellInputs,sub_states,sub_idx, ...
    handoff,'OutputDirectory',out_dir,'RFOptions',rf_opts);
check('every untracked cell is reconstructed, none silently falls back', ...
    rep1.n_rf == numel(rf_cells) && rep1.n_fallback == 0 && rep1.n_failed == 0, ...
    sprintf('replay %d, rf %d, fallback %d, failed %d', rep1.n_replay, rep1.n_rf, ...
    rep1.n_fallback, rep1.n_failed));
check('tracked cells keep their replayed state', rep1.n_replay == numel(sub_idx) && ...
    all(rep1.method_per_cell(sub_idx) == "replay"));
c_rf = rf_cells(1);
seed1 = load(fullfile(out_dir,sprintf('rf_seed_%04d.mat',c_rf)));
check('a reconstructed cell hands on its own soil moisture and names its template', ...
    isequal(seed1.O_out,handoff(c_rf).O_out) && ...
    ismember(seed1.reconstruction_template_cell_id,sub_idx) && ...
    seed1.reconstruction_method == "rf");
check('every state is stamped with its method and a failure flag', ...
    all(isfield(S1,{'reconstruction_method','reconstruction_failed'})) && ...
    isequal(string({S1.reconstruction_method}).', rep1.method_per_cell) && ...
    ~any([S1.reconstruction_failed]));
check('every state carries its own cell id', isequal([S1.cell_id].', (1:N).'));
check('no checkpoint is deleted to separate runs (provenance does it)', isfile(other_run));
L1 = read_calls(log_file);
check('every spin-up call carries a provenance key', ...
    height(L1) == numel(rf_cells) && all(strlength(L1.provenance) == 64));
check('a reconstructed cell runs only its handoff gate', ...
    all(L1.seeded == 1) && all(L1.handoff_only == 1));

%% 3, 4: a seeded failure falls back; a double failure is flagged
k_fb = rf_cells(2); k_dead = rf_cells(3);
setenv('TC_STUB_FAIL_SEEDED', sprintf('%d %d', k_fb, k_dead));
setenv('TC_STUB_FAIL_COLD', sprintf('%d', k_dead));
if isfile(log_file); delete(log_file); end
warning('off','Reconstruct_Cell_Pools_RF:CellFallback');
warning('off','Reconstruct_Cell_Pools_RF:CellFailed');
warning('off','Reconstruct_Cell_Pools_RF:UnreconstructedCells');
ok_run = true;
try
    [S2,rep2] = Reconstruct_Cell_Pools_RF(topo,resolved,CellInputs,sub_states,sub_idx, ...
        handoff,'OutputDirectory',out_dir,'RFOptions',rf_opts);
catch ME
    ok_run = false; fprintf('    %s\n', ME.message);
end
warning('on','Reconstruct_Cell_Pools_RF:CellFallback');
warning('on','Reconstruct_Cell_Pools_RF:CellFailed');
warning('on','Reconstruct_Cell_Pools_RF:UnreconstructedCells');
check('a double failure does not abort the reconstruction of the field', ok_run);
if ok_run
    L2 = read_calls(log_file);
    cold_for_fb = any(L2.cell == k_fb & L2.seeded == 0);
    seeded_key = L2.provenance(L2.cell == k_fb & L2.seeded == 1);
    cold_key = L2.provenance(L2.cell == k_fb & L2.seeded == 0);
    first_key = L1.provenance(L1.cell == k_fb);
    check('the reconstruction and the cold run of a cell carry different provenance', ...
        isscalar(seeded_key) && isscalar(cold_key) && seeded_key ~= cold_key);
    check('a re-run of the same reconstruction carries the same provenance', ...
        isscalar(first_key) && seeded_key == first_key);
    check('a failed reconstruction falls back to a cold start that runs', ...
        rep2.method_per_cell(k_fb) == "fallback" && cold_for_fb && ...
        S2(k_fb).reconstruction_method == "fallback");
    check('a double failure is flagged failed, with the cell''s own id', ...
        rep2.method_per_cell(k_dead) == "failed" && S2(k_dead).reconstruction_failed && ...
        S2(k_dead).cell_id == k_dead && ...
        strcmp(char(string(S2(k_dead).cell_agronomic_identity)), ident{k_dead}));
    check('a failed cell is marked unconverged and not handoff-ready', ...
        ~S2(k_dead).spinup_converged && ~S2(k_dead).spinup_erw_handoff_ready);
    check('the other cells are unaffected by the two failures', ...
        rep2.n_rf == numel(rf_cells) - 2 && rep2.n_failed == 1 && rep2.n_fallback == 1);
end
setenv('TC_STUB_FAIL_SEEDED',''); setenv('TC_STUB_FAIL_COLD','');

%% 6: CellIndex
k1 = rf_cells(4);
[S3,rep3] = Reconstruct_Cell_Pools_RF(topo,resolved,CellInputs,sub_states,sub_idx, ...
    handoff,'OutputDirectory',out_dir,'RFOptions',rf_opts,'CellIndex',k1);
check('CellIndex returns one reconstructed state for that cell, equal to the full run''s', ...
    numel(S3) == 1 && S3.cell_id == k1 && rep3.method_per_cell == "rf" && ...
    max(abs(S3.B_out(end,:) - S1(k1).B_out(end,:))) < 1e-9);

%% 7: options
ok = false;
try
    Reconstruct_Cell_Pools_RF(topo,resolved,CellInputs,sub_states,sub_idx, ...
        handoff,'OutputDirectory',out_dir,'ShortSpinupCap',50);
catch ME
    ok = strcmp(ME.identifier,'Reconstruct_Cell_Pools_RF:UnknownOption');
end
check('an unknown option is rejected', ok);

%% 8: bin coarsening
n_target = max(1, round(0.08*N));
[~,strata3] = Select_Topographic_Subsample(topo, [], 0.08, 'Bins', 3);
[idx_auto,strata_auto] = Select_Topographic_Subsample(topo, [], 0.08);
n3 = numel(unique(strata3)); na = numel(unique(strata_auto));
check('automatic binning coarsens when three bins exceed the budget', ...
    n3 <= n_target || na < n3, sprintf('3 bins: %d strata, auto: %d, budget %d', n3, na, n_target));
check('the coarsened subsample still covers every stratum', ...
    isequal(sort(unique(strata_auto(idx_auto))), sort(unique(strata_auto))));

fprintf('  %d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0; error('test_rf_reconstruction_paths: %d check(s) failed', n_fail); end
end

%% =======================================================================
function st = mock_state(B, cell_id, identity, K)
%% A state that load_bg_state and Validate_BG_State_For_Transient accept.
st = struct();
st.B_out = B; st.O_out = K.O_out; st.NavlI_out = K.NavlI_out;
st.NupI_L_out = K.NupI_L_out; st.RexmyI_out = K.RexmyI_out;
st.cell_id = cell_id; st.cell_agronomic_identity = identity;
st.spinup_converged = true; st.spinup_complete = true; st.spinup_final_pass = true;
st.spinup_erw_handoff_evaluated = true; st.spinup_erw_handoff_ready = true;
st.spinup_nutrient_budget_vers = 12;
b = struct('closed',true,'residual',0,'tolerance',1);
st.spinup_nutrient_budgets = struct('N',b,'P',b,'K',b,'Ca',b,'Mg',b,'Si',b);
st.spinup_biogeochemical_soil_mass_kg_m2 = 1295*300/1000;
st.spinup_water_loss_mean_l_m2_d = 0.5;
end

%% =======================================================================
function write_stub(stub_dir)
lines = {
'function bg_file = Run_Automatic_Site_Spinup(S)'
'%% Test stub: writes a validated state from the seed (seeded run) or from'
'%% the converged template (cold run); fails for the listed cell ids.'
'fail_seeded = str2num(getenv(''TC_STUB_FAIL_SEEDED'')); %#ok<ST2NM>'
'fail_cold   = str2num(getenv(''TC_STUB_FAIL_COLD''));   %#ok<ST2NM>'
'seeded = isfield(S,''SPINUP_SEED_FILE'');'
'handoff_only = isfield(S,''SPINUP_HANDOFF_ONLY'') && S.SPINUP_HANDOFF_ONLY;'
'fid = fopen(getenv(''TC_STUB_LOG''),''a''); fprintf(fid,''%d,%d,%d,%s\n'',S.CELL_ID,seeded,handoff_only,S.SPINUP_PROVENANCE); fclose(fid);'
'if seeded && ismember(S.CELL_ID,fail_seeded); error(''stub:seeded'',''stub seeded failure''); end'
'if ~seeded && ismember(S.CELL_ID,fail_cold); error(''stub:cold'',''stub cold failure''); end'
'if seeded; src = load(S.SPINUP_SEED_FILE); else; src = load(getenv(''TC_STUB_TEMPLATE'')); end'
'st = struct();'
'st.B_out = src.B_out(end,:); st.O_out = src.O_out(end,:); st.NavlI_out = src.NavlI_out(end,:);'
'st.NupI_L_out = src.NupI_L_out; st.RexmyI_out = src.RexmyI_out(end,:);'
'st.cell_id = S.CELL_ID; st.cell_agronomic_identity = char(string(S.CELL_AGRONOMIC_IDENTITY));'
'st.spinup_converged = true; st.spinup_complete = true; st.spinup_final_pass = true;'
'st.spinup_erw_handoff_evaluated = true; st.spinup_erw_handoff_ready = true;'
'st.spinup_nutrient_budget_vers = 12;'
'b = struct(''closed'',true,''residual'',0,''tolerance'',1);'
'st.spinup_nutrient_budgets = struct(''N'',b,''P'',b,''K'',b,''Ca'',b,''Mg'',b,''Si'',b);'
'st.spinup_biogeochemical_soil_mass_kg_m2 = 1295*300/1000;'
'st.spinup_water_loss_mean_l_m2_d = 0.5;'
'save(S.BG_OUTPUT_FILE,''-struct'',''st'');'
'bg_file = S.BG_OUTPUT_FILE;'
'end'};
fid = fopen(fullfile(stub_dir,'Run_Automatic_Site_Spinup.m'),'w');
fprintf(fid,'%s\n',lines{:}); fclose(fid);
end

function L = read_calls(log_file)
%% The stub's call log: cell id, seeded and handoff-only flags, provenance key.
L = readtable(log_file,'FileType','text','Delimiter',',','ReadVariableNames',false, ...
    'TextType','string','Format','%f%f%f%s');
L.Properties.VariableNames = {'cell','seeded','handoff_only','provenance'};
end

function restore(old_dir, old_path, stub_dir, out_dir)
cd(old_dir); path(old_path); rehash; clear Run_Automatic_Site_Spinup
setenv('TC_STUB_LOG',''); setenv('TC_STUB_TEMPLATE','');
setenv('TC_STUB_FAIL_SEEDED',''); setenv('TC_STUB_FAIL_COLD','');
try; rmdir(stub_dir,'s'); catch; end
try; rmdir(out_dir,'s'); catch; end
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
