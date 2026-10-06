function test_field_spinup_warm_start()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_field_spinup_warm_start
%%
%% The warm start of Run_Field_Cell_Spinup ('WarmStartStates'): a field
%% spun up once is spun up again from its own converged states once site
%% observations replace its provisional soil inputs.  The expensive
%% Run_Automatic_Site_Spinup and Run_Field_Lateral_Tracking are replaced by
%% stubs (functions of the same names in the current folder) that log what
%% they are asked to do and write validated states:
%%   1  a missing states file, a cell the file lacks, a state from another
%%      grid and an output folder holding the earlier states are refused
%%      before any work
%%   2  per_cell: a cell's warm run is seeded with its own earlier state,
%%      runs as a warm start (SPINUP_WARM_START) under the observed inputs,
%%      carries a provenance key unlike the cold run's, and is reused on a
%%      rerun
%%   3  flux_tracking: no class seed is spun up; the tracking starts from
%%      every cell's own earlier state; the replays run as warm starts with
%%      their recorded exchange; and the forest stage accepts the warm
%%      replays (one provenance label where a replay runs and where it is
%%      read)
%% Uses the Control 2 diagnostic grid (15 m, 50 % coverage rule), the real
%% DEM tiles, and the archived converged Kantola state for the pools.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== test_field_spinup_warm_start ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code, fullfile(code,'Parameters'), fullfile(v2,'spinup'), ...
    fullfile(v2,'interpolation'), fullfile(v2,'run'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
    function ok = raises(fn,id)
        ok = false;
        try; fn(); catch ME
            ok = strcmp(ME.identifier,id);
            if ~ok; fprintf('    raised %s: %s\n',ME.identifier,ME.message); end
        end
    end

shp = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp');
dem = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','DEM');
h = 15; cov = 0.5; feats = [1 2 3]; poly = 1;
setups = Build_Shapefile_Grid_Setups(shp,h,feats,cov); topo = setups(poly);
N = topo.N_cells;
tmp = tempname; mkdir(tmp);

%% the provisional inputs of the first spin-up, and the observed ones
T = table(topo.cell_id,topo.cell_xc_projected_m,topo.cell_yc_projected_m, ...
    repmat(6.1,N,1),nan(N,1),false(N,1), ...
    repmat(0.10,N,1),repmat(0.35,N,1),nan(N,1),nan(N,1),false(N,1), ...
    repmat(0.025,N,1),nan(N,1),false(N,1), ...
    repmat(8.6,N,1),nan(N,1),false(N,1), ...
    nan(N,1),nan(N,1),nan(N,1),false(N,1), ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction', ...
    'measured_sand_fraction','measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC', ...
    'measured_bs_Ca_fraction','measured_bs_Mg_fraction','measured_bs_K_fraction', ...
    'use_measured_base_saturation'});
agro0 = fullfile(tmp,'provisional.csv'); writetable(T,agro0);
T.measured_bulk_pH(:) = 5.6 + 0.4*(0:N-1).'/max(N-1,1); T.use_measured_bulk_pH(:) = true;
T.measured_CEC_cmolc_kg(:) = 6 + 6*(0:N-1).'/max(N-1,1); T.use_measured_CEC(:) = true;
agro1 = fullfile(tmp,'observed.csv'); writetable(T,agro1);
pc = {'FieldFeatureIndices',feats,'FieldPolygonIndex',poly,'GridCellSizeM',h, ...
    'MinCoverage',cov,'Method','per_cell'};
[~,CI0] = Run_Field_Cell_Spinup(shp,agro0,pc{:},'OutputFile',fullfile(tmp,'x.mat'),'DryRun',true);
[~,CI1] = Run_Field_Cell_Spinup(shp,agro1,pc{:},'OutputFile',fullfile(tmp,'x.mat'),'DryRun',true);

%% the earlier field: converged states under the provisional inputs, the
%% organic stocks varying across the field so the forest has a signal
K = load(fullfile(code,'Inputs','Data_Kantola_BG.mat'));
B0 = K.B_out(end,:);
    function st = earlier_state(c)
        st = struct();
        B = B0; B(9:47) = B0(9:47)*(0.9 + 0.2*(c-1)/max(N-1,1));
        st.B_out = B; st.O_out = K.O_out(end,:); st.NavlI_out = K.NavlI_out(end,:);
        st.NupI_L_out = K.NupI_L_out; st.RexmyI_out = K.RexmyI_out(end,:);
        st.cell_id = topo.cell_id(c);
        st.cell_agronomic_identity = string(CI0(c).agronomic_identity);
        st.spinup_converged = true; st.spinup_complete = true; st.spinup_final_pass = true;
        st.spinup_erw_handoff_evaluated = true; st.spinup_erw_handoff_ready = true;
        st.spinup_nutrient_budget_vers = 12;
        b = struct('closed',true,'residual',0,'tolerance',1);
        st.spinup_nutrient_budgets = struct('N',b,'P',b,'K',b,'Ca',b,'Mg',b,'Si',b);
        st.spinup_biogeochemical_soil_mass_kg_m2 = 1295*300/1000;
        st.spinup_water_loss_mean_l_m2_d = 0.5;
        st.spinup_provenance = sprintf('earlier-%d',topo.cell_id(c));
    end
earlier = arrayfun(@(c) earlier_state(c),(1:N).');
    function f = states_file(name,CellBGStates) %#ok<INUSD>
        f = fullfile(tmp,[name '.mat']); save(f,'CellBGStates','-v7.3');
    end
earlier_file = states_file('earlier',earlier);

%% --- the stubs, in the current folder (they take precedence over the path)
stub_dir = fullfile(tmp,'stub'); mkdir(stub_dir);
template_file = fullfile(stub_dir,'template_state.mat');
save(template_file,'-struct','K');
write_stubs(stub_dir);
log_file = fullfile(stub_dir,'calls.log');
old_dir = pwd; old_path = path;
cleanup = onCleanup(@() restore(old_dir,old_path,tmp)); %#ok<NASGU>
cd(stub_dir); rehash; clear Run_Automatic_Site_Spinup Run_Field_Lateral_Tracking
setenv('TC_STUB_LOG',log_file); setenv('TC_STUB_TEMPLATE',template_file);
setenv('TC_STUB_DIR',stub_dir);

%% --- 1. refusals -----------------------------------------------------------
c = 2;
cell_stage = @(outdir,varargin) Run_Field_Cell_Spinup(shp,agro1,pc{:}, ...
    'OutputFile',fullfile(outdir,'CellBGStates.mat'),'Stage','cell','CellIndex',c,varargin{:});
check('a missing warm-start file is refused', raises(@() cell_stage(fullfile(tmp,'r1'), ...
    'WarmStartStates',fullfile(tmp,'none.mat')),'Run_Field_Cell_Spinup:WarmStatesMissing'));
check('warm-start states lacking a cell are refused', raises(@() cell_stage(fullfile(tmp,'r2'), ...
    'WarmStartStates',states_file('short',earlier(1:N-1))),'Run_Field_Cell_Spinup:WarmStateMissing'));
moved = earlier;
moved(N).cell_agronomic_identity = regexprep(moved(N).cell_agronomic_identity,'x=[^;]*', ...
    sprintf('x=%.17g',topo.cell_xc_projected_m(N)+1));
check('a warm-start state from another grid is refused', raises(@() cell_stage(fullfile(tmp,'r3'), ...
    'WarmStartStates',states_file('moved',moved)),'Run_Field_Cell_Spinup:WarmStateGrid'));
check('a warm start into the folder of the states it starts from is refused', ...
    raises(@() cell_stage(tmp,'WarmStartStates',earlier_file), ...
    'Run_Field_Cell_Spinup:WarmStartSameDirectory'));
check('no spin-up ran for a refused warm start', ~isfile(log_file));

%% --- 2. per_cell -----------------------------------------------------------
warm_dir = fullfile(tmp,'warm'); cold_dir = fullfile(tmp,'cold');
cell_stage(warm_dir,'WarmStartStates',earlier_file);
L = read_calls(log_file);
check('the warm cell runs once, seeded, as a warm start, without a replay', ...
    height(L) == 1 && L.cell(1) == topo.cell_id(c) && L.seeded(1) && L.warm(1) && ...
    ~L.handoff_only(1) && ~L.replay(1));
seed = load(fullfile(warm_dir,'warm',sprintf('warm_%04d.mat',topo.cell_id(c))));
check('its seed is its own earlier state', isequal(seed.B_out,earlier(c).B_out) && ...
    isequal(seed.O_out,earlier(c).O_out) && strcmp(seed.spinup_provenance,earlier(c).spinup_provenance));
S = load(fullfile(warm_dir,sprintf('cellstate_%04d.mat',topo.cell_id(c))));
check('the warm state runs under the observed inputs', ...
    strcmp(char(S.st.cell_agronomic_identity),CI1(c).agronomic_identity) && ...
    ~strcmp(CI1(c).agronomic_identity,CI0(c).agronomic_identity) && S.method == "full");
cell_stage(cold_dir);
L = read_calls(log_file);
check('the cold run of the same cell is a different run', height(L) == 2 && ...
    ~L.seeded(2) && ~L.warm(2) && ~strcmp(L.provenance(1),L.provenance(2)));
cell_stage(warm_dir,'WarmStartStates',earlier_file);
check('a rerun of the warm cell reuses its state', height(read_calls(log_file)) == 2);

%% --- 3. flux_tracking ------------------------------------------------------
ft_dir = fullfile(tmp,'flux'); mkdir(ft_dir);
ft = {'FieldFeatureIndices',feats,'FieldPolygonIndex',poly,'GridCellSizeM',h, ...
    'MinCoverage',cov,'OutputFile',fullfile(ft_dir,'CellBGStates.mat'), ...
    'FieldDEM',dem,'RFOptions',{'KFold',3},'WarmStartStates',earlier_file};
plan = readtable(Run_Field_Cell_Spinup(shp,agro1,ft{:},'Stage','plan'),'TextType','string');
tracked = find(logical(plan.is_tracked)); untracked = find(~logical(plan.is_tracked));
n_calls = height(read_calls(log_file));
seed_out = Run_Field_Cell_Spinup(shp,agro1,ft{:},'Stage','class_seed','ClassIndex',1);
check('a warm start spins up no class seed', isempty(seed_out) && ...
    ~isfile(fullfile(ft_dir,'classseed_01.mat')) && height(read_calls(log_file)) == n_calls);
Run_Field_Cell_Spinup(shp,agro1,ft{:},'Stage','track');
G = load(fullfile(stub_dir,'tracking_seeds.mat'));
B_seeds = vertcat(G.seeds.B_out); B_earlier = vertcat(earlier.B_out);
check('the tracking starts from every cell''s own earlier state', ...
    numel(G.seeds) == N && isequal(B_seeds,B_earlier(:,:)) && ...
    isequal(vertcat(G.seeds.O_out),vertcat(earlier.O_out)));
for k = tracked(:).'
    Run_Field_Cell_Spinup(shp,agro1,ft{:},'Stage','replay','CellIndex',k);
end
L = read_calls(log_file); L = L(n_calls+1:end,:);
check('every tracked cell replays its exchange as a warm start', ...
    height(L) == numel(tracked) && all(L.seeded) && all(L.warm) && all(L.replay) && ...
    isequal(sort(L.cell),sort(topo.cell_id(tracked))));
if exist('TreeBagger','file') == 2
    k = untracked(1);
    ok = true;
    try
        Run_Field_Cell_Spinup(shp,agro1,ft{:},'Stage','reconstruct','CellIndex',k);
    catch ME
        ok = false; fprintf('    raised %s: %s\n',ME.identifier,ME.message);
    end
    R = load(fullfile(ft_dir,sprintf('cellstate_%04d.mat',topo.cell_id(k))));
    check('the forest stage accepts the warm replays and reconstructs a cell', ...
        ok && R.method == "rf");
else
    fprintf('  SKIP  forest stage: Statistics and Machine Learning Toolbox not available\n');
end

fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_field_spinup_warm_start: %d check(s) failed',n_fail); end
end

%% =======================================================================
function write_stubs(stub_dir)
spin = {
'function bg_file = Run_Automatic_Site_Spinup(S)'
'%% Test stub: logs the run and writes a validated state from its seed'
'%% (seeded run) or from the converged template (cold run).'
'seeded = isfield(S,''SPINUP_SEED_FILE'');'
'warm = isfield(S,''SPINUP_WARM_START'') && S.SPINUP_WARM_START;'
'handoff_only = isfield(S,''SPINUP_HANDOFF_ONLY'') && S.SPINUP_HANDOFF_ONLY;'
'replay = isfield(S,''LATERAL_REPLAY_FILE'');'
'fid = fopen(getenv(''TC_STUB_LOG''),''a'');'
'fprintf(fid,''%d,%d,%d,%d,%d,%s\n'',S.CELL_ID,seeded,warm,handoff_only,replay,S.SPINUP_PROVENANCE); fclose(fid);'
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
'st.spinup_provenance = S.SPINUP_PROVENANCE;'
'save(S.BG_OUTPUT_FILE,''-struct'',''st'');'
'bg_file = S.BG_OUTPUT_FILE;'
'end'};
track = {
'function out = Run_Field_Lateral_Tracking(topo,CellInputs,cell_ids,seeds,track_idx,output_directory,case_number) %#ok<INUSL>'
'%% Test stub: keeps the seeds it is given and writes the records, the'
'%% replay seeds and the handoff states a tracking run leaves.'
'save(fullfile(getenv(''TC_STUB_DIR''),''tracking_seeds.mat''),''seeds'');'
'track_dir = fullfile(output_directory,''track''); if ~isfolder(track_dir); mkdir(track_dir); end'
'key = sprintf(''stub-%.17g'',sum(arrayfun(@(s) sum(s.B_out),seeds)));'
'for c = track_idx(:).'''
'    record = struct(''cell_id'',cell_ids(c),''identity'',sprintf(''record-%d'',cell_ids(c)),''tracking_key'',key); %#ok<NASGU>'
'    save(fullfile(track_dir,sprintf(''record_%04d.mat'',cell_ids(c))),''record'');'
'    seed = seeds(c); seed.cell_id = cell_ids(c); seed.tracking_key = key;'
'    save(fullfile(track_dir,sprintf(''ringend_%04d.mat'',cell_ids(c))),''-struct'',''seed'');'
'end'
'handoff_states = seeds; key_saved = key; %#ok<NASGU>'
'out = struct(''key'',key,''handoff_file'',fullfile(track_dir,''handoff_states.mat''));'
'save(out.handoff_file,''handoff_states'',''cell_ids'',''key_saved'');'
'end'};
write_lines(fullfile(stub_dir,'Run_Automatic_Site_Spinup.m'),spin);
write_lines(fullfile(stub_dir,'Run_Field_Lateral_Tracking.m'),track);
end

function write_lines(file_name,lines)
fid = fopen(file_name,'w'); fprintf(fid,'%s\n',lines{:}); fclose(fid);
end

function L = read_calls(log_file)
%% The stub's call log: cell id, seeded, warm, handoff-only and replay flags,
%% provenance key.
L = readtable(log_file,'FileType','text','Delimiter',',','ReadVariableNames',false, ...
    'TextType','string','Format','%f%f%f%f%f%s');
L.Properties.VariableNames = {'cell','seeded','warm','handoff_only','replay','provenance'};
L.seeded = logical(L.seeded); L.warm = logical(L.warm);
L.handoff_only = logical(L.handoff_only); L.replay = logical(L.replay);
end

function restore(old_dir,old_path,tmp)
cd(old_dir); path(old_path); rehash; clear Run_Automatic_Site_Spinup Run_Field_Lateral_Tracking
setenv('TC_STUB_LOG',''); setenv('TC_STUB_TEMPLATE',''); setenv('TC_STUB_DIR','');
try; rmdir(tmp,'s'); catch; end
end
