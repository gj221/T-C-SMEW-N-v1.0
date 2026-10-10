function varargout = test_field_transient(n_days,workers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_field_transient
%%
%% The field transient (Run_Field_Transient) and the segment carry it rests
%% on (Field_Transient_Carry, Field_Transient_Restore):
%%   1  an ERW column carried to a longer segment keeps its last two
%%      days' ten-minute rows as its first, every later row at the array's
%%      value at the end, its last three daily rows and the CO2 stamp before
%%      the new axis; an array grown row by row keeps its rows to the
%%      anchor, and the CO2 record laid on the axis by the first ERW day is
%%      laid on the new one; a column with particle classes (rock), a segment
%%      that ran its final flush, a segment that does not start on the
%%      carried days and one that holds only them are refused, and so is a
%%      time-indexed driver array the carry does not list
%%   2  four Hamden cells (two map units' pH), each handed off from a short
%%      Hamden class seed (Tests/Fixtures/hamden_probe_class_seed.mat) by
%%      its own ERW handoff gate (SPINUP_HANDOFF_ONLY), run n_days (5) in
%%      2-D with the DEM when it is present: in one segment, in segments of
%%      2, 2 and 1 days and in segments of 3 and 2 days the records and the
%%      end state agree bit for bit; a CEC scenario runs every cell at its
%%      CEC from the same states, keeping their exchanger fractions
%%   3  a rerun reuses every segment and writes the same result; a
%%      finished segment of other code is refused
%%   4  the records: one per ERW day plus the start and the end of the
%%      window, at the right instants; pH, base saturation and the
%%      exchanger's fractions in range; leaching and drainage non-negative
%%   5  a state handed off at another instant, a subset naming a cell the
%%      grid lacks, and a site row with rock are refused
%% It runs for about an hour on a workstation (the handoffs and four short
%% 2-D runs of the four cells).
%%
%% Usage
%%   test_field_transient        % 5 days
%%   test_field_transient(6)
%%   test_field_transient(5,2)   % the runs on 2 thread workers
%%   out = test_field_transient(...)   % the runs' results and end states,
%%                                     % to compare one worker count's with
%%                                     % another's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(n_days); n_days = 5; end
if nargin < 2 || isempty(workers); workers = 0; end
if n_days < 5
    error('test_field_transient:BadDays','Segments of 2, 2 and 1 days need at least 5 days.');
end
fprintf('\n=== test_field_transient (%d days) ===\n',n_days);
here = fileparts(mfilename('fullpath'));
code = fileparts(here); v2 = fileparts(code);
addpath(code,fullfile(code,'Parameters'),fullfile(code,'Inputs'),fullfile(v2,'spinup'), ...
    fullfile(v2,'run'));
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
old_dir = pwd; restore_dir = onCleanup(@() cd(old_dir));
tmp = tempname; mkdir(tmp); cleanup = onCleanup(@() rmdir(tmp,'s'));
old_force = getenv('TC_CKPT_FORCE_RESUME');
setenv('TC_CKPT_FORCE_RESUME','');
restore_env = onCleanup(@() setenv('TC_CKPT_FORCE_RESUME',old_force));

%% 1  the ERW column carry ---------------------------------------------
%% A 3-day segment hands its last two days to a 4-day one.
NN = 72; NNd = 4; res = 6; n_t = NN*res; keep = 48*res;
E = struct('init',true,'init_done_j2',true,'res',res,'n_t',n_t,'dtw',1/(24*res), ...
    't_end',NN/24,'dm_end',(NN-23)*res,'hist',struct('names',{{}}));
E.t = E.dtw:E.dtw:E.t_end; E.th = (0:NN)'; E.td = E.th(1:24:end); E.tm = ((0:n_t-1)')/res;
Ca = 410 + (1:NN)'/100;
E.CO2_h = [Ca(1); Ca(:)]/22.41;
E.bg = struct('pH',(1:n_t)','f_CEC',[(1:n_t)' 2*(1:n_t)'],'pH_bulk',nan(n_t,1));
E.bg.pH(E.dm_end+1:end) = 0;                % rows the day did not reach
E.bg.pH_bulk(6:6:E.dm_end) = 5;
E.bg.K_CEC = (1:E.dm_end)' + (1:5);         % grown row by row to the anchor
E.bg.CO2_atm = nan(n_t,1);                  % laid on the axis by the first day
E.const = struct('T_ref',298.15*ones(n_t,1),'K1',1e-6);
E.TC = struct('s',(1:NNd)','Tw',ones(NN,1));
Ct = erw_carry(E,NN,NNd);
NN2 = 96; n_t2 = NN2*res;                   % the 4-day segment after it
Ca_next = 411 + (1:NN2)'/100;
E2 = Field_Transient_Carry('erw',Ct,NN2,5,Ca_next);
check('the carried days'' ten-minute rows become the first rows', ...
    isequal(E2.bg.pH(1:keep),E.bg.pH(n_t-keep+1:n_t)) && ...
    isequal(E2.bg.f_CEC(1:keep,:),E.bg.f_CEC(n_t-keep+1:n_t,:)));
check('every later row holds the array''s value at the end (its preallocation)', ...
    all(E2.bg.pH(keep+1:end) == 0) && all(isnan(E2.bg.pH_bulk(keep+1:end))) && ...
    all(E2.const.T_ref == 298.15) && size(E2.bg.f_CEC,1) == n_t2 && E2.const.K1 == 1e-6);
check('the ERW day resumes from the carried anchor, on the new axes', ...
    E2.dm_end == 25*res && E2.n_t == n_t2 && isequal(E2.tm,((0:n_t2-1)')/res) && ...
    isequal(E2.th,(0:NN2)') && numel(E2.t) == n_t2 && E2.t_end == 4);
check('the last three daily rows lead the daily series, zeros after', ...
    isequal(E2.TC.s,[2;3;4;0;0]) && isequal(E2.TC.Tw,zeros(NN2,1)));
check('the CO2 stamp at 00:00 of the carried days leads the new hourly axis', ...
    E2.CO2_h(1) == E.CO2_h(NN-47) && isequal(E2.CO2_h(2:end),Ca_next/22.41) && ...
    ~isfield(E2,'field_transient_carry'));
check('an array the ERW days grow row by row keeps its rows to the anchor', ...
    isequal(E2.bg.K_CEC,E.bg.K_CEC(n_t-keep+1:E.dm_end,:)));
check('the CO2 record the first ERW day lays on its axis is laid on the new one', ...
    numel(E2.bg.CO2_atm) == n_t2 && ...
    max(abs(E2.bg.CO2_atm(1:res:end) - E2.CO2_h(1:NN2)),[],'all') < 1e-12);
rock = E; rock.hist.names = {'psd'};
check('a column with particle classes (rock) is refused', ...
    raises(@() erw_carry(rock,NN,NNd),'Field_Transient_Carry:RockColumn'));
flushed = E; flushed.dm_end = NN*res;
check('a segment that ran its final flush is refused', ...
    raises(@() erw_carry(flushed,NN,NNd),'Field_Transient_Carry:FinalFlush'));
carry = struct('overlap_hours',48,'hour_dn',datenum(2020,1,2)+(0:47)'/24,'N_col',1);
check('a segment that does not start on the carried days is refused', ...
    raises(@() Field_Transient_Carry('check',carry,datenum(2020,1,3)+(0:71)'/24,1), ...
    'Field_Transient_Carry:Misaligned'));
check('a segment of the carried days alone is refused', ...
    raises(@() Field_Transient_Carry('check',carry,datenum(2020,1,2)+(0:47)'/24,1), ...
    'Field_Transient_Carry:ShortSegment'));
names = Field_Transient_Carry('names');
check('the carry lists every hourly and daily array once', ...
    numel(unique(names.hourly)) == numel(names.hourly) && ...
    numel(unique(names.daily)) == numel(names.daily) && ...
    isempty(intersect(names.hourly,names.daily)) && ...
    all(ismember({'O','V','Tdp','SWE','Lk','T_H','Ci_sunH'},names.hourly)) && ...
    all(ismember({'P','B_H','Bfac_dayH','PARI_L','NavlI','LEAK_Ca','PH_bulk_d'},names.daily)));
check('a time-indexed driver array the carry does not list is refused', ...
    raises(@() carry_with_stray_array(names),'Field_Transient_Carry:UncarriedArray'));

%% 2  four Hamden cells, one segment and several ------------------------
shp = fullfile(v2,'Spatial inputs','Hamden Blue Trail','Hamden_Blue_Trail_Parcel.shp');
dem = fullfile(v2,'Spatial inputs','Hamden Blue Trail','DEM');
if ~isfolder(dem) || isempty(dir(fullfile(dem,'*.tif')))
    fprintf('  (no DEM window: the cells run as a flat field)\n');
    dem = '';
end
grid_dir = fullfile(code,'Inputs','Field_Spinup','Hamden');
agro = fullfile(grid_dir,'Cell_Agronomic_Hamden_Blue_Trail_Parcel_polygon_01.csv');
M = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_MOD_columns.csv'));
Ecol = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_ERW_columns.csv'));
grid = {'GridCellSizeM',5,'MinCoverage','majority'};
[~,CI] = Run_Field_Cell_Spinup(shp,agro,grid{:},'MODParamColumns',M,'ERWParamColumns',Ecol, ...
    'OutputFile',fullfile(tmp,'unused.mat'),'DryRun',true);
%% A cell of the Cheshire unit (pH 5.3) with a Wilbraham neighbour (5.5)
%% and two more neighbours: a field of four cells joined by their faces.
setups = Build_Shapefile_Grid_Setups(shp,5,[],'majority'); topo = setups(1);
pH = arrayfun(@(ci) ci.MOD_PARAM.PH,CI);
ids = [];
for a = 1:topo.N_cells
    nb = find(topo.adj(a,:));
    if numel(nb) < 3 || pH(a) ~= 5.3; continue; end
    other = nb(pH(nb) == 5.5);
    if isempty(other); continue; end
    rest = setdiff(nb,other(1),'stable');
    ids = sort(topo.cell_id([a other(1) rest(1:2)])).';
    break
end
check('four cells of both map units, joined by their faces', ...
    numel(ids) == 4 && numel(unique(pH(ismember(topo.cell_id,ids)))) == 2, mat2str(ids));
seed0 = load(fullfile(here,'Fixtures','hamden_probe_class_seed.mat'));
states = cell(numel(ids),1);
gate_ok = true;
for k = 1:numel(ids)
    c = find(topo.cell_id == ids(k));
    seed = seed0;
    seed.cell_id = ids(k);
    seed.cell_agronomic_identity = CI(c).agronomic_identity;
    seed_file = fullfile(tmp,sprintf('seed_%04d.mat',ids(k)));
    save(seed_file,'-struct','seed','-v7.3');
    out_file = fullfile(tmp,sprintf('cell_%04d_BG.mat',ids(k)));
    settings = struct('NN_case_study_ERW_MOD',5,'RUN_MODE','spinup','RUN_CHEMISTRY','erw', ...
        'RUN_DOMAIN','single','SAVE_RESULTS',0,'MOD_PARAM_override',CI(c).MOD_PARAM, ...
        'ERW_PARAM_override',CI(c).ERW_PARAM,'CELL_ID',ids(k), ...
        'CELL_AGRONOMIC_IDENTITY',CI(c).agronomic_identity,'BG_OUTPUT_FILE',out_file, ...
        'SPINUP_SEED_FILE',seed_file,'SPINUP_HANDOFF_ONLY',true, ...
        'SPINUP_PROVENANCE','test_field_transient');
    log_text = evalc('Run_Automatic_Site_Spinup(settings);');
    gate_ok = gate_ok && contains(log_text,'ERW handoff readiness: pass');
    states{k} = load_bg_state(out_file);
end
check('every cell''s own ERW handoff passes its gate',gate_ok);
state_file = fullfile(tmp,'CellBGStates.mat');
CellBGStates = harmonise(states); %#ok<NASGU>
save(state_file,'CellBGStates','-v7.3');
cd(old_dir);
common = [grid {'MODParamColumns',M,'ERWParamColumns',Ecol,'FieldDEM',dem, ...
    'CellIds',ids,'NDays',n_days,'CheckpointMinutes',0,'ParallelWorkers',workers}];
t0 = tic;
f_one = Run_Field_Transient(shp,agro,state_file,common{:},'SegmentDays',n_days, ...
    'OutputDirectory',fullfile(tmp,'one'));
t_one = toc(t0);
f_two = Run_Field_Transient(shp,agro,state_file,common{:},'SegmentDays',2, ...
    'OutputDirectory',fullfile(tmp,'two'));
f_three = Run_Field_Transient(shp,agro,state_file,common{:},'SegmentDays',3, ...
    'OutputDirectory',fullfile(tmp,'three'));
R1 = load(f_one); R1 = R1.result;
R2 = load(f_two); R2 = R2.result;
R3 = load(f_three); R3 = R3.result;
fprintf('  one segment: %.1f min for %d cells x %d days\n',t_one/60,numel(ids),n_days);
fields = [R1.state_names R1.flux_names {'time','hours','precipitation_mm'}];
differs = @(R) fields(~cellfun(@(f) isequaln(R1.(f),R.(f)),fields));
check(sprintf(['segments of 2 days, the middle one carried in and handed on, reproduce ' ...
    'the single run bit for bit (%d records)'],numel(R1.time)), ...
    isempty(differs(R2)),strjoin(differs(R2),', '));
check('segments of 3 and 2 days reproduce it too',isempty(differs(R3)),strjoin(differs(R3),', '));
n2 = ceil(n_days/2); n3 = ceil(n_days/3);
S1 = load(fullfile(tmp,'one','segment_001.mat'),'info');
S2 = load(fullfile(tmp,'two',sprintf('segment_%03d.mat',n2)),'info');
S3 = load(fullfile(tmp,'three',sprintf('segment_%03d.mat',n3)),'info');
parts_differing = unique([part_differences(S1,S2) part_differences(S1,S3)]);
check('the end state is the same, whole: pools, soil water and heat, plants, ERW', ...
    strcmp(S1.info.end_state_digest,S2.info.end_state_digest) && ...
    strcmp(S1.info.end_state_digest,S3.info.end_state_digest), ...
    strjoin(parts_differing(:).',', '));

%% A CEC scenario: every cell at one CEC, from the same states.
f_cec = Run_Field_Transient(shp,agro,state_file,common{:},'NDays',2,'SegmentDays',2, ...
    'CEC',6.4,'OutputDirectory',fullfile(tmp,'cec'));
Rc = load(f_cec); Rc = Rc.result;
check('a CEC scenario runs every cell at its CEC, from the same states, into its own result', ...
    contains(f_cec,'_CEC_6p4') && all(Rc.cells.CEC_cmolc_kg == 6.4) && ...
    max(abs(Rc.CEC_umolc_m2(1,:)./R1.CEC_umolc_m2(1,:) - 6.4/10.8)) < 1e-9);
check('it keeps the handoff''s exchanger fractions: the same base saturation at the start', ...
    max(abs(Rc.base_saturation(1,:) - R1.base_saturation(1,:))) < 1e-3, ...
    mat2str([Rc.base_saturation(1,:); R1.base_saturation(1,:)],5));

%% 3  rerun and resume ---------------------------------------------------
f_again = '';    % declared: evalc cannot add a variable to this workspace
log_text = evalc(['f_again = Run_Field_Transient(shp,agro,state_file,common{:},' ...
    '''SegmentDays'',2,''OutputDirectory'',fullfile(tmp,''two''));']);
Ra = load(f_again); Ra = Ra.result;
check('a rerun reuses every finished segment and writes the same records', ...
    count(log_text,'reusing') == n2 && ~contains(log_text,'=== 2D driver') && ...
    all(cellfun(@(f) isequaln(R2.(f),Ra.(f)),fields)));
seg_file = fullfile(tmp,'two','segment_001.mat');
P = load(seg_file);
P.code_identity = 'other code';
save(seg_file,'-struct','P','-v7.3');
check('a finished segment of other code is refused', raises(@() Run_Field_Transient( ...
    shp,agro,state_file,common{:},'SegmentDays',2,'OutputDirectory',fullfile(tmp,'two')), ...
    'Run_Field_Transient:CodeChanged'));
check('a segment shorter than the two days the next repeats is refused', raises(@() ...
    Run_Field_Transient(shp,agro,state_file,common{:},'SegmentDays',1, ...
    'OutputDirectory',fullfile(tmp,'short')),'Run_Field_Transient:BadOptions'));

%% 4  the records ----------------------------------------------------------
t = R1.time;
start = datenum(char(M.transient_start),'yyyy-mm-dd');
expected = [start; start + 1/24; start + (1:n_days-1)' + 1/24; start + n_days];
check('one record at the start, at 01:00 of every day and at the end of the window', ...
    numel(t) == n_days + 2 && max(abs(t - expected)) < 1e-9 && ...
    isequal(R1.hours,[0; 1; 24*ones(n_days-1,1); 23]), mat2str(R1.hours.'));
pH_ok = all(R1.pH_bulk(:) > 3 & R1.pH_bulk(:) < 8) && all(R1.pH_pore(:) > 3 & R1.pH_pore(:) < 8);
f = cat(3,R1.f_Ca,R1.f_Mg,R1.f_K,R1.f_Na,R1.f_Al,R1.f_H);
check('pH in range, the exchanger''s fractions in [0,1] and summing to one', ...
    pH_ok && all(f(:) >= 0 & f(:) <= 1) && max(abs(sum(f,3)-1),[],'all') < 1e-6 && ...
    max(abs(R1.base_saturation - sum(f(:,:,1:4),3)),[],'all') < 1e-12);
check('the first record is each cell''s handoff', ...
    max(abs(R1.pH_bulk(1,:) - arrayfun(@(s) s.spinup_erw_handoff_diagnostic.pH_bulk_initial, ...
    harmonise(states).'))) < 1e-6);
flux = [R1.leach_Ca_g_m2(2:end-1,:); R1.leach_Mg_g_m2(2:end-1,:); R1.leach_K_g_m2(2:end-1,:); ...
    R1.leach_Al_umol_m2(2:end,:); R1.leach_Na_umol_m2(2:end,:); R1.drainage_mm(2:end,:)];
check('leaching and drainage are finite and non-negative; the end has no soil-biogeochemistry step', ...
    all(isfinite(flux(:))) && all(flux(:) >= 0) && all(isnan(R1.leach_Ca_g_m2(end,:))) && ...
    all(isnan(R1.leach_Ca_g_m2(1,:))));
check('field means over the cells'' areas', ...
    abs(R1.field.base_saturation(end) - R1.base_saturation(end,:)*R1.field_weights(:)) < 1e-12 && ...
    abs(sum(R1.field_weights) - 1) < 1e-12 && isfile(strrep(f_one,'.mat','_field_mean.csv')));

%% 5  refusals ---------------------------------------------------------------
M_late = M; M_late.transient_start = "2020-01-02";
check('a state handed off at another instant is refused', raises(@() Run_Field_Transient( ...
    shp,agro,state_file,grid{:},'MODParamColumns',M_late,'ERWParamColumns',Ecol, ...
    'CellIds',ids,'NDays',1,'OutputDirectory',fullfile(tmp,'late')), ...
    'Load_Field_Spinup_States:HandoffMismatch'));
check('a subset naming a cell the grid lacks is refused', raises(@() Run_Field_Transient( ...
    shp,agro,state_file,common{:},'CellIds',[ids 99999],'OutputDirectory',fullfile(tmp,'x')), ...
    'Run_Field_Transient:UnknownCell'));
E_rock = Ecol; E_rock.M_rock_in = 5;
check('a site row with rock is refused', raises(@() Run_Field_Transient( ...
    shp,agro,state_file,grid{:},'MODParamColumns',M,'ERWParamColumns',E_rock, ...
    'CellIds',ids,'NDays',1,'OutputDirectory',fullfile(tmp,'rock')), ...
    'Run_Field_Transient:Rock'));

fprintf('\n=== test_field_transient: %d passed, %d failed ===\n',n_pass,n_fail);
if n_fail > 0
    error('test_field_transient:Failed','%d check(s) failed.',n_fail);
end
if nargout > 0
    varargout{1} = struct('workers',workers,'one',R1,'two',R2,'three',R3,'cec',Rc, ...
        'end_state_digest',{{S1.info.end_state_digest,S2.info.end_state_digest, ...
        S3.info.end_state_digest}});
end
end

%% =======================================================================
function Ct = erw_carry(E,NN,NNd)
%% One ERW column through the capture, as a segment of NN hours hands it on.
W = struct('Date',datenum(2020,1,1)+(0:NN-1)'/24,'NN',NN,'NNd',NNd,'j',NNd, ...
    'ERW',{{E}},'lateral_diag',struct(),'jDay',1:NNd,'L_day',1:NNd);
names = Field_Transient_Carry('names');
for k = names.hourly; W.(k{1}) = zeros(NN,1); end
for k = names.daily;  W.(k{1}) = zeros(NNd,1); end
for k = names.column; W.(k{1}) = 0; end
C = Field_Transient_Carry('capture',W);
Ct = C.ERW{1};
end

function d = part_differences(A,B)
%% the parts of two end states that differ (Field_Transient_Carry 'parts')
a = A.info.end_state_parts; b = B.info.end_state_parts;
[common,ia,ib] = intersect(a(:,1),b(:,1),'stable');
d = [setxor(a(:,1),b(:,1)).' common(~strcmp(a(ia,2),b(ib,2))).'];
end

function carry_with_stray_array(names)
NN = 48; NNd = 3; N_col = 2;
W = struct('Date',datenum(2020,1,1)+(0:NN-1)'/24,'NN',NN,'NNd',NNd,'j',NNd, ...
    'ERW',{{struct(),struct()}},'lateral_diag',struct(),'jDay',1:NNd,'L_day',1:NNd);
for k = names.hourly; W.(k{1}) = zeros(NN,1,N_col); end
for k = names.daily;  W.(k{1}) = zeros(NNd,1,N_col); end
for k = names.column; W.(k{1}) = 0; end
W.whos = struct('name','Something_New','size',[NN 1 N_col],'class','double');
Field_Transient_Carry('capture',W);
end

function states = harmonise(cells)
names = {};
for c = 1:numel(cells); names = union(names,fieldnames(cells{c}),'stable'); end
for c = 1:numel(cells)
    for k = 1:numel(names)
        if ~isfield(cells{c},names{k}); cells{c}.(names{k}) = []; end
    end
    cells{c} = orderfields(cells{c},names);
end
states = vertcat(cells{:});
end
