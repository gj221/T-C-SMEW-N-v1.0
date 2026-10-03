function test_lateral_flux_tracking(n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_lateral_flux_tracking
%%
%% The flux-tracking machinery of the field spin-up (Lian et al. 2025):
%% record a cell's lateral exchange in the 2-D field run, replay it in the
%% plot-scale column.
%%
%%   A  Topographic_Radiation_Projection: flat ground keeps the station
%%      values exactly; a south-facing slope gets more noon beam than a
%%      north-facing one.
%%   B  Lateral_Track_Record composes the hourly upwind exchange into the
%%      day's affine map exactly, splits a face's two directions layer by
%%      layer, and keeps an empty pool's export proportional to the water.
%%   C  Lateral_Replay_*: a spin-up year outside the record is refused, a
%%      seed-check window outside it runs with no exchange, the water bounds
%%      are kept with the refused part booked, and the solute ledger is
%%      exact (in - out = applied change).
%%   D  The decisive check: a small topographic 2-D run (ERW off, as in
%%      every spin-up phase, with a storm so run-on and lateral solute
%%      exchange are substantial) tracks two cells; each is then run as a
%%      plot-scale column replaying its record.  The replay must follow the
%%      2-D cell (hourly water contents, daily mobile pools) far more
%%      closely than the same column with the exchange left out.
%%
%% Uses the forcing bundled in the equivalence baseline
%%   v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(n_days); n_days = 3; end
fprintf('\n=== test_lateral_flux_tracking (%d days) ===\n', n_days);
here = fileparts(mfilename('fullpath'));
v2   = fileparts(fileparts(here));
code = fullfile(v2,'T&C_Code');
addpath(code, fullfile(code,'Inputs'), fullfile(code,'Parameters'), here);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end
    function ok = raises(fn,id)
        ok = false;
        try; fn(); catch ME; ok = strcmp(ME.identifier,id); if ~ok; disp(ME.message); end; end
    end

%% --- A. radiation projection --------------------------------------------
Dm = [2016 6 21 13]; Lat = 36.7; Lon = -78.5; dG = round(Lon/15);
[a1,a2,a3,a4,a5,a6] = Topographic_Radiation_Projection(500,300,120,80,250,90, ...
    Dm,dG,Lon,Lat,1,0,0,180,1);
check('flat ground keeps every station value exactly', ...
    isequal([a1 a2 a3 a4 a5 a6],[500 300 120 80 250 90]));
south = Topographic_Radiation_Projection(500,300,120,80,250,90,Dm,dG,Lon,Lat,1,0,0.15,180,0.98);
north = Topographic_Radiation_Projection(500,300,120,80,250,90,Dm,dG,Lon,Lat,1,0,0.15,0,0.98);
check('a south-facing slope gets more noon beam than a north-facing one', south > north, ...
    sprintf('south %.1f, north %.1f W m-2',south,north));

%% --- B. track composition ------------------------------------------------
%% Three cells in a row, two layers, two species (pool 31 and 32).
topo3 = struct('N_cells',3,'edge_a',[1;2],'edge_b',[2;3]);
species = {'NH4','NO3'}; pool_idx = {31,32};
lf = [0.5 1; 0.5 1; 0.5 1];
tr = Lateral_Track_Record('init',[],[2 3],3,4,3,2,species,pool_idx,lf,(0:3)'/24);
tr = Lateral_Track_Record('day',tr,2,1.04);
P = zeros(3,79); P(:,31) = [2;4;6]; P(:,32) = [1;0;3];
%% edge 1: layer 1 carries 0.4 g a->b (cell 1 -> 2), layer 2 0.1 g b->a
%% (cell 2 -> 1); edge 2: 0.2 g of species 1 from cell 2 -> 3, nothing else
dM = zeros(2,2,2);
dM(1,1,1) = 0.4; dM(1,2,1) = -0.1; dM(2,1,1) = 0.2;
Q = zeros(2,2); Q(1,1) = 1; Q(1,2) = -0.5; Q(2,1) = 1;
S = struct('topo',topo3,'P_cells',P,'dtheta_water',[0 0;0.01 -0.01;0.02 0], ...
    'dM_edge',dM,'Q_solute',Q,'theta_before',0.3*ones(3,2),'A_cell',[1;1;1], ...
    'Zlayer_m',[0.1 0.1],'active_fraction',ones(3,2),'dth',1);
tr = Lateral_Track_Record('subsurface',tr,2,2,S);
%% cell 2 (k=1), species 1: out = 0.1 + 0.2 = 0.3 g of T = 4 -> retain 0.925,
%% in = 0.4 g; cell 3 (k=2): out 0, in 0.2
check('subsurface step: retained share and inflow from the edge masses', ...
    abs(tr.retain(2,1,1) - 0.925) < 1e-12 && abs(tr.input(2,1,1) - 0.4) < 1e-12 && ...
    tr.retain(2,1,2) == 1 && abs(tr.input(2,1,2) - 0.2) < 1e-12);
check('opposite directions in two layers of one face are both counted', ...
    abs(tr.input(2,1,1) - 0.4) < 1e-12 && abs((1-tr.retain(2,1,1))*4 - 0.3) < 1e-12);
check('an empty pool keeps the water-borne export share', ...
    tr.retain(2,2,1) < 1 && tr.retain(2,2,1) >= 0, sprintf('%.4g',tr.retain(2,2,1)));
check('the water exchange is recorded per tracked cell and layer', ...
    isequal(squeeze(tr.dtheta(2,:,:)),[0.01 0.02; -0.01 0]));
%% a second step composes: T -> a2*(a1*T + b1) + b2
S2 = S; S2.P_cells(2,31) = 4*0.925 + 0.4;
S2.dM_edge = zeros(2,2,2); S2.dM_edge(2,1,1) = 0.5;  % cell 2 -> 3
S2.Q_solute = zeros(2,2); S2.Q_solute(2,1) = 1;
A1 = tr.retain(2,1,1); B1 = tr.input(2,1,1);
tr = Lateral_Track_Record('subsurface',tr,3,2,S2);
a2 = 1 - 0.5/S2.P_cells(2,31);
check('hourly maps compose exactly into the day''s map', ...
    abs(tr.retain(2,1,1) - a2*A1) < 1e-12 && abs(tr.input(2,1,1) - a2*B1) < 1e-12);
rec = Lateral_Track_Record('cell',tr,1);
check('a cell''s record is extracted with its pools', ...
    isequal(rec.pool_idx,[31 32]) && isequal(size(rec.dtheta),[4 2]));

%% --- C. replay hooks ---------------------------------------------------
NN = 49; Date = datenum(2016,6,1) + (0:NN-1)'/24;
[YE,MO,DA,HO] = datevec(Date); Datam = [YE MO DA HO];
NNd = ceil(NN/24)+1; ms = 2;
R = struct('hour_dn',Date(2:NN),'q_runon',0.5*ones(NN-1,1),'dtheta',repmat([0.02 -0.05],NN-1,1), ...
    'day_dn',floor(Date([2 26])),'retain',[0.9 1; 0.8 1],'input',[0.1 0; 0.2 0.3], ...
    'pool_idx',[31 32],'required',true);
rad = num2cell(100*ones(NN,6),1);
rep = Lateral_Replay_Prepare(R,Date,Datam,NN,NNd,ms,rad{:},1,dG,Lon,Lat);
check('the record maps onto the run''s steps and days', ...
    all(rep.q_runon(2:NN) == 0.5) && rep.q_runon(1) == 0 && ...
    isequal(rep.retain(2,:),[0.9 1]) && isequal(rep.retain(3,:),[0.8 1]));
R_short = R; R_short.hour_dn = R.hour_dn(1:end-5); R_short.q_runon = R.q_runon(1:end-5);
R_short.dtheta = R.dtheta(1:end-5,:);
check('a required run outside the record is refused', ...
    raises(@() Lateral_Replay_Prepare(R_short,Date,Datam,NN,NNd,ms,rad{:},1,dG,Lon,Lat), ...
    'Lateral_Replay:NotCovered'));
R_short.required = false;
rep_s = Lateral_Replay_Prepare(R_short,Date,Datam,NN,NNd,ms,rad{:},1,dG,Lon,Lat);
check('a seed-check window outside the record runs with no exchange there', ...
    all(rep_s.q_runon(end-4:end) == 0) && all(rep_s.dtheta(end-4:end,:) == 0,'all'));
[Oi,Vi,rep] = Lateral_Replay_Water(rep,2,[0.30 0.06],[0 0],[0.05 0.05],[0.45 0.45],[100 200]);
check('water moves within the solver''s bounds, the refused part booked', ...
    isequal(Oi,[0.32 0.05]) && abs(rep.water_mm(2) - (2 - 2)) < 1e-12 && ...
    abs(rep.water_clip_mm(2) - (-0.04*200)) < 1e-12 && isequal(Vi,(Oi - 0.05).*[100 200]));
Pj = zeros(1,79); Pj(31) = 2; Pj(32) = 1;
[Pn,rep] = Lateral_Replay_Solute(rep,3,Pj);
check('the solute map is applied and its ledger is exact', ...
    abs(Pn(31) - (0.8*2 + 0.2)) < 1e-12 && abs(Pn(32) - 1.3) < 1e-12 && ...
    max(abs((rep.solute_in_g(3,:) - rep.solute_out_g(3,:)) - (Pn([31 32]) - Pj([31 32])))) < 1e-12);

%% --- D. 2-D field run vs plot-scale replays ---------------------------------
base_file = fullfile(v2,'SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat');
if ~isfile(base_file)
    error('test_lateral_flux_tracking:NoForcing', ...
        'Forcing bundle not found: %s (see test_2d_driver_equivalence).', base_file);
end
R1 = load(base_file);
n_hours = n_days*24;
topo = Build_Cell_Topology(60,5,3,3);
z = 100 + 0.06*(max(topo.cell_yc) - topo.cell_yc) + 0.5*sin(topo.cell_xc/15);
topo = Attach_Topography_To_Topology(topo, struct('cell_z_m',z));
%% two tracked cells: the highest and the lowest of the main grid
main = find(strcmp(topo.cell_region,'main'));
[~,hi] = max(topo.cell_z_m(main)); [~,lo] = min(topo.cell_z_m(main));
tracked = [main(hi); main(lo)];
W = run_field(R1, code, topo, n_hours, tracked);
check('the tracked 2-D run completes and records every step', ...
    isstruct(W.track) && size(W.track.q_runon,1) == n_hours);
moved = max(abs(W.track.dtheta(:)));
check('the tracked cells exchange water laterally', moved > 0, sprintf('max |dtheta| %.3g',moved));
for k = 1:numel(tracked)
    c = tracked(k);
    rec = Lateral_Track_Record('cell',W.track,k);
    rec.hour_dn = rec.hour_dn(2:n_hours); rec.q_runon = rec.q_runon(2:n_hours);
    rec.dtheta = rec.dtheta(2:n_hours,:);
    days = 2:W.j;
    rec.day_dn = rec.day_dn(days); rec.retain = rec.retain(days,:); rec.input = rec.input(days,:);
    rec.slope = topo.cell_slope(c); rec.aspect_deg = topo.cell_aspect_deg(c);
    rec.svf = topo.cell_svf(c); rec.solar_window = [1 0]; rec.required = true;
    none = rec; none.q_runon(:) = 0; none.dtheta(:) = 0; none.retain(:) = 1; none.input(:) = 0;
    rep_run = run_column(R1, code, n_hours, rec);
    none_run = run_column(R1, code, n_hours, none);
    O2 = W.O(:,:,c);
    e_rep = max(abs(rep_run.O(2:end,:) - O2(2:end,:)),[],'all');
    e_none = max(abs(none_run.O(2:end,:) - O2(2:end,:)),[],'all');
    check(sprintf('cell %d: the replay follows the 2-D water contents (%.2g vs %.2g without exchange)', ...
        c,e_rep,e_none), e_rep < 2e-3 && e_rep < 0.2*e_none);
    %% each mobile pool's daily error against its own largest 2-D value
    idx = rec.pool_idx;
    P2 = W.P(2:W.j,idx,c);
    scale = max(max(abs(P2),[],1),1e-9);
    p_rep = max(max(abs(rep_run.P(2:W.j,idx) - P2),[],1)./scale);
    p_none = max(max(abs(none_run.P(2:W.j,idx) - P2),[],1)./scale);
    check(sprintf('cell %d: the replay follows the 2-D mobile pools (rel %.2g vs %.2g without)', ...
        c,p_rep,p_none), p_rep < 0.02 && p_rep < 0.2*p_none);
    ledger = sum(rep_run.lat.solute_in_g - rep_run.lat.solute_out_g,1);
    check(sprintf('cell %d: the replay''s solute ledger is finite and its water fits the bounds', c), ...
        all(isfinite(ledger)) && sum(abs(rep_run.lat.water_clip_mm)) < 1e-3*max(1,sum(abs(rep_run.lat.water_mm))));
end

fprintf('  %d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_lateral_flux_tracking: %d check(s) failed', n_fail);
end
end

%% =======================================================================
function W = run_field(R1, code, topo, n_hours, tracked)
%% The topographic 2-D driver with the ERW chemistry off and a storm, the
%% context rebuilt as test_2d_topographic_smoke does.
[NN,NNd,ms,cc,Zs,dt,dth,Datam,Date,Pr,Ta,Ds,Ws,ea,N,Pre,Tdew,SAB1,SAB2,SAD1,SAD2, ...
    PARB,PARD,Ca,Oa,Lat,Lon,DeltaGMT,t_bef,t_aft,id_location,NN_case_study_ERW_MOD, ...
    MOD_PARAM_input,ERW_PARAM_input] = forcing(R1,n_hours); %#ok<ASGLU>
OPT_ERW = 0; OPT_SoilBiogeochemistry = 1;                     %#ok<NASGU>
if isfield(R1,'B_out');      B_out = R1.B_out;           end  %#ok<NASGU>
if isfield(R1,'O_out');      O_out = R1.O_out;           end  %#ok<NASGU>
if isfield(R1,'NavlI_out');  NavlI_out = R1.NavlI_out;   end  %#ok<NASGU>
if isfield(R1,'RexmyI_out'); RexmyI_out = R1.RexmyI_out; end  %#ok<NASGU>
if isfield(R1,'NupI_L_out'); NupI_L_out = R1.NupI_L_out; end  %#ok<NASGU>
IrD = zeros(NN,1);                                            %#ok<NASGU>
Parameter_Directory = fullfile(code,'Parameters');
ColCfg = repmat(struct( ...
    'PARAM_IC', fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']), ...
    'PFT_OVERRIDE','','apply_rock',false), topo.N_cells, 1); %#ok<NASGU>
LATERAL_TRACK = struct('cells',tracked);                       %#ok<NASGU>
%% Declared here, or after the script j would be the imaginary unit.
j = 1;
old = cd(code); cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_2D;
%% field by field: struct() would spread a cell-valued field into an array
W = struct(); W.track = lateral_track; W.O = O; W.P = P; W.j = j;
end

function C = run_column(R1, code, n_hours, record)
%% The plot-scale column of the same forcing, replaying one record.
[NN,NNd,ms,cc,Zs,dt,dth,Datam,Date,Pr,Ta,Ds,Ws,ea,N,Pre,Tdew,SAB1,SAB2,SAD1,SAD2, ...
    PARB,PARD,Ca,Oa,Lat,Lon,DeltaGMT,t_bef,t_aft,id_location,NN_case_study_ERW_MOD, ...
    MOD_PARAM_input,ERW_PARAM_input] = forcing(R1,n_hours); %#ok<ASGLU>
OPT_ERW = 0; OPT_SoilBiogeochemistry = 1;                     %#ok<NASGU>
if isfield(R1,'B_out');      B_out = R1.B_out;           end  %#ok<NASGU>
if isfield(R1,'O_out');      O_out = R1.O_out;           end  %#ok<NASGU>
if isfield(R1,'NavlI_out');  NavlI_out = R1.NavlI_out;   end  %#ok<NASGU>
if isfield(R1,'RexmyI_out'); RexmyI_out = R1.RexmyI_out; end  %#ok<NASGU>
if isfield(R1,'NupI_L_out'); NupI_L_out = R1.NupI_L_out; end  %#ok<NASGU>
IrD = zeros(NN,1);                                            %#ok<NASGU>
Parameter_Directory = fullfile(code,'Parameters');
PARAM_IC = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']); %#ok<NASGU>
ERW_PARAM_IC = fullfile(Parameter_Directory,'ERW_MOD_PARAM_excel.m'); %#ok<NASGU>
LATERAL_REPLAY = record;                                       %#ok<NASGU>
old = cd(code); cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_10mins;
C = struct(); C.O = O; C.P = P; C.lat = lat_rep;
end

function [NN,NNd,ms,cc,Zs,dt,dth,Datam,Date,Pr,Ta,Ds,Ws,ea,N,Pre,Tdew,SAB1,SAB2,SAD1,SAD2, ...
    PARB,PARD,Ca,Oa,Lat,Lon,DeltaGMT,t_bef,t_aft,id_location,NN_case_study_ERW_MOD, ...
    MOD_PARAM_input,ERW_PARAM_input] = forcing(R1,n_hours)
%% The bundled forcing truncated to n_hours, with a 25 mm/h storm in hours
%% 6-9 so overland run-on and lateral exchange are substantial.
NN = n_hours; NNd = ceil(NN/24)+1;
ms = R1.ms; cc = R1.cc; Zs = R1.Zs; dt = R1.dt; dth = R1.dth;
cr = @(a) a(:, 1:min(size(a,2),n_hours));
cm = @(a) a(1:min(size(a,1),n_hours), :);
Datam = cm(R1.Datam); Date = cr(R1.Date);
Pr = cr(R1.Pr); Ta = cr(R1.Ta);
hrs = 6:min(9,numel(Pr));
Pr(hrs) = Pr(hrs) + 25; Ta(hrs) = max(Ta(hrs),4);
Ds = cr(R1.Ds); Ws = cr(R1.Ws); ea = cr(R1.ea); N = cr(R1.N); Pre = cr(R1.Pre);
Tdew = cr(R1.Tdew); SAB1 = cr(R1.SAB1); SAB2 = cr(R1.SAB2); SAD1 = cr(R1.SAD1);
SAD2 = cr(R1.SAD2); PARB = cr(R1.PARB); PARD = cr(R1.PARD); Ca = cr(R1.Ca); Oa = R1.Oa;
Lat = R1.Lat; Lon = R1.Lon; DeltaGMT = R1.DeltaGMT;
t_bef = R1.t_bef; t_aft = R1.t_aft; id_location = R1.id_location;
NN_case_study_ERW_MOD = R1.NN_case_study_ERW_MOD;
MOD_PARAM_input = R1.MOD_PARAM_input;
ERW_PARAM_input = R1.ERW_PARAM_input;
end
