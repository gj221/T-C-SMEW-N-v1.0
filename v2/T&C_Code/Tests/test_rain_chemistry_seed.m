function test_rain_chemistry_seed()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's rain-chemistry record (NADP NTN) and the ERW seed it feeds.
%%   1  Load_Rain_Chemistry: no record without a rain_chem_file; a named
%%      record is read from Inputs/Rain_Chem_<id_location>.mat; a missing
%%      or invalid record is refused
%%   2  Hamden's row names CT15's record: rain pH 5.20, SO4 + Cl
%%      44.1 umolc m-2 d-1, Na 23.7 umol m-2 d-1
%%   3  ERW_Mean_Water_Loss is the drivers' Tw + Lw, per day
%%   4  the drivers hand the record and the spin-up's long-term water loss
%%      on: the spin-up records the ring's water loss and hands its mean
%%      off with the state, the 1-D driver and the 2-D frame and step pass
%%      both to the seed (the 2-D frame the column's state, else the
%%      site's), and a record without the water loss is refused
%% (test_kantola_seed_charge runs the seed itself on a record.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_directory = fileparts(mfilename('fullpath'));
code = fileparts(test_directory);
v2 = fileparts(code);
addpath(code);

%% --- 1
T = table("X",'VariableNames',{'Location'});
S = Load_Rain_Chemistry(T,tempdir,'X');
assert(isempty(S.source_file) && isnan(S.rain_pH_vwm) && ...
    S.rain_Na_flux_umol_m2_d == 0 && S.rain_strong_anion_flux_umolc_m2_d == 0, ...
    'Without rain_chem_file there is no record.');
for blank = {"", "none", NaN}
    Tb = T; Tb.rain_chem_file = blank{1};
    assert(isempty(Load_Rain_Chemistry(Tb,tempdir,'X').source_file), ...
        'An empty or ''none'' rain_chem_file names no record.');
end
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d,'s'));
Tr = T; Tr.rain_chem_file = "Rain_Chem_X.mat";
assert_error(@() Load_Rain_Chemistry(Tr,d,'X'),'TC_SMEW_N_run:MissingRainChemistry');
Rain_Chem = struct('pH_vwm',5.1,'strong_anion_ueq_L',20, ...
    'strong_anion_deposition_umolc_m2',1e4,'strong_anion_flux_umolc_m2_d',50, ...
    'Ca_flux_g_m2_d',1e-4,'Mg_flux_g_m2_d',5e-5,'K_flux_g_m2_d',4e-5, ...
    'Na_flux_g_m2_d',5e-4,'Na_flux_umol_m2_d',21.7);
save(fullfile(d,'Rain_Chem_X.mat'),'Rain_Chem');
S = Load_Rain_Chemistry(Tr,d,'X');
assert(strcmp(S.source_file,'Rain_Chem_X.mat') && S.rain_pH_vwm == 5.1 && ...
    S.rain_Na_flux_umol_m2_d == 21.7 && S.rain_strong_anion_flux_umolc_m2_d == 50, ...
    'A named record is read from Rain_Chem_<id_location>.mat.');
Rain_Chem.Na_flux_umol_m2_d = NaN;
save(fullfile(d,'Rain_Chem_X.mat'),'Rain_Chem');
assert_error(@() Load_Rain_Chemistry(Tr,d,'X'),'TC_SMEW_N_run:InvalidRainChemistry');

%% --- 2
grid_dir = fullfile(code,'Inputs','Field_Spinup','Hamden');
M = Read_Parameter_Columns(fullfile(grid_dir,'Hamden_MOD_columns.csv'));
H = Load_Rain_Chemistry(struct2table(M),fullfile(code,'Inputs'),M.Location);
assert(strcmp(H.source_file,'Rain_Chem_Hamden.mat') && ...
    abs(H.rain_pH_vwm-5.203) < 1e-3 && ...
    abs(H.rain_strong_anion_flux_umolc_m2_d-44.14) < 0.01 && ...
    abs(H.rain_Na_flux_umol_m2_d-23.74) < 0.01, ...
    'Hamden''s row must name CT15''s record.');

%% --- 3
Bio_Zs = [0.4 0.6 0 0];                         % the layer: the top two
RfH_Zs = [0.5 0.3 0.2 0; 0.1 0.2 0.3 0.4];      % two crowns
RfL_Zs = [0.7 0.3 0 0; 0.25 0.25 0.25 0.25];
T_H = [0 0; 0.10 0.20; 0.05 0.00];              % [mm/h]; row 1 initial
T_L = [0 0; 0.02 0.04; 0.00 0.06];
Lk = [9; 0.01; 0.03];
cfH = (sum((ones(2,1)*(Bio_Zs>0)).*RfH_Zs,2))'; % as MAIN_FRAME_BG_ERW_10mins
cfL = (sum((ones(2,1)*(Bio_Zs>0)).*RfL_Zs,2))';
Tw = sum(cfH.*T_H+cfL.*T_L,2)/1000;             % [m/h]
Lw = Lk/1000;
expected = mean(Tw(2:end)+Lw(2:end))*24*1000;   % [L m-2 d-1]
assert(abs(ERW_Mean_Water_Loss(T_H,T_L,Lk,Bio_Zs,RfH_Zs,RfL_Zs)-expected) < 1e-12, ...
    'ERW_Mean_Water_Loss must be the drivers'' Tw + Lw per day.');

%% --- 4
spinup = fileread(fullfile(code,'Run_BG_Spinup.m'));
site = fileread(fullfile(v2,'spinup','Run_Automatic_Site_Spinup.m'));
frame_1d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_10mins.m'));
frame_2d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_2D.m'));
step_2d = fileread(fullfile(code,'ERW_Daily_Step_2D.m'));
subfield = fileread(fullfile(v2,'run','Run_One_Subfield.m'));
assert(contains(spinup,'p3_water_loss_by_ring(ring_pos) = ERW_Mean_Water_Loss(') && ...
    contains(spinup,'bg_state.p3_water_loss_by_ring = p3_water_loss_by_ring;') && ...
    contains(spinup,'bg_state_no_erw.seed_water_loss_l_m2_d = ERW_seed_water_loss_l_m2_d;') && ...
    contains(spinup,'Run_BG_Spinup:NoSeedWaterLoss'), ...
    'The spin-up records the ring''s water loss and hands its mean off.');
assert(~contains(site,'function S = Load_Rain_Chemistry') && ...
    contains(site,'ERW_seed_water_loss_l_m2_d = bg_state.seed_water_loss_l_m2_d;') && ...
    contains(site,'TC_SMEW_N_run:NoSeedWaterLoss') && ...
    contains(frame_1d,'ERW_TC_const.seed_water_loss_l_m2_d = ERW_seed_water_loss_l_m2_d;'), ...
    'The 1-D driver reads the shared record and seeds with the state''s water loss.');
assert(contains(frame_2d,'PAR{c}.rain_Na_flux_umol_m2_d = ColCfg(c).RAIN_CHEM.rain_Na_flux_umol_m2_d;') && ...
    contains(frame_2d,'PAR{c}.seed_water_loss_l_m2_d = ColCfg(c).BG_STATE.seed_water_loss_l_m2_d;') && ...
    contains(frame_2d,'PAR{c}.seed_water_loss_l_m2_d = seed_water_loss_l_m2_d;') && ...
    contains(frame_2d,'MAIN_FRAME_BG_ERW_2D:NoSeedWaterLoss') && ...
    contains(frame_2d,'sprintf('';rain=%s''') && ...
    contains(step_2d,'C.rain_Na_flux_umol_m2_d = p.rain_Na_flux_umol_m2_d;') && ...
    contains(step_2d,'C.rain_strong_anion_flux_umolc_m2_d = p.rain_strong_anion_flux_umolc_m2_d;') && ...
    contains(step_2d,'C.seed_water_loss_l_m2_d = p.seed_water_loss_l_m2_d;') && ...
    contains(subfield,'''RAIN_CHEM'',rain_chemistry'), ...
    'The 2-D frame and step take a column''s record and its state''s water loss.');

fprintf('test_rain_chemistry_seed: all checks passed\n');
end

function assert_error(f,id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,id),'Expected %s, got %s.',id,ME.identifier);
    return
end
error('test_rain_chemistry_seed:NoError','Expected %s.',id);
end
