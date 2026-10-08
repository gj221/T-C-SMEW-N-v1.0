%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PFT_Hamden_Deciduous_Grass.m
%%
%% Vegetation of the Hamden (145 Blue Trail, CT) site: temperate deciduous
%% broadleaf trees over grass that is cut once a year, no crops, no
%% fertiliser, no liming.  Run as a script AFTER the base parameter file
%% (MOD_PARAM_Kantola_Maize_Soy.m), named by the row's pft_override column,
%% so it replaces only the vegetation and inherits soil hydraulics, snow
%% and interception from the base file (the pattern of PFT_Prairie_Margin).
%%
%% CROWNS (cover per cell from the row: Configure_Cell_Agronomic_Inputs
%% writes Ccrown = [tree, open grass], Cbare, Curb, Cwat, Crock from the
%% land-cover map or a manual setting; this file does not touch them)
%%   crown 1  HIGH  deciduous trees     T&C site-PFT table row 181, Harvard
%%                                      Forest MA ("BDT temperate")
%%            LOW   grass beneath them  row 102, Fendt DE ("C3 grass")
%%   crown 2  HIGH  none
%%            LOW   open grass          row 102, Fendt DE
%% The table is the legacy T&C-SMEW lookup (Table_Sites_PFT_LOOKUP_latest,
%% copied with the site inputs in Inputs/Hamden/Vegetation).  Its root
%% shape, PAR threshold and leaf-shed temperature are blank for these rows
%% and stay NaN, as in the base and prairie files.
%%
%% MANAGEMENT  the grass (both crowns) is cut on 1 September every year to
%% LAI 0.5; 90 % of the cut leaf and standing dead is removed, 10 % left
%% as litter.  Trees are not cut.  No sowing, harvest or fertiliser.
%%
%% PERENNIAL  Trees and grass live from one year to the next.  The stand
%% below starts a cold start: a winter-dormant stand typical of Harvard
%% Forest -- live wood ~110 Mg C/ha (sapwood 30, heartwood 80), fine roots
%% 3, carbohydrate reserve 8 Mg C/ha, canopy 24 m, leafless -- over an
%% established sward, its roots sized so the legacy leaf-to-root ratio and
%% specific leaf area admit a mown lawn's peak LAI of about 2.5 (literature
%% priors, not Hamden measurements).  From Phase 3 on, every year of the
%% spin-up, its handoff and the transient start from the previous run's
%% end state (PFT_perennial, Perennial_Vegetation_Apply; Run_BG_Spinup
%% carries it), so roots, wood, reserves and the sward grow into the stand
%% the legacy parameters and the site sustain; leaf mass is capped at LtR
%% times the fine roots.  Phase 1, which has no nutrient uptake, restarts
%% each year from this stand.
%%
%% DEPOSITION  replaces the built-in global maps for N, K, Ca and Mg, year
%% by year (Inputs/Deposition_Hamden.mat, Build_Site_Deposition): NADP CT15
%% wet deposition arrives with each day's precipitation (the year's
%% concentration times the day's rain or snow) and TDep dry deposition
%% every day (Deposition_Day).  P, which neither measures, keeps the global
%% map.  No fertiliser, and no fertiliser years.
%%
%% WEATHERING  the base file's primary-mineral supply of Ca and Mg (the
%% Kantola uplift, Tup_Ca 0.003 and Tup_Mg 0.002 g/m2/d) keeps its charge,
%% 1.15 keq/ha/yr, but not its split, which is the Energy Farm's and gives
%% more Mg than Ca by moles.  The biogeochemistry exchanges, fixes, leaches
%% and takes up Ca and Mg alike by mass, so its soil solution settles at the
%% molar ratio of their supply; Hamden's is set so that weathering and
%% deposition supply Ca and Mg at 1.58 to 1, the dissolved ratio in which
%% the measured exchanger (Ca/Mg 3.28, Representative_soil_chemistry.csv)
%% is at equilibrium under SCEPTER's Ca/Mg selectivity (10^0.316): Tup_Ca
%% 0.0038876 and Tup_Mg 0.0014617 g/m2/d (14.2 and 5.3 kg/ha/yr).  A cold
%% start's primary-mineral pools of K, Ca and Mg begin at the steady state of
%% their supply (Tup/K3), so the release is the supply from the first day,
%% not that of the base file's placeholder inventories (502 g/m2 for each).
%%
%% EXCHANGER  Hamden's own Appelo constants replace the base file's (the
%% Energy Farm's), fitted the way those were (v2/calibration/
%% Fit_Exchange_Affinity): SCEPTER's logKH shifted by z*delta at alpha 3.4,
%% with the pinned Al(OH)3 phase, so that the exchanger emerging from the
%% class seed's spun-up soil solution (Inputs/Hamden/Calibration) holds the
%% measured base saturation, 0.598 of the CEC at pH 7, and KCl Al, 0.032, at
%% the pH they were measured at, 5.22 (Inputs/Hamden/Soils/
%% Representative_soil_chemistry.csv).  The phase lies between crystalline
%% gibbsite (log *Kso 8.04) and the microcrystalline phase (9.66); the
%% amorphous phase the base file pins (10.8) cannot hold the bases at any
%% shift.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(numel(Ccrown) == 2 && sum(Ccrown) > 0)
    error('PFT_Hamden_Deciduous_Grass:NoCover', ...
        ['The row gives no vegetated cover (Ccrown = %s): set the cell''s land ' ...
         'cover (Cell_Land_Cover) or Ccrown = [tree open-grass].'],mat2str(Ccrown));
end

%% --- Roots [mm]
ZR95_H  = [1000 0];
ZR50_H  = [NaN  NaN];
ZRmax_H = [NaN  NaN];
ZR95_L  = [250 250];
ZR50_L  = [NaN NaN];
ZRmax_L = [NaN NaN];
Zbio = MOD_PARAM_input.Zbio;
if ~(isscalar(Zbio) && isfinite(Zbio) && Zbio > 0)
    error('PFT_Hamden_Deciduous_Grass:NoActiveDepth', ...
        'The caller must supply a single positive MOD_PARAM_input.Zbio.');
end
Bio_Zs = Evaporation_layers(Zs,Zbio);
[RfH_Zs,RfL_Zs] = Root_Fraction_General(Zs,CASE_ROOT, ...
    ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%% --- Photosynthesis (both C3)
CT_H = [3 3];          CT_L = [3 3];
FI_H = [0.081 0.081];  FI_L = [0.081 0.081];
Do_H = [1000 1000];    Do_L = [1000 1000];
a1_H = [7 7];          a1_L = [6 6];
go_H = [0.01 0.01];    go_L = [0.01 0.01];
DSE_H = [0.649 0.649]; DSE_L = [0.65 0.65];
Ha_H = [76 76];        Ha_L = [55 55];
gmes_H = [Inf Inf];    gmes_L = [Inf Inf];
rjv_H = [2.4 2.4];     rjv_L = [2.2 2.2];
Vmax_H = [75 0];       %% no tree in crown 2
Vmax_L = [50 50];

%% --- Leaf traits
Sl_H = [0.026 0.026];  Sl_L = [0.030 0.030];   %% [m2/gC]
Nl_H = [28 28];        Nl_L = [23 23];         %% [gC/gN]
[Stoich_H(1)] = Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)] = Veg_Stoichiometric_Parameter(Nl_H(2));
[Stoich_L(1)] = Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)] = Veg_Stoichiometric_Parameter(Nl_L(2));
KnitH = [0.25 0.25];   KnitL = [0.15 0.15];
mSl_H = [0 0];         mSl_L = [0 0];
d_leaf_H = [4 4];      d_leaf_L = [0.8 0.8];   %% [cm]

%% --- Plant type: 1 = deciduous (seasonal) tree, 2 = perennial grass
aSE_H = [1 1];         aSE_L = [2 2];

%% --- Respiration
r_H = [0.028 0.028];   r_L = [0.06 0.06];      %% [gC/gN d]
gR_H = [0.25 0.25];    gR_L = [0.25 0.25];

%% --- Turnover [1/d]
drn_H = [1/1300 1/1300];   drn_L = [1/450 1/450];   %% fine roots
dsn_H = [0.001 0.001];     dsn_L = [1/365 1/365];   %% sapwood to heartwood
Klf_H = [1/28 1/28];       Klf_L = [1/30 1/30];     %% dead-leaf fall
Mf_H = [0.0125 0.0125];    Mf_L = [0.02 0.02];
%% Heartwood turnover, as the legacy coupled T&C-BG-2D forest runs (37.6 and
%% 40.1 years; the site table's 0 comes from water-carbon runs, where
%% heartwood is inert).  Without it heartwood grows for ever and buries the
%% nitrogen of its sapwood, 2 g N/m2 a year at Hamden, three times the
%% deposition, and the canopy declines; with it dead wood returns to the soil.
Wm_H = [1/(40*365) 1/(40*365)];   Wm_L = [0 0];

%% --- Mortality and stress
dd_max_H = [0.01 0.01];        dd_max_L = [1/45 1/45];
dc_C_H = [0.09863 0.09863];    dc_C_L = [0.019178 0.019178];
Tcold_H = [2.5 2.5];           Tcold_L = [-3 -3];

%% --- Phenology
Tlo_H = [3.5 3.5];             Tlo_L = [1.5 1.5];
Tls_H = [NaN NaN];             Tls_L = [NaN NaN];
LDay_min_H = [13.15 13.15];    LDay_min_L = [10.7 10.7];
LDay_cr_H = [11.7 11.7];       LDay_cr_L = [10.7 10.7];
mjDay_H = [250 250];           mjDay_L = [250 250];
age_cr_H = [120 120];          age_cr_L = [180 180];
dmg_H = [30 30];               dmg_L = [20 20];
Bfac_lo_H = [0.95 0.95];       Bfac_lo_L = [0.99 0.99];
Bfac_ls_H = [NaN NaN];         Bfac_ls_L = [-7 -7];
PAR_th_H = [NaN NaN];          PAR_th_L = [NaN NaN];
LAI_min_H = [0.01 0.01];       LAI_min_L = [0.1 0.1];
Trr_H = [3 3];                 Trr_L = [2.5 2.5];
eps_ac_H = [1 1];              eps_ac_L = [0.2 0.2];

%% --- Allocation
LtR_H = [0.3 0.3];             LtR_L = [0.6 0.6];
fab_H = [0.74 0.74];  fbe_H = [0.26 0.26];
fab_L = [0 0];        fbe_L = [1 1];
ff_r_H = [0.1 0.1];            ff_r_L = [0.1 0.1];

%% --- Crop-only parameters: not applicable
soCrop_H = [NaN NaN];      soCrop_L = [NaN NaN];
Sl_emecrop_H = [NaN NaN];  Sl_emecrop_L = [NaN NaN];
MHcrop_H = [NaN NaN];      MHcrop_L = [NaN NaN];

%% --- Water relations [MPa]
Psi_sto_00_H = [-0.8 -0.8];  Psi_sto_50_H = [-2.5 -2.5];
Psi_sto_00_L = [-0.5 -0.5];  Psi_sto_50_L = [-3 -3];
PsiL00_H = [-1 -1];          PsiL50_H = [-3 -3];
PsiL00_L = [-0.9 -0.9];      PsiL50_L = [-4 -4];
Kleaf_max_H = [10 10];       Kleaf_max_L = [5 5];
Cl_H = [1200 1200];          Cl_L = [1200 1200];
Axyl_H = [15 0];             Axyl_L = [0 0];
Kx_max_H = [80000 80000];    Kx_max_L = [80000 80000];
PsiX50_H = [-4.5 -4.5];      PsiX50_L = [-4.5 -4.5];
Cx_H = [150 150];            Cx_L = [150 150];
PsiG50_H = [-0.8 -0.8];      PsiG99_H = [-2.5 -2.5];
PsiG50_L = [-3 -3];          PsiG99_L = [-4 -4];
gcoef_H = [4.5 4.5];         gcoef_L = [3.5 3.5];

%% --- Optical properties: 7 broadleaf deciduous, 13 C3 grass
[PFT_opt_H(1)] = Veg_Optical_Parameter(7);
[PFT_opt_H(2)] = Veg_Optical_Parameter(7);
[PFT_opt_L(1)] = Veg_Optical_Parameter(13);
[PFT_opt_L(2)] = Veg_Optical_Parameter(13);
OM_H = [1 1];  OM_L = [1 1];

%% --- Interception
Sp_LAI_H_In = [0.2 0.2];  Sp_LAI_L_In = [0.2 0.2];

%% --- Exudation; no biological N fixation
[ParEx_H(1)] = Exudation_Parameter(0);
[ParEx_H(2)] = Exudation_Parameter(0);
[ParEx_L(1)] = Exudation_Parameter(0);
[ParEx_L(2)] = Exudation_Parameter(0);

%% --- Management: no sowing or harvest; the crop machinery stays idle.
%% Fields are overwritten on the structs the base file built, as in
%% PFT_Prairie_Margin (a fresh Vegetation_Management_Parameter no longer
%% matches the base file's structs field for field).
for islot = 1:2
    Mpar_H(islot).Date_sowing = [];      Mpar_L(islot).Date_sowing = [];
    Mpar_H(islot).Date_harvesting = [];  Mpar_L(islot).Date_harvesting = [];
    Mpar_H(islot).Crop_B = [0 0];        Mpar_L(islot).Crop_B = [0 0];
    Mpar_H(islot).fract_left = 1;        Mpar_L(islot).fract_left = 1;
    Mpar_H(islot).fract_left_fr = 1;     Mpar_L(islot).fract_left_fr = 1;
    Mpar_H(islot).fract_left_AB = 1;     Mpar_L(islot).fract_left_AB = 1;
    Mpar_H(islot).fract_left_BG = 1;     Mpar_L(islot).fract_left_BG = 1;
    Mpar_H(islot).jDay_cut = [];
end
Mpar_H(1).Crop_crown = ones(1,12);    %% trees in crown 1, all year
Mpar_H(2).Crop_crown = zeros(1,12);   %% no trees in crown 2
Mpar_L(1).Crop_crown = ones(1,12);    %% grass beneath the trees
Mpar_L(2).Crop_crown = ones(1,12);    %% open grass

%% Annual cut of the grass on 1 September: LAI capped at 0.5; 90 % of the
%% removed leaf and standing dead taken off site (fract_left 0.1).  Roots
%% and reserves are untouched, so the sward regrows.  jDay_cut holds
%% datenums (VEGGIE_UNIT compares against datenum).
for islot = 1:2
    Mpar_L(islot).jDay_cut = datenum(1975:1:2100,9,1);
    Mpar_L(islot).LAI_cut = 0.5;
    Mpar_L(islot).fract_left = 0.1;
    Mpar_L(islot).fract_left_fr = 0.1;
end
%% Residue nutrient content [N P K Ca Mg Si] (deciduous leaf litter; C3 grass)
Mpar_H(1).NPK_res_ini = [9 0.6 3 12 1.8 1];
Mpar_H(2).NPK_res_ini = [0 0 0 0 0 0];
Mpar_L(1).NPK_res_ini = [18 2.5 18 5 2 10];
Mpar_L(2).NPK_res_ini = [18 2.5 18 5 2 10];

%% --- No fertiliser.  The borrowed maize-soy row lists its maize years in
%% FertN_years, and the plot-scale frame (MAIN_FRAME_BG_ERW_10mins, every
%% spin-up stage) refills FertN with 202 kg N/ha in each of them; a stand
%% that is never fertilised lists none.
B_IO.FertN = zeros(1,366);  B_IO.FertP = zeros(1,366);  B_IO.FertK = zeros(1,366);
B_IO.FertCa = zeros(1,366); B_IO.FertMg = zeros(1,366); B_IO.FertSi = zeros(1,366);
B_IO.FertN_years = [];
%% --- Deposition in place of the global maps for N, K, Ca and Mg: NADP CT15
%% wet deposition with each day's precipitation and TDep dry deposition
%% every day, both year by year (Build_Site_Deposition, Deposition_Day);
%% the constant rates are the 2015-2025 means.
deposition_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'Inputs','Deposition_Hamden.mat');
if ~isfile(deposition_file)
    error('PFT_Hamden_Deciduous_Grass:MissingDeposition', ...
        'Hamden deposition not found: %s (Build_Site_Deposition).',deposition_file);
end
hamden_dep = load(deposition_file,'Deposition');
B_IO.DepN = hamden_dep.Deposition.mean_g_m2_d(1);
B_IO.DepK = hamden_dep.Deposition.mean_g_m2_d(2);
B_IO.DepCa = hamden_dep.Deposition.mean_g_m2_d(3);
B_IO.DepMg = hamden_dep.Deposition.mean_g_m2_d(4);
B_IO.Dep_events = hamden_dep.Deposition.events;

%% --- Weathering (header, WEATHERING)
B_IO.Tup_Ca = 0.0038876;   % [g Ca/m2 d]
B_IO.Tup_Mg = 0.0014617;   % [g Mg/m2 d]
if exist('OPT_SoilBiogeochemistry','var') && OPT_SoilBiogeochemistry == 1 && ...
        ~exist('B_out','var')
    hamden_bgp = Biogeochemistry_Parameter(1);
    P(1,[55 63 71]) = [B_IO.Tup_K/hamden_bgp.K3k B_IO.Tup_Ca/hamden_bgp.K3ca ...
        B_IO.Tup_Mg/hamden_bgp.K3mg];
end

%% --- Plant reserves: each element at the storage capacity of the stand's
%% carbohydrate reserve, reserve x (1/leaf - 1/storage) -- the state
%% Crop_Planting_Harvest gives a sown crop, so every relative-concentration
%% factor starts neutral (PFT_Prairie_Margin explains why the base file's
%% legacy 1000 g/m2 reserves must not stand).
tree_B = [0 3000 300 800 0 8000 0 0];      %% [gC/m2]: leaves sapwood fineroot reserve fruit heartwood deadleaf idle
grass_B = [0 0 140 93 0 0 0 0];            %% an established sward: roots 2.5/(Sl_L*LtR_L),
                                           %% reserve at its cap of 2/3 the roots
reserve_at = @(St,B4) max(0,B4 .* ...
    ([1/St.Nl 1/St.Phol 1/St.Kpotl 1/St.Call 1/St.Magl 1/St.Sill] - ...
     [1/St.Ns 1/St.Phos 1/St.Kpots 1/St.Cals 1/St.Mags 1/St.Sils]));
res_H = reserve_at(Stoich_H(1),tree_B(4));
res_L = reserve_at(Stoich_L(1),grass_B(4));
Nreserve_H(1,:)  = [res_H(1) 0];   Nreserve_L(1,:)  = [res_L(1) res_L(1)];
Preserve_H(1,:)  = [res_H(2) 0];   Preserve_L(1,:)  = [res_L(2) res_L(2)];
Kreserve_H(1,:)  = [res_H(3) 0];   Kreserve_L(1,:)  = [res_L(3) res_L(3)];
Careserve_H(1,:) = [res_H(4) 0];   Careserve_L(1,:) = [res_L(4) res_L(4)];
Mgreserve_H(1,:) = [res_H(5) 0];   Mgreserve_L(1,:) = [res_L(5) res_L(5)];
Sireserve_H(1,:) = [res_H(6) 0];   Sireserve_L(1,:) = [res_L(6) res_L(6)];

%% --- Initial stand: dormant (1 January)
B_H(1,1,:) = tree_B;       B_H(1,2,:) = zeros(1,8);
B_L(1,1,:) = grass_B;      B_L(1,2,:) = grass_B;
LAI_H(1,:) = [0 0];        LAI_L(1,:) = [0 0];
hc_H(1,:) = [24 0];        hc_L(1,:) = [0.05 0.05];   %% [m]
SAI_H(1,:) = [0.5 0];      SAI_L(1,:) = [0.001 0.001];
PHE_S_H(1,:) = [1 0];      PHE_S_L(1,:) = [1 1];
AgeL_H(1,:) = [0 0];       AgeL_L(1,:) = [0 0];
e_rel_H(1,:) = [1 0];      e_rel_L(1,:) = [1 1];

%% --- Perennial: a run that follows another starts from its end state.
PFT_perennial = true;
Perennial_Vegetation_Apply;

%% --- Growth without T&C's sink limitation, as the legacy T&C-BG-2D forest
%% runs (OPT_EnvLimitGrowth = 0, T&C's own default).  The v0.9 plot-scale
%% driver switches it on for every site before running this file; with it,
%% NPP above the environmental growth limit goes to the carbohydrate
%% reserve past its cap of 2/3 of the sapwood (VEGETATION_DYNAMIC, Add_AR),
%% so a stand that lives from year to year piles carbon into its reserve
%% while its canopy stalls.
OPT_EnvLimitGrowth = 0;

%% --- Exchanger (header, EXCHANGER)
if ~exist('bulk_pH_params','var') || ~isstruct(bulk_pH_params)
    bulk_pH_params = struct();
end
hamden_delta = -1.042050;
bulk_pH_params.appelo = struct('logKH',[5.9 4.8 10.47 10.786 16.47] + ...
    [1 1 2 2 3]*hamden_delta,'alpha',3.4);
bulk_pH_params.Al_pin_logKso = 9.203355;

%% --- Reference height of the forcing [m]: 2 m above the canopy, T&C's
%% convention over tall vegetation (the driver's own guard adds 2 m when a
%% canopy reaches it).  The base file's 8 m suits maize; under a 24 m canopy
%% (zero-plane displacement ~16 m) the aerodynamic resistance would have no
%% height to work over.  The forcing's near-surface air and wind are
%% applied at it.
zatm = max([hc_H(1,:) hc_L(1,:)]) + 2;
fprintf(['[PFT_Hamden_Deciduous_Grass] trees (Harvard Forest, ZR95 %d mm) over grass ' ...
    '(Fendt, ZR95 %d mm); cover %s; grass cut 1 Sep (90%% removed); CT15 wet + TDep dry deposition by event; ' ...
    'reference height %g m; exchange shift %.4f per charge, Al(OH)3 log *Kso %.4f\n'], ...
    ZR95_H(1),ZR95_L(1),mat2str(Ccrown,4),zatm,hamden_delta,bulk_pH_params.Al_pin_logKso);
clear islot tree_B grass_B reserve_at res_H res_L hamden_dep deposition_file hamden_delta hamden_bgp
