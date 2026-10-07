%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  MULTI-COLUMN MAIN FRAME  --  T&C-SMEW-N sub-plot domains %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Drives one soil-vegetation column per sub-plot cell of a
%% Build_Cell_Topology partition, coupling the columns laterally every
%% hydrologic timestep.
%%
%% RELATION TO THE SINGLE-COLUMN MAIN FRAME
%%   The science is unchanged: this script calls exactly the same unit
%%   functions as MAIN_FRAME_BG_ERW_10mins (HYDROLOGIC_UNIT, BIOGEO_UNIT,
%%   VEGGIE_UNIT, Plant_Exports, CropHeightType, ERW_function) with the
%%   same arguments.  What changes is that every state array carries a
%%   trailing column dimension and the columns exchange water and
%%   dissolved solute and heat after each hourly step.  With N_col = 1 and a
%%   single-cell topology it reproduces the single-column result, which
%%   is how it is verified (test_2d_driver_equivalence).
%%
%% WHAT IS RESOLVED IN 2D, AND WHAT IS NOT
%%   Water     lateral Buckingham-Darcy exchange between adjacent cells,
%%             per soil layer, every hour (Darcy_Lateral_Topology).
%%   Mass      advective transport of the aqueous phase with that water,
%%             upwind (Solute_Lateral_Transport).  Sorbed and solid
%%             phases do not move laterally. ERW-only Na, Al, residual
%%             anions and DIC use the same water flux in the active layer.
%%   Energy    lateral Fourier conduction and liquid-water sensible heat
%%             advection couple the soil thermal profiles.
%%   Gas       soil-air CO2 diffuses across cell faces with the existing
%%             Millington-Quirk effective diffusivity.
%%
%% METEOROLOGY
%%   The plot is flat and small, so all columns share one meteorological
%%   forcing series.  Sub-domains diverge through vegetation, soil state
%%   and treatment, not through their driving weather.
%%
%% EXPECTED IN THE CALLER WORKSPACE
%%   topo                Build_Cell_Topology struct (N_col = topo.N_cells)
%%   ColCfg              [N_col x 1] struct: .PARAM_IC .PFT_OVERRIDE
%%                       .apply_rock (logical)
%%   Parameter_Directory path holding the parameter scripts
%%   MOD_PARAM_input, ERW_PARAM_input, NN, NNd, ms, Zs, cc, dt, dth,
%%   Datam, Date, and the meteorological series, as for the 1-D frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_col = topo.N_cells;
if ~exist('parallel_workers','var') || isempty(parallel_workers)
    parallel_workers = 0;
end
fprintf('\n=== 2D driver: %d columns over %.0f m2 ===\n', ...
    N_col, sum(topo.cell_area_m2));

if ~exist('OPT_SoilBiogeochemistry','var'); OPT_SoilBiogeochemistry = 1; end
if ~exist('OPT_ERW','var');                 OPT_ERW = 1;                 end
if ~exist('IrD','var') || isempty(IrD);     IrD = zeros(NN,1);           end

NCP = 8;    %% number of carbon pools

%% ========================================================================
%% PER-COLUMN PARAMETERS
%%   Each column's parameter script runs in its own private workspace, so
%%   columns may differ in PFT, management, rooting depth and soil.
%% ========================================================================
%% Build field by field: struct('f',M) with a non-scalar M would create a
%% struct ARRAY the size of M, not a scalar struct holding M.
ctx = struct();
ctx.MOD_PARAM_input = MOD_PARAM_input;
ctx.ERW_PARAM_input = ERW_PARAM_input;
ctx.NN = NN;   ctx.NNd = NNd;   ctx.ms = ms;
ctx.Datam = Datam;
ctx.Lat = Lat; ctx.Lon = Lon;   ctx.DeltaGMT = DeltaGMT;
ctx.t_bef = t_bef; ctx.t_aft = t_aft;
ctx.Ta = Ta;   ctx.Ca = Ca;
ctx.OPT_SoilBiogeochemistry = OPT_SoilBiogeochemistry;
ctx.NN_case_study_ERW_MOD = NN_case_study_ERW_MOD;
if exist('B_out','var');      ctx.B_out = B_out;           end
if exist('O_out','var');      ctx.O_out = O_out;           end
if exist('NavlI_out','var');  ctx.NavlI_out = NavlI_out;   end
if exist('RexmyI_out','var'); ctx.RexmyI_out = RexmyI_out; end
if exist('NupI_L_out','var'); ctx.NupI_L_out = NupI_L_out; end
%%% v1.0: a perennial stand's end state (Perennial_Vegetation_Apply).
if exist('Veg_out','var');    ctx.Veg_out = Veg_out;       end

PAR = cell(N_col,1);
biogeo_exchange_scheme = repmat({'TC_GT'},N_col,1);
biogeo_KNmaxd = cell(N_col,1);
biogeo_overrides = cell(N_col,1);
for c = 1:N_col
    cell_ctx = ctx;
    if isfield(ColCfg(c),'MOD_PARAM') && istable(ColCfg(c).MOD_PARAM) && ...
            height(ColCfg(c).MOD_PARAM) == 1
        cell_ctx.MOD_PARAM_input = ColCfg(c).MOD_PARAM;
    end
    if isfield(ColCfg(c),'ERW_PARAM') && istable(ColCfg(c).ERW_PARAM) && ...
            height(ColCfg(c).ERW_PARAM) == 1
        cell_ctx.ERW_PARAM_input = ColCfg(c).ERW_PARAM;
    end
    if isfield(ColCfg(c),'BG_STATE') && isstruct(ColCfg(c).BG_STATE) && ...
            ~isempty(fieldnames(ColCfg(c).BG_STATE))
        cell_ctx = Cell_BG_State_Context(cell_ctx,ColCfg(c).BG_STATE);
    end
    ovr = '';
    if isfield(ColCfg(c),'PFT_OVERRIDE') && ~isempty(ColCfg(c).PFT_OVERRIDE)
        ovr = fullfile(Parameter_Directory, [char(ColCfg(c).PFT_OVERRIDE) '.m']);
    end
    PAR{c} = Load_Cell_Parameters(char(ColCfg(c).PARAM_IC), ovr, cell_ctx);
    if isfield(PAR{c},'bulk_pH_params') && ...
            isfield(PAR{c}.bulk_pH_params,'scheme')
        biogeo_exchange_scheme{c} = PAR{c}.bulk_pH_params.scheme;
    end
    if isfield(ColCfg(c),'KNmaxd')
        biogeo_KNmaxd{c} = ColCfg(c).KNmaxd;
    end
    if isfield(ColCfg(c),'biogeo_overrides')
        biogeo_overrides{c} = ColCfg(c).biogeo_overrides;
    end
    if isfield(ColCfg(c),'appelo_alpha') && ~isempty(ColCfg(c).appelo_alpha)
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'AppeloAlpha requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.appelo.alpha = ColCfg(c).appelo_alpha;
    end
    if isfield(ColCfg(c),'cec_dpH') && ~isempty(ColCfg(c).cec_dpH)
        %%% Variable charge is a property of the soil, not of the exchange
        %%% scheme, but the sensitivity that reproduces a given titration
        %%% curve depends on the selectivity constants it is fitted with.
        %%% A scalar supplied here applies to whichever scheme is active;
        %%% per-scheme values are declared in the parameter file instead.
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'CECpHSensitivity requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.dCEC_dpH = ColCfg(c).cec_dpH;
    end
    if isfield(ColCfg(c),'appelo_logKH') && ~isempty(ColCfg(c).appelo_logKH)
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'AppeloLogKH requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.appelo.logKH = ColCfg(c).appelo_logKH;
    end
    %%% v1.0: a site's rain-chemistry record (Load_Rain_Chemistry) gives the
    %%% column's ERW its rain pH and its Na and strong-anion supply, as in the
    %%% 1-D driver, and the column's spin-up state -- its own, or else the
    %%% site's -- the long-term water loss the rain Na seed divides by.
    %%% Inert without a record.
    if isfield(ColCfg(c),'RAIN_CHEM') && isstruct(ColCfg(c).RAIN_CHEM) && ...
            isfield(ColCfg(c).RAIN_CHEM,'source_file') && ...
            ~isempty(strtrim(char(ColCfg(c).RAIN_CHEM.source_file)))
        if isfield(ColCfg(c),'BG_STATE') && isstruct(ColCfg(c).BG_STATE) && ...
                isfield(ColCfg(c).BG_STATE,'seed_water_loss_l_m2_d')
            PAR{c}.seed_water_loss_l_m2_d = ColCfg(c).BG_STATE.seed_water_loss_l_m2_d;
        elseif exist('seed_water_loss_l_m2_d','var')
            PAR{c}.seed_water_loss_l_m2_d = seed_water_loss_l_m2_d;
        else
            error('MAIN_FRAME_BG_ERW_2D:NoSeedWaterLoss', ...
                ['Column %d: the rain record %s seeds Na from the spin-up''s ' ...
                 'long-term water loss, which its spin-up state does not carry.'], ...
                c,char(ColCfg(c).RAIN_CHEM.source_file));
        end
        PAR{c}.rain_pH_vwm = ColCfg(c).RAIN_CHEM.rain_pH_vwm;
        PAR{c}.rain_Na_flux_umol_m2_d = ColCfg(c).RAIN_CHEM.rain_Na_flux_umol_m2_d;
        PAR{c}.rain_strong_anion_flux_umolc_m2_d = ...
            ColCfg(c).RAIN_CHEM.rain_strong_anion_flux_umolc_m2_d;
    end
    fprintf('  column %2d  %-14s  %s\n', c, topo.cell_name{c}, ...
        ternary(isempty(ovr),'base parameters','PFT override'));
end
if exist('run_provenance','var') && isstruct(run_provenance)
    appelo_alpha_effective = nan(N_col,1);
    appelo_logKH_effective = nan(N_col,5);
    for c = 1:N_col
        if strcmpi(biogeo_exchange_scheme{c},'appelo')
            appelo_params = [];
            if isfield(PAR{c},'bulk_pH_params') && ...
                    isfield(PAR{c}.bulk_pH_params,'appelo')
                appelo_params = PAR{c}.bulk_pH_params.appelo;
            end
            [~,~,~,~,~,appelo_alpha_effective(c), ...
                appelo_logKH_effective(c,:)] = ...
                Resolve_Appelo_Constants(appelo_params);
        end
    end
    cec_dpH_effective = nan(N_col,1);
    cec_ref_pH_effective = nan(N_col,1);
    for c = 1:N_col
        if isfield(PAR{c},'bulk_pH_params')
            [cec_dpH_effective(c),cec_ref_pH_effective(c)] = ...
                ERW_Resolve_CEC_pH_Sensitivity(PAR{c}.bulk_pH_params, ...
                    char(biogeo_exchange_scheme{c}));
        end
    end
    run_provenance.exchange_scheme = biogeo_exchange_scheme;
    run_provenance.AppeloAlpha = appelo_alpha_effective;
    run_provenance.AppeloLogKH = appelo_logKH_effective;
    run_provenance.CECpHSensitivity = cec_dpH_effective;
    run_provenance.CECrefpH = cec_ref_pH_effective;
end

%% Soil layering must agree across columns: the lateral solver exchanges
%% water layer by layer, which is only defined on a shared discretisation.
for c = 2:N_col
    if ~isequal(PAR{c}.Zs, PAR{1}.Zs)
        error('MAIN_FRAME_BG_ERW_2D:LayerMismatch', ...
            ['Column %d uses a different soil layer mesh from column 1. ' ...
             'Lateral exchange is computed layer by layer and is undefined ' ...
             'across differing meshes.'], c);
    end
end
if any(abs(cellfun(@(p) p.Zbio,PAR)-PAR{1}.Zbio) > 1e-9)
    error('MAIN_FRAME_BG_ERW_2D:ActiveDepthMismatch', ...
        'All columns must use the same active biogeochemical depth.');
end
Zlayer_m = diff(PAR{1}.Zs)/1000;      %% [m] layer thickness, shared

%% Topography enters through the cell topology (Attach_Topography_To_Topology
%% adds cell_z_m and the derived attributes), NOT through the per-column
%% parameter files.  The FV lateral solver (Darcy_Lateral_Topology) drives
%% subsurface flow on total hydraulic head H = psi + z, and the surface
%% router moves ponded runoff down the DEM; the vertical column stays 1-D
%% (Slo_top = 0 in the HYDROLOGIC_UNIT call below), so its internal
%% kinematic Qi_out cannot double-count the graph exchange.  A non-zero
%% parameter-file Slo_top would re-introduce that double count and is
%% rejected.
topography_attached = isfield(topo,'topography_attached') && topo.topography_attached;
if N_col > 1
    for c = 1:N_col
        if abs(PAR{c}.Slo_top) > 1e-12 || any(abs(PAR{c}.Slo_pot) > 1e-12)
            error('MAIN_FRAME_BG_ERW_2D:SlopedDomainUnsupported', ...
                ['Column %d carries a non-zero parameter-file Slo_top/Slo_pot. ' ...
                 'Put field topography in the cell topology via ' ...
                 'Attach_Topography_To_Topology (topo.cell_z_m), not in the ' ...
                 'parameter file.'], c);
        end
    end
    if topography_attached
        fprintf(['  topography attached: %s; slope %.1f-%.1f%%, ' ...
            'relief %.1f m, %d outlet cell(s)\n'], ...
            char(string(topo.topography_source)), ...
            100*min(topo.cell_slope), 100*max(topo.cell_slope), ...
            max(topo.cell_z_m)-min(topo.cell_z_m), sum(topo.cell_is_outlet));
    else
        fprintf('  no topography attached: flat-field lateral coupling\n');
    end
end

bio_active_fraction = zeros(N_col,ms);
for c = 1:N_col
    active_thickness_m = PAR{c}.Bio_Zs*(PAR{c}.Zbio/1000);
    bio_active_fraction(c,:) = active_thickness_m./Zlayer_m;
end
if any(~isfinite(bio_active_fraction(:))) || ...
        any(bio_active_fraction(:) < -1e-12) || ...
        any(bio_active_fraction(:) > 1+1e-12)
    error('MAIN_FRAME_BG_ERW_2D:BadActiveLayer', ...
        'Biogeochemical active-layer fractions must lie within each soil layer.');
end
bio_active_fraction = min(1,max(0,bio_active_fraction));

%% Per-column soil parameters for the lateral solver
soilpar_cell = Soil_Param_To_Routing(PAR{1}.Soil_Param, PAR{1}.SPAR, PAR{1}.aR);
soilpar_cell = repmat(soilpar_cell, N_col, 1);
for c = 1:N_col
    soilpar_cell(c) = Soil_Param_To_Routing(PAR{c}.Soil_Param, PAR{c}.SPAR, PAR{c}.aR);
end
A_cell = topo.cell_area_m2(:);             %% [m2]
A_tot  = sum(A_cell);

%% ========================================================================
%% STATE ALLOCATION   (trailing dimension = column)
%% ========================================================================
z2 = @(a,b) zeros(a,b,N_col);
z3 = @(a,b,d) zeros(a,b,d,N_col);

V=z2(NN,ms); O=z2(NN,ms); Qi_out=z2(NN,ms); WTR=z2(NN,ms); POT=z2(NN,ms);
Tdp=z2(NN,ms); Vice=z2(NN,ms); Oice=z2(NN,ms);
OF=z2(NN,1); OS=z2(NN,1); ZWT=z2(NN,1);
Pr_sno=z2(NN,1); Pr_liq=z2(NN,1); er=z2(NN,1);
Ts=z2(NN,1); ra=z2(NN,1); r_soil=z2(NN,1); b_soil=z2(NN,1); Tice=z2(NN,1);
Rn=z2(NN,1); H=z2(NN,1); QE=z2(NN,1); Qv=z2(NN,1); Lpho=z2(NN,1);
EG=z2(NN,1); ESN=z2(NN,1); ESN_In=z2(NN,1); EWAT=z2(NN,1);
EIn_urb=z2(NN,1); EIn_rock=z2(NN,1); dw_SNO=z2(NN,1);
G=z2(NN,1); SWE=z2(NN,1); SND=z2(NN,1); ros=z2(NN,1);
In_SWE=z2(NN,1); SP_wc=z2(NN,1); WR_SP=z2(NN,1); U_SWE=z2(NN,1); NIn_SWE=z2(NN,1);
dQ=z2(NN,1); Qfm=z2(NN,1); t_sls=z2(NN,1); DQ=z2(NN,1); DT=z2(NN,1);
In_urb=z2(NN,1); In_rock=z2(NN,1); SE_rock=z2(NN,1); SE_urb=z2(NN,1);
Csno=z2(NN,1); WIS=z2(NN,1); Lk=z2(NN,1); f=z2(NN,1); Rh=z2(NN,1); Rd=z2(NN,1);
NDVI=z2(NN,1); alp_soil=z2(NN,1); tau_sno=z2(NN,1); e_sno=z2(NN,1); ALB=z2(NN,1);
EK=z2(NN,1); dQVEG=z2(NN,1); TsV=z2(NN,1); HV=z2(NN,1); QEV=z2(NN,1);
Cice=z2(NN,1); Lk_wat=z2(NN,1); Lk_rock=z2(NN,1); EICE=z2(NN,1);
WAT=z2(NN,1); ICE=z2(NN,1); ICE_D=z2(NN,1); WR_IP=z2(NN,1); NIce=z2(NN,1);
Cicew=z2(NN,1); IP_wc=z2(NN,1); Csnow=z2(NN,1); FROCK=z2(NN,1);
Imelt=z2(NN,1); Smelt=z2(NN,1); Tdamp=z2(NN,1); Gfin=z2(NN,1);
In_Litter=z2(NN,1); ELitter=z2(NN,1); Ws_under=z2(NN,1);
r_litter=z2(NN,cc);
In_H=z2(NN,cc); In_L=z2(NN,cc); OH=z2(NN,cc); OL=z2(NN,cc);
Psi_s_H=z2(NN,cc); Psi_s_L=z2(NN,cc); Tdp_H=z2(NN,cc); Tdp_L=z2(NN,cc);
rap_H=z2(NN,cc); rap_L=z2(NN,cc); Dr_H=z2(NN,cc); Dr_L=z2(NN,cc);
T_H=z2(NN,cc); T_L=z2(NN,cc); EIn_H=z2(NN,cc); EIn_L=z2(NN,cc);
rb_H=z2(NN,cc); rb_L=z2(NN,cc);
An_H=z2(NN,cc); Rdark_H=z2(NN,cc); An_L=z2(NN,cc); Rdark_L=z2(NN,cc);
Ci_sunH=z2(NN,cc); Ci_sunL=z2(NN,cc); Ci_shdH=z2(NN,cc); Ci_shdL=z2(NN,cc);
rs_sunH=z2(NN,cc); rs_sunL=z2(NN,cc); rs_shdH=z2(NN,cc); rs_shdL=z2(NN,cc);
Vx_H=z2(NN,cc); Vl_H=z2(NN,cc); Vx_L=z2(NN,cc); Vl_L=z2(NN,cc);
Psi_x_H=z2(NN,cc); Psi_l_H=z2(NN,cc); Psi_x_L=z2(NN,cc); Psi_l_L=z2(NN,cc);
gsr_H=z2(NN,cc); Jsx_H=z2(NN,cc); Jxl_H=z2(NN,cc); Kleaf_H=z2(NN,cc); Kx_H=z2(NN,cc);
gsr_L=z2(NN,cc); Jsx_L=z2(NN,cc); Jxl_L=z2(NN,cc); Kleaf_L=z2(NN,cc); Kx_L=z2(NN,cc);
fapar_H=z2(NN,cc); fapar_L=z2(NN,cc); SIF_H=z2(NN,cc); SIF_L=z2(NN,cc);
Tdeb = zeros(NN, max(1,length(PAR{1}.Zs_deb)-1), N_col);

%% Daily vegetation and biogeochemistry state
LAI_L=z2(NNd,cc); B_L=z3(NNd,cc,NCP); NPP_L=z2(NNd,cc); Rg_L=z2(NNd,cc);
RA_L=z2(NNd,cc); Rmc_L=z2(NNd,cc); ANPP_L=z2(NNd,cc);
Rms_L=z2(NNd,cc); Rmr_L=z2(NNd,cc); PHE_S_L=z2(NNd,cc); dflo_L=z2(NNd,cc);
AgeL_L=z2(NNd,cc); e_rel_L=ones(NNd,cc,N_col); SAI_L=z2(NNd,cc); hc_L=z2(NNd,cc);
e_relN_L=ones(NNd,cc,N_col); BA_L=z2(NNd,cc);
LAI_H=z2(NNd,cc); B_H=z3(NNd,cc,NCP); NPP_H=z2(NNd,cc); Rg_H=z2(NNd,cc);
RA_H=z2(NNd,cc); Rms_H=z2(NNd,cc); Rmr_H=z2(NNd,cc); ANPP_H=z2(NNd,cc);
PHE_S_H=z2(NNd,cc); dflo_H=z2(NNd,cc); Rmc_H=z2(NNd,cc);
AgeL_H=z2(NNd,cc); e_rel_H=ones(NNd,cc,N_col); SAI_H=z2(NNd,cc); hc_H=z2(NNd,cc);
e_relN_H=ones(NNd,cc,N_col); BA_H=z2(NNd,cc);
Rrootl_H=z2(NNd,cc); Rrootl_L=z2(NNd,cc);
Bfac_dayH=ones(NNd,cc,N_col); Bfac_weekH=ones(NNd,cc,N_col);
NPPI_H=z2(NNd,cc); TdpI_H=z2(NNd,cc);
Bfac_dayL=ones(NNd,cc,N_col); Bfac_weekL=ones(NNd,cc,N_col);
NPPI_L=z2(NNd,cc); TdpI_L=z2(NNd,cc);
Sr_H=z2(NNd,cc); Slf_H=z2(NNd,cc); Sfr_H=z2(NNd,cc); Swm_H=z2(NNd,cc); Sll_H=z2(NNd,cc);
Sr_L=z2(NNd,cc); Slf_L=z2(NNd,cc); Sfr_L=z2(NNd,cc); Swm_L=z2(NNd,cc); Sll_L=z2(NNd,cc);
LAIdead_H=z2(NNd,cc); LAIdead_L=z2(NNd,cc);
Rexmy_H=z3(NNd,cc,3); AgeDL_H=z2(NNd,cc);
Nreserve_H=z2(NNd,cc); Preserve_H=z2(NNd,cc); Kreserve_H=z2(NNd,cc);
Careserve_H=z2(NNd,cc); Mgreserve_H=z2(NNd,cc); Sireserve_H=z2(NNd,cc);
rNc_H=z2(NNd,cc); rPc_H=z2(NNd,cc); rKc_H=z2(NNd,cc);
rCac_H=z2(NNd,cc); rMgc_H=z2(NNd,cc); rSic_H=z2(NNd,cc);
Rexmy_L=z3(NNd,cc,3); AgeDL_L=z2(NNd,cc);
Nreserve_L=z2(NNd,cc); Preserve_L=z2(NNd,cc); Kreserve_L=z2(NNd,cc);
Careserve_L=z2(NNd,cc); Mgreserve_L=z2(NNd,cc); Sireserve_L=z2(NNd,cc);
rNc_L=z2(NNd,cc); rPc_L=z2(NNd,cc); rKc_L=z2(NNd,cc);
rCac_L=z2(NNd,cc); rMgc_L=z2(NNd,cc); rSic_L=z2(NNd,cc);
NBLeaf_H=z2(NNd,cc); PARI_H=z3(NNd,cc,3); NBLI_H=z2(NNd,cc);
NBLeaf_L=z2(NNd,cc); PARI_L=z3(NNd,cc,3); NBLI_L=z2(NNd,cc);
NupI_H=z3(NNd,cc,6); NupI_L=z3(NNd,cc,6); NavlI=z2(NNd,6);
RB_L=z3(NNd,cc,7); RB_H=z3(NNd,cc,7);
NuLit_H=z3(NNd,cc,6); NuLit_L=z3(NNd,cc,6);
BLit=z2(NNd,cc);
Nuptake_H=z2(NNd,cc); NH4uptake_H=z2(NNd,cc); NO3uptake_H=z2(NNd,cc);
Puptake_H=z2(NNd,cc); Kuptake_H=z2(NNd,cc);
Cauptake_H=z2(NNd,cc); Mguptake_H=z2(NNd,cc); Siuptake_H=z2(NNd,cc);
FNC_H=z2(NNd,cc);
TexC_H=z2(NNd,cc); TexN_H=z2(NNd,cc); TexP_H=z2(NNd,cc); TexK_H=z2(NNd,cc);
TNIT_H=z2(NNd,cc); TPHO_H=z2(NNd,cc); TPOT_H=z2(NNd,cc);
SupN_H=z2(NNd,cc); SupP_H=z2(NNd,cc); SupK_H=z2(NNd,cc);
TexCa_H=z2(NNd,cc); TexMg_H=z2(NNd,cc); TexSi_H=z2(NNd,cc);
TCAL_H=z2(NNd,cc); TMAG_H=z2(NNd,cc); TSIL_H=z2(NNd,cc);
SupCa_H=z2(NNd,cc); SupMg_H=z2(NNd,cc); SupSi_H=z2(NNd,cc);
Nuptake_L=z2(NNd,cc); NH4uptake_L=z2(NNd,cc); NO3uptake_L=z2(NNd,cc);
Puptake_L=z2(NNd,cc); Kuptake_L=z2(NNd,cc);
Cauptake_L=z2(NNd,cc); Mguptake_L=z2(NNd,cc); Siuptake_L=z2(NNd,cc);
FNC_L=z2(NNd,cc);
TexC_L=z2(NNd,cc); TexN_L=z2(NNd,cc); TexP_L=z2(NNd,cc); TexK_L=z2(NNd,cc);
TNIT_L=z2(NNd,cc); TPHO_L=z2(NNd,cc); TPOT_L=z2(NNd,cc);
SupN_L=z2(NNd,cc); SupP_L=z2(NNd,cc); SupK_L=z2(NNd,cc);
TexCa_L=z2(NNd,cc); TexMg_L=z2(NNd,cc); TexSi_L=z2(NNd,cc);
TCAL_L=z2(NNd,cc); TMAG_L=z2(NNd,cc); TSIL_L=z2(NNd,cc);
SupCa_L=z2(NNd,cc); SupMg_L=z2(NNd,cc); SupSi_L=z2(NNd,cc);
ISOIL_H=z3(NNd,cc,27); ISOIL_L=z3(NNd,cc,27);
P=z2(NNd,79);
R_litter=z2(NNd,1); R_microbe=z2(NNd,1); R_litter_sur=z2(NNd,1); R_ew=z2(NNd,1);
VOL=z2(NNd,1); N2flx=z2(NNd,1); N2Oflx=z2(NNd,1); ...
    NITflx=z2(NNd,1); DENITflx=z2(NNd,1); N2Oflx_nit=z2(NNd,1); ...
    FertN_d=z2(NNd,1); DepN_d=z2(NNd,1); PH_d=z2(NNd,1); ...
    PH_pore_d=z2(NNd,1); PH_bulk_d=z2(NNd,1);
Dep_d=z2(NNd,6); %% the day's deposition [N P K Ca Mg Si] (Deposition_Day)
Min_N=z2(NNd,1); Min_P=z2(NNd,1); R_bacteria=z2(NNd,1);
RmycAM=z2(NNd,1); RmycEM=z2(NNd,1); Prod_B=z2(NNd,1); Prod_F=z2(NNd,1);
BfixN=z2(NNd,1); LitFirEmi=z2(NNd,2);
LEAK_NH4=z2(NNd,1); LEAK_NO3=z2(NNd,1); LEAK_P=z2(NNd,1); LEAK_K=z2(NNd,1);
LEAK_DOC=z2(NNd,1); LEAK_DON=z2(NNd,1); LEAK_DOP=z2(NNd,1);
LEAK_Ca=z2(NNd,1); LEAK_Mg=z2(NNd,1); LEAK_Si=z2(NNd,1);
RexmyI=z2(NNd,3);
AgrHarNut=z2(NNd,6);
SeedInput=z2(NNd,7); %% [C N P K Ca Mg Si], [g element m-2 ground]
Ccrown_t=ones(NNd,cc,N_col);
AgePl_H=z2(NNd,cc); AgePl_L=z2(NNd,cc); Tden_H=z2(NNd,cc); Tden_L=z2(NNd,cc);
TBio_Ht=z2(NNd,cc); TBio_Lt=z2(NNd,cc);
q_runon=z2(NN,1); Qi_in=z2(NN,ms);
ManIH=zeros(cc,N_col); ManIL=zeros(cc,N_col);

%% Per-column scalars carried between steps
Ccrown_c = zeros(N_col,cc);  Cbare_c = zeros(N_col,1);
ZR95_H_c = zeros(N_col,cc);  ZR95_L_c = zeros(N_col,cc);
RfH_Zs_c = cell(N_col,1);    RfL_Zs_c = cell(N_col,1);
PH_c     = zeros(N_col,1);   PH_bulk_c = zeros(N_col,1); Tstm0_c = zeros(N_col,1);
Bam_c    = zeros(N_col,1);   Bem_c    = zeros(N_col,1);
snow_alb_c = cell(N_col,1);
B_IO_c   = cell(N_col,1);
TBio_H_c = zeros(N_col,cc);  TBio_L_c = zeros(N_col,cc);

for c = 1:N_col
    p = PAR{c};
    Ccrown_c(c,:) = p.Ccrown;   Cbare_c(c) = p.Cbare;
    ZR95_H_c(c,:) = p.ZR95_H;   ZR95_L_c(c,:) = p.ZR95_L;
    RfH_Zs_c{c} = p.RfH_Zs;     RfL_Zs_c{c} = p.RfL_Zs;
    %% Cold-start pH seed. p.PH is the OBSERVED BULK-SOIL pH from the
    %% parameter file, the quantity BIOGEO_UNIT consumes, so seeding
    %% PH_bulk_c from it is correct. The pore-water pH the ERW chemistry
    %% starts from is only known once ERW_MOD_PARAM_excel has inverted that
    %% observation (bulk_pH_params.initial_pH_mode = 'bulk_observed',
    %% Invert_Bulk_pH), which happens inside the first ERW_Daily_Step_2D
    %% call; the daily records written before that call are backfilled
    %% from the ERW initial condition there ("cold-start backfill" below).
    %% Before the first ERW call PH_c is only recorded (PH_pore_d), never
    %% consumed, so seeding it from the bulk value is harmless.
    PH_c(c) = p.PH(1);          PH_bulk_c(c) = p.PH(1); snow_alb_c{c} = p.snow_alb;
    if isfield(p,'ERW_biogeochem') && isfield(p.ERW_biogeochem,'pH_bulk') && ...
            isfinite(p.ERW_biogeochem.pH_bulk(1))
        PH_bulk_c(c) = p.ERW_biogeochem.pH_bulk(1);
    end
    B_IO_c{c} = p.B_IO;
    TBio_H_c(c,:) = p.TBio_H;   TBio_L_c(c,:) = p.TBio_L;

    %% Seed row 1 of every state array from this column's parameter file
    O(1,:,c) = p.O(1,:);  V(1,:,c) = p.V(1,:);
    Ts(1,c) = p.Ts(1);  Tstm0_c(c) = p.Ts(1);
    Tdp(1,:,c) = p.Tdp(1,:);  Tdamp(1,c) = p.Tdamp(1);
    Tdp_H(1,:,c) = p.Tdp_H(1,:);  Tdp_L(1,:,c) = p.Tdp_L(1,:);
    TdpI_H(1,:,c) = p.TdpI_H(1,:); TdpI_L(1,:,c) = p.TdpI_L(1,:);
    SWE(1,c)=p.SWE(1); SND(1,c)=p.SND(1); ros(1,c)=p.ros(1);
    In_SWE(1,c)=p.In_SWE(1); SP_wc(1,c)=p.SP_wc(1); t_sls(1,c)=p.t_sls(1);
    e_sno(1,c)=p.e_sno(1); tau_sno(1,c)=p.tau_sno(1); EK(1,c)=p.EK(1);
    WAT(1,c)=p.WAT(1); ICE(1,c)=p.ICE(1); IP_wc(1,c)=p.IP_wc(1);
    ICE_D(1,c)=p.ICE_D(1); FROCK(1,c)=p.FROCK(1); Ws_under(1,c)=p.Ws_under(1);
    In_L(1,:,c)=p.In_L(1,:); In_H(1,:,c)=p.In_H(1,:);
    In_urb(1,c)=p.In_urb(1); In_rock(1,c)=p.In_rock(1); In_Litter(1,c)=p.In_Litter(1);
    Ci_sunL(1,:,c)=p.Ci_sunL(1,:); Ci_sunH(1,:,c)=p.Ci_sunH(1,:);
    Ci_shdL(1,:,c)=p.Ci_shdL(1,:); Ci_shdH(1,:,c)=p.Ci_shdH(1,:);
    LAI_H(1,:,c)=p.LAI_H(1,:); LAI_L(1,:,c)=p.LAI_L(1,:);
    B_H(1,:,:,c)=p.B_H(1,:,:);  B_L(1,:,:,c)=p.B_L(1,:,:);
    Rrootl_H(1,:,c)=p.Rrootl_H(1,:); Rrootl_L(1,:,c)=p.Rrootl_L(1,:);
    PHE_S_H(1,:,c)=p.PHE_S_H(1,:); PHE_S_L(1,:,c)=p.PHE_S_L(1,:);
    dflo_H(1,:,c)=p.dflo_H(1,:); dflo_L(1,:,c)=p.dflo_L(1,:);
    AgeL_H(1,:,c)=p.AgeL_H(1,:); AgeL_L(1,:,c)=p.AgeL_L(1,:);
    e_rel_H(1,:,c)=p.e_rel_H(1,:); e_rel_L(1,:,c)=p.e_rel_L(1,:);
    hc_H(1,:,c)=p.hc_H(1,:); hc_L(1,:,c)=p.hc_L(1,:);
    SAI_H(1,:,c)=p.SAI_H(1,:); SAI_L(1,:,c)=p.SAI_L(1,:);
    Nreserve_H(1,:,c)=p.Nreserve_H(1,:); Preserve_H(1,:,c)=p.Preserve_H(1,:);
    Kreserve_H(1,:,c)=p.Kreserve_H(1,:);
    Careserve_H(1,:,c)=p.Careserve_H(1,:); Mgreserve_H(1,:,c)=p.Mgreserve_H(1,:);
    Sireserve_H(1,:,c)=p.Sireserve_H(1,:);
    Nreserve_L(1,:,c)=p.Nreserve_L(1,:); Preserve_L(1,:,c)=p.Preserve_L(1,:);
    Kreserve_L(1,:,c)=p.Kreserve_L(1,:);
    Careserve_L(1,:,c)=p.Careserve_L(1,:); Mgreserve_L(1,:,c)=p.Mgreserve_L(1,:);
    Sireserve_L(1,:,c)=p.Sireserve_L(1,:);
    FNC_H(1,:,c)=p.FNC_H(1,:); FNC_L(1,:,c)=p.FNC_L(1,:);
    NupI_H(1,:,:,c)=p.NupI_H(1,:,:); NupI_L(1,:,:,c)=p.NupI_L(1,:,:);
    RexmyI(1,:,c)=p.RexmyI(1,:); NavlI(1,:,c)=p.NavlI(1,:);
    Vx_H(1,:,c)=p.Vx_H; Vl_H(1,:,c)=p.Vl_H;
    Vx_L(1,:,c)=p.Vx_L; Vl_L(1,:,c)=p.Vl_L;
    BLit(1,:,c)=p.BLit(1);
    if OPT_SoilBiogeochemistry == 1
        P(1,:,c) = p.P(1,:);
        PH_d(1,c) = PH_bulk_c(c);
        PH_pore_d(1,c) = p.PH(1);
        PH_bulk_d(1,c) = PH_bulk_c(c);
    end
    Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,p.Cbare,p.Ccrown);
end
CcrownFIX = Ccrown_c;
[~,~,lateral_species,lateral_pool_idx,lateral_fraction_fields] = Biogeo_Lateral_Pools( ...
    reshape(P(1,:,:),79,N_col).');
N_lateral_species = numel(lateral_species);
i_lateral_NO3 = find(strcmp(lateral_species,'NO3'),1);
lateral_fraction = zeros(N_col,N_lateral_species);
BiogeoPar_lateral = Biogeochemistry_Parameter(1);
for s = 1:N_lateral_species
    lateral_fraction(:,s) = BiogeoPar_lateral.(lateral_fraction_fields{s});
end
if any(~isfinite(lateral_fraction(:))) || any(lateral_fraction(:) < 0) || ...
        any(lateral_fraction(:) > 1)
    error('MAIN_FRAME_BG_ERW_2D:BadMobileFraction', ...
        'Biogeo_Leakage mobility fractions must be finite and within [0,1].');
end
%% Flux tracking for the field spin-up (Lian et al. 2025, Fig. 3c): the
%% lateral exchange of the cells named in LATERAL_TRACK.cells is recorded
%% (Lateral_Track_Record), so the spin-up can replay each of them in a
%% plot-scale column.  No other run defines LATERAL_TRACK.
lateral_track_on = N_col > 1 && exist('LATERAL_TRACK','var') && ...
    isstruct(LATERAL_TRACK) && isfield(LATERAL_TRACK,'cells') && ...
    ~isempty(LATERAL_TRACK.cells);
if lateral_track_on
    lateral_track = Lateral_Track_Record('init',[],LATERAL_TRACK.cells,N_col, ...
        NN,NNd,ms,lateral_species,lateral_pool_idx,lateral_fraction,Date);
end

%% ========================================================================
%% NUMERICAL OPTIONS  (identical to the single-column frame)
%% ========================================================================
Opt_CR = optimset('TolFun',1);
Opt_ST = optimset('TolFun',0.1);
OPT_SM  = odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD  = odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',0.02);
OPT_PH  = odeset('AbsTol',0.01);
OPT_FR_SOIL = 1; OPT_VegSnow = 1; OPT_SoilTemp = 1;
if ~exist('OPT_VCA','var');  OPT_VCA = 0; OPT_ALLOME = 0; end
if ~exist('OPT_PlantHydr','var');      OPT_PlantHydr = 0;      end
if ~exist('OPT_EnvLimitGrowth','var'); OPT_EnvLimitGrowth = 1; end
%%% v1.0: a column whose vegetation file sets OPT_EnvLimitGrowth keeps it.
OPT_EnvLimitGrowth_c = repmat(OPT_EnvLimitGrowth,1,N_col);
for c = 1:N_col
    if isfield(PAR{c},'OPT_EnvLimitGrowth'); OPT_EnvLimitGrowth_c(c) = PAR{c}.OPT_EnvLimitGrowth; end
end
OPT_WET = 0;

%% ERW per-column state
%% Cell array, not a struct array: each column's ERW state gains fields
%% after its first call, which a struct array cannot absorb elementwise.
ERW = cell(N_col,1);
for c = 1:N_col; ERW{c} = struct('init',false,'has_update',false); end
conv_mol = 1e6; conv_Al = 1e3;
Tdp_bio_c = zeros(N_col,1); VSUM_c = zeros(N_col,1); VTSUM_c = zeros(N_col,1);
lateral_diag = struct('water_buf_to_main', zeros(NNd,1), ...
                      'mass_buf_to_main',  zeros(NNd,N_lateral_species), ...
                      'NO3_buf_to_main',   zeros(NNd,1), ...
                      'heat_buf_to_main_J',zeros(NNd,1), ...
                      'heat_residual_J',   zeros(NNd,1), ...
                      'heat_residual_relative',zeros(NNd,1), ...
                      'heat_sub_steps_max',zeros(NNd,1), ...
                      'erw_mass_buf_to_main_umol',zeros(NNd,4), ...
                      'erw_CO2_buf_to_main_umol',zeros(NNd,1), ...
                      'erw_mass_residual_umol',zeros(NNd,4), ...
                      'erw_mass_abs_residual_umol',zeros(NNd,4), ...
                      'erw_mass_throughput_umol',zeros(NNd,4), ...
                      'sub_steps_max',     zeros(NNd,1), ...
                      'species_names',     {lateral_species}, ...
                      'erw_species_names', {{'Na','Al','An','DIC'}});

%% Domain-boundary accounting. The subsurface lateral perimeter is closed;
%% bottom leakage leaves through Leakage_Bottom.  Surface runoff leaves the
%% plot at each cell (flat, no topography) OR is routed downslope and
%% leaves only at the topographic outlet(s) (Surface_Runoff_Routing).
exit_water_vert_m3 = zeros(N_col,1);
exit_water_surface_m3 = zeros(N_col,1);
exit_mass_surface_g = zeros(N_lateral_species,1);   %% solute exported at outlets
exit_mass_surface_erw_umol = zeros(4,1);            %% Na, Al, An, DIC
exit_mass_surface_erw_field_umol = zeros(4,1);
erw_pond_umol = zeros(N_col,4);                     %% ponded ERW run-on load
erw_pond_field_umol = zeros(N_col,4);               %% its field-origin part
%% ------------------------------------------------------------------
%% SOLAR WINDOW FOR THE TOPOGRAPHIC BEAM PROJECTION
%%   SetSunVariables averages the sun over [HR-t_bef, HR+t_aft].  The run
%%   scripts pass t_bef=-0.67, t_aft=1.67, i.e. [HR+0.67, HR+1.67]; the
%%   ERA5-Land forcing labels each record with the END of its accumulation
%%   hour, so record HR covers [HR-1, HR].  The two differ by ~1.7 h, which
%%   is ~25 deg of azimuth at midday -- enough to reverse the east/west
%%   aspect contrast that the topographic projection exists to resolve.
%%   The 1-D energy balance is far less sensitive to this and its
%%   calibration is a separate question, so the global t_bef/t_aft is left
%%   alone and only the per-cell BEAM PROJECTION uses the forcing's own
%%   convention.  Override with TOPO_SOLAR_WINDOW = [t_bef t_aft].
if ~exist('TOPO_SOLAR_WINDOW','var') || isempty(TOPO_SOLAR_WINDOW)
    TOPO_SOLAR_WINDOW = [1 0];          %% hour-ending, matching ERA5-Land
end
topo_t_bef = TOPO_SOLAR_WINDOW(1);
topo_t_aft = TOPO_SOLAR_WINDOW(2);
if topography_attached
    %% Report how well the chosen window tracks the observed shortwave, on
    %% a subsample, so a forcing with a different time convention is
    %% visible in the run log instead of silently rotating every aspect.
    probe = 1:max(1,floor(NN/500)):NN;
    Rsw_probe = SAB1(probe)+SAB2(probe)+SAD1(probe)+SAD2(probe);
    sh_topo = zeros(numel(probe),1); sh_run = zeros(numel(probe),1);
    for pk = 1:numel(probe)
        sh_topo(pk) = max(0,sin(SetSunVariables(Datam(probe(pk),:), ...
            DeltaGMT,Lon,Lat,topo_t_bef,topo_t_aft)));
        sh_run(pk) = max(0,sin(SetSunVariables(Datam(probe(pk),:), ...
            DeltaGMT,Lon,Lat,t_bef,t_aft)));
    end
    if std(Rsw_probe(:)) > 0 && std(sh_topo) > 0
        c_topo = corr(sh_topo,Rsw_probe(:));
        c_run  = corr(sh_run,Rsw_probe(:));
        fprintf(['  topographic solar window [HR%+.2f, HR%+.2f]: sun-vs-' ...
            'shortwave correlation %.3f (run-script window %.3f)\n'], ...
            -topo_t_bef,topo_t_aft,c_topo,c_run);
        if c_topo < c_run
            warning('MAIN_FRAME_BG_ERW_2D:SolarWindowMismatch', ...
                ['The forcing tracks the run-script solar window (%.3f) ' ...
                 'better than the topographic one (%.3f); set ' ...
                 'TOPO_SOLAR_WINDOW to match this forcing time convention ' ...
                 'or the per-cell beam projection will be rotated.'], ...
                c_run,c_topo);
        end
    end
end
%% Polygon-overlap fraction per cell, needed inside the hourly loop to
%% carry field-origin provenance through the overland router.
[~,~,field_fraction] = Field_Aggregation_Areas(topo);
%% Overland routing carries this step's ponded runoff to the NEXT step's
%% re-infiltration opportunity (q_runon).  T&C's operator split is
%% 1-step-lagged by design; overland travel times over a farm field are
%% far below one hydrologic step (test_lateral_routing #14).
q_runon_next = zeros(N_col,1);
%% Solute that arrived with routed run-on and is still ponded on the
%% surface, waiting for the next vertical step to infiltrate it or run it
%% off again [gX/m2].  Crediting it straight to the soil pool would strand
%% it: the water it rode in on has not infiltrated yet, and the part that
%% leaves again as infiltration-excess runoff must take its share with it.
pond_mass_area = zeros(N_col,N_lateral_species);
%% Field-origin provenance of the pond store: the water depth [mm] and the
%% solute [gX/m2] that came from inside the polygon.  Carried with the
%% store so a stop-over on a partial perimeter cell does not re-label the
%% load with that cell's own overlap fraction.
pond_field_water_mm = zeros(N_col,1);
pond_field_mass_area = zeros(N_col,N_lateral_species);
exit_water_field_surface_m3 = 0;    %% field-origin part of the outlet export
exit_mass_surface_field_g = zeros(N_lateral_species,1);
lateral_diag.surface_water_out_m3 = zeros(NNd,1);
lateral_diag.surface_mass_out_g   = zeros(NNd,N_lateral_species);
lateral_diag.surface_reinfilt_m3  = zeros(NNd,1);
%% Subsurface lateral throughput per day: the sum over faces, layers and
%% hours of |Q_edge| [m3].  The buffer/main ledgers above are zero on a
%% whole-field polygon grid (no buffer), so this is the quantity that
%% says how much water the topography actually moved between cells.
lateral_diag.water_edge_abs_m3 = zeros(NNd,1);
%% Net redistribution per day [m3]: half the sum over cells of each cell's
%% net lateral gain or loss, i.e. the water that actually changed cells,
%% counted once.  The throughput above also counts exchange back and forth.
lateral_diag.water_net_m3 = zeros(NNd,1);

j = 1;
tic;

%% ========================================================================
%% CHECKPOINT / RESUME
%% ========================================================================
%% A full column-year costs tens of minutes, so a run that aborts late
%% (or a laptop that hibernates) should not have to start from day 1.
%% Every ERW_ckpt_minutes of wall time the whole workspace is written to
%% one rolling file, atomically (temp + rename), so an interrupted save
%% can never leave a half-written checkpoint behind.
%%
%% Resuming is only safe when the checkpoint was produced by the SAME
%% code and the same run configuration: a checkpoint carries the
%% trajectory that the old solver produced, and replaying it under
%% changed chemistry would silently splice two model versions together.
%% The fingerprint below covers the run identity, the horizon, and the
%% mtime+size of the files that decide the chemistry.  On a mismatch the
%% run starts clean and says so.  Set ERW_ckpt_force_resume = true to
%% override deliberately (useful when debugging a late-timestep failure
%% where the earlier trajectory is known to be unaffected).
if ~exist('ERW_ckpt_minutes','var') || isempty(ERW_ckpt_minutes)
    ERW_ckpt_minutes = 10;
end
if ~exist('ERW_ckpt_file','var') || isempty(ERW_ckpt_file)
    ERW_ckpt_file = '';
end
if ~exist('ERW_ckpt_force_resume','var') || isempty(ERW_ckpt_force_resume)
    ERW_ckpt_force_resume = false;
end
if ~exist('ERW_cell_input_dependencies','var') || isempty(ERW_cell_input_dependencies)
    ERW_cell_input_dependencies = {};
end
ERW_runtime_parts = cell(1,N_col);
for c = 1:N_col
    appelo_alpha_c = [];
    appelo_logKH_c = [];
    knmaxd_c = [];
    cec_dpH_c = [];
    ERW_runtime_param_c = ERW_PARAM_input;
    MOD_runtime_param_c = MOD_PARAM_input;
    if isfield(ColCfg(c),'ERW_PARAM')
        ERW_runtime_param_c = ColCfg(c).ERW_PARAM;
    end
    if isfield(ColCfg(c),'MOD_PARAM') && istable(ColCfg(c).MOD_PARAM) && ...
            height(ColCfg(c).MOD_PARAM) == 1
        MOD_runtime_param_c = ColCfg(c).MOD_PARAM;
    end
    if isfield(ColCfg(c),'appelo_alpha')
        appelo_alpha_c = ColCfg(c).appelo_alpha;
    end
    if isfield(ColCfg(c),'appelo_logKH')
        appelo_logKH_c = ColCfg(c).appelo_logKH;
    end
    if isfield(ColCfg(c),'KNmaxd')
        knmaxd_c = ColCfg(c).KNmaxd;
    end
    if isfield(ColCfg(c),'cec_dpH')
        cec_dpH_c = ColCfg(c).cec_dpH;
    end
    last_app_c = [];
    if isfield(ColCfg(c),'last_application_year')
        last_app_c = ColCfg(c).last_application_year;
    end
    %% Replaced biogeochemistry literals and the saturation exponent are
    %% part of the run's identity: a checkpoint of one member must not
    %% resume another.
    bgovr_c = '';
    if isfield(ColCfg(c),'biogeo_overrides') && isstruct(ColCfg(c).biogeo_overrides)
        bgovr_names = sort(fieldnames(ColCfg(c).biogeo_overrides));
        bgovr_pairs = cellfun(@(n) sprintf('%s=%s',n, ...
            mat2str(ColCfg(c).biogeo_overrides.(n),17)),bgovr_names, ...
            'UniformOutput',false);
        bgovr_c = strjoin(bgovr_pairs,',');
    end
    sexp_c = [];
    if isfield(ColCfg(c),'saturation_exponent')
        sexp_c = ColCfg(c).saturation_exponent;
    end
    pft_override_c = '';
    if isfield(ColCfg(c),'PFT_OVERRIDE') && ~isempty(ColCfg(c).PFT_OVERRIDE)
        pft_override_c = char(ColCfg(c).PFT_OVERRIDE);
    end
    agronomic_identity_c = '';
    if isfield(ColCfg(c),'agronomic_identity')
        agronomic_identity_c = char(ColCfg(c).agronomic_identity);
    end
    ERW_runtime_parts{c} = sprintf( ...
        ['c%d:column_param=%s;pft=%s;apply_rock=%d;mod_param=%s;' ...
         'erw_param=%s;agronomic=%s;alpha=%s;logKH=%s;pcla=%s;psan=%s;' ...
         'knmaxd=%s;kbot=%s;cecdpH=%s;lastapp=%s;bgovr=%s;sexp=%s'], ...
        c,char(ColCfg(c).PARAM_IC),pft_override_c,logical(ColCfg(c).apply_rock), ...
        ERW_Runtime_Parameter_Identity(MOD_runtime_param_c), ...
        ERW_Runtime_Parameter_Identity(ERW_runtime_param_c), ...
        agronomic_identity_c, ...
        mat2str(appelo_alpha_c,17), ...
        mat2str(appelo_logKH_c,17),mat2str(PAR{c}.Pcla,17), ...
        mat2str(PAR{c}.Psan,17),mat2str(knmaxd_c,17), ...
        mat2str(PAR{c}.Kbot,17), ...
        mat2str(cec_dpH_c,17),mat2str(last_app_c,17),bgovr_c, ...
        mat2str(sexp_c,17));
    if isfield(PAR{c},'rain_Na_flux_umol_m2_d')
        ERW_runtime_parts{c} = [ERW_runtime_parts{c} sprintf(';rain=%s', ...
            mat2str([PAR{c}.rain_pH_vwm PAR{c}.rain_Na_flux_umol_m2_d ...
            PAR{c}.rain_strong_anion_flux_umolc_m2_d PAR{c}.seed_water_loss_l_m2_d],17))];
    end
end
ERW_forcing_identity = sprintf('date_start=%.17g;date_end=%.17g',Date(1),Date(end));
topology_overlap_area = topo.cell_area_m2;
if isfield(topo,'cell_overlap_area_m2')
    topology_overlap_area = topo.cell_overlap_area_m2;
end
topology_rotation = [];
topology_origin = [];
if isfield(topo,'grid_rotation_degrees'); topology_rotation = topo.grid_rotation_degrees; end
if isfield(topo,'projected_grid_origin_m'); topology_origin = topo.projected_grid_origin_m; end
ERW_topology_identity = sprintf([ ...
    'topology_areas=%s;topology_overlap_areas=%s;' ...
    'topology_rotation=%s;topology_origin=%s;' ...
    'topology_edge_a=%s;topology_edge_b=%s;' ...
    'topology_edge_width=%s;topology_edge_distance=%s;' ...
    'topology_edge_distance_a=%s;topology_edge_distance_b=%s;' ...
    'topology_edge_active=%s;topology_outer_width=%s'], ...
    mat2str(topo.cell_area_m2,17),mat2str(topology_overlap_area,17), ...
    mat2str(topology_rotation,17),mat2str(topology_origin,17), ...
    mat2str(topo.edge_a,17), ...
    mat2str(topo.edge_b,17),mat2str(topo.edge_face_width_m,17), ...
    mat2str(topo.edge_face_distance_m,17), ...
    mat2str(topo.edge_face_distance_a_m,17), ...
    mat2str(topo.edge_face_distance_b_m,17), ...
    mat2str(topo.edge_active,17),mat2str(topo.outer_face_width_m,17));
%% Topographic attributes only extend the identity when a DEM has been
%% attached, so a flat run's fingerprint is unchanged from before.
if isfield(topo,'topography_attached') && topo.topography_attached
    %% Attach_Topography_To_Topology writes all of these together, so a
    %% missing field here means a malformed topology and should error.
    topo_field = @(name) mat2str(reshape(double(topo.(name)),[],1),17);
    ERW_topology_identity = sprintf(['%s;' ...
        'topology_source=%s;topology_cell_z=%s;topology_cell_slope=%s;' ...
        'topology_cell_aspect=%s;topology_cell_svf=%s;topology_cell_twi=%s;' ...
        'topology_cell_curv_prof=%s;topology_cell_curv_plan=%s;' ...
        'topology_edge_dz=%s;topology_edge_surface_distance=%s;' ...
        'topology_overland_ab=%s;topology_overland_ba=%s;topology_is_outlet=%s'], ...
        ERW_topology_identity, ...
        char(string(topo.topography_source)), ...
        topo_field('cell_z_m'),topo_field('cell_slope'), ...
        topo_field('cell_aspect_deg'),topo_field('cell_svf'), ...
        topo_field('cell_twi'),topo_field('cell_curv_prof'), ...
        topo_field('cell_curv_plan'),topo_field('edge_dz_m'), ...
        topo_field('edge_face_distance_surface_m'), ...
        topo_field('edge_overland_weight_ab'), ...
        topo_field('edge_overland_weight_ba'),topo_field('cell_is_outlet'));
end
ERW_numerical_identity = strtrim(getenv('TC_NUMERICAL_IDENTITY'));
if isempty(ERW_numerical_identity); ERW_numerical_identity = 'unset'; end
%% The static TC_NUMERICAL_IDENTITY string is what the submission script
%% asked for, not what MATLAB actually loaded.  A new module build can pass
%% campaign validation while an older checkpoint still resumes, splicing
%% two numerical runtimes into one trajectory, so the identity carries the
%% interpreter, platform, BLAS and thread count that are really in force.
ERW_execution_identity = sprintf( ...
    ['parallel_workers=%d;numerical_identity=%s;matlab=%s;computer=%s;' ...
     'blas=%s;max_threads=%d'], ...
    parallel_workers,ERW_numerical_identity,version,computer, ...
    regexprep(strtrim(evalc('version -blas')),'\s+',' '),maxNumCompThreads);
ERW_ckpt_runtime_identity = strjoin( ...
    [{ERW_forcing_identity,ERW_topology_identity,ERW_execution_identity} ...
     ERW_runtime_parts],'|');
if lateral_track_on
    %% A tracking run records its own set of cells, with the ERW chemistry
    %% off; it may resume no other run's checkpoint, nor another its own.
    ERW_ckpt_runtime_identity = sprintf('%s|lateral_track=%s;opt_erw=%d', ...
        ERW_ckpt_runtime_identity,mat2str(LATERAL_TRACK.cells(:).'),OPT_ERW);
end
ERW_ckpt_fingerprint = ERW_Checkpoint_Fingerprint(NN, N_col, ...
    ERW_ckpt_runtime_identity, ...
    fullfile(cd, '*.m'), ...
    fullfile(cd, 'Parameters', '*.m'), ...
    fullfile(cd, 'BIOGEOCHEMISTRY_DYNAMIC3.m'), ...
    fullfile(cd, 'BIOGEO_UNIT.m'), ...
    fullfile(cd, 'Biogeochemistry_Parameter.m'), ...
    fullfile(cd, 'Biogeo_Cation_Partition_Rates.m'), ...
    fullfile(cd, 'Biogeo_Mobile_Water_Volume.m'), ...
    fullfile(cd, 'Biogeo_Nitrification_PH_Response.m'), ...
    fullfile(cd, 'CropHeightType.m'), ...
    fullfile(cd, 'Plant_Exports.m'), ...
    fullfile(cd, 'Plant_Reserve_Fraction.m'), ...
    fullfile(cd, 'Plant_Nutrient_Suppression.m'), ...
    fullfile(cd, 'Plant_Harvest_Reserve_Export.m'), ...
	fullfile(cd, 'Crop_Planting_Harvest.m'), ...
	fullfile(cd, 'Plant_Establishment_Element_Boundary.m'), ...
	fullfile(cd, 'Plant_Seed_Input_Accumulator.m'), ...
    fullfile(cd, 'ERW_function.m'), fullfile(cd, 'ERW_Daily_Step_2D.m'), ...
    fullfile(cd, 'ERW_Particle_History.m'), ...
    fullfile(cd, 'ERW_Initialize_TC_Cation_Pools.m'), ...
    fullfile(cd, 'ERW_Normalize_TC_Daily_Deltas.m'), ...
    fullfile(cd, 'ERW_Should_Apply_TC_Daily_Delta.m'), ...
    fullfile(cd, 'ERW_Return_Weathered_Phosphorus.m'), ...
    fullfile(cd, 'ERW_Appelo_Mass_Conserving_State.m'), ...
    fullfile(cd, 'ERW_Appelo_State_At_H.m'), ...
    fullfile(cd, 'ERW_CEC_At_H.m'), ...
    fullfile(cd, 'ERW_AlOH3_LogK.m'), ...
    fullfile(cd, 'ERW_Calcite_Kinetics.m'), ...
    fullfile(cd, 'ERW_Resolve_CEC.m'), ...
    fullfile(cd, 'ERW_Resolve_CEC_pH_Sensitivity.m'), ...
    fullfile(cd, 'ERW_Seed_Water_Loss.m'), ...
    fullfile(cd, 'ERW_Calcite_Reaction_Limit.m'), ...
    fullfile(cd, 'ERW_Feedstock_Pathway.m'), ...
    fullfile(cd, 'ERW_Particle_Mass_Fractions.m'), ...
    fullfile(cd, 'ERW_Mineral_DIC_Stoichiometry.m'), ...
    fullfile(cd, 'ERW_Mineral_Molar_Mass_Row.m'), ...
    fullfile(cd, 'ERW_Mineral_Saturation.m'), ...
    fullfile(cd, 'ERW_Total_Calcite_Solid.m'), ...
    fullfile(cd, 'ERW_pre_processing.m'), fullfile(cd, 'MAIN_FRAME_BG_ERW_2D.m'), ...
    fullfile(cd, 'Parameters', 'ERW_MOD_PARAM_excel.m'), ...
    fullfile(cd, 'Parameters', 'MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(cd, 'Inputs', 'ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(cd, 'Inputs', 'Data_Kantola.mat'), ...
    fullfile(cd, 'Inputs', 'Data_Kantola_BG.mat'), ...
    fullfile(cd, 'Inputs', 'Ca_Data.mat'), ...
    fullfile(cd, 'Inputs', 'All_deposition_data.mat'), ...
    fullfile(cd, 'Inputs', 'blueridge_CDF_data.json'), ...
    fullfile(cd, '..', 'run', 'Run_One_Subfield.m'), ...
    fullfile(cd, '..', 'run', 'Run_Matrix_12.m'), ...
    fullfile(cd, 'Resolve_Appelo_Constants.m'), ...
    fullfile(cd, 'Invert_Bulk_pH.m'), ...
    fullfile(cd, 'Porewater_to_Bulk_pH.m'), ...
    ERW_cell_input_dependencies{:});
i_resume = 2;
if ~isempty(ERW_ckpt_file) && isfile(ERW_ckpt_file)
    ckpt_probe = load(ERW_ckpt_file, 'ERW_ckpt_fingerprint', 'i');
    ckpt_ok = isfield(ckpt_probe,'ERW_ckpt_fingerprint') && ...
        strcmp(ckpt_probe.ERW_ckpt_fingerprint, ERW_ckpt_fingerprint);
    if ckpt_ok || ERW_ckpt_force_resume
        keep_ckpt_file = ERW_ckpt_file;
        keep_ckpt_minutes = ERW_ckpt_minutes;
        keep_ckpt_fp = ERW_ckpt_fingerprint;
        %% The caller's options describe this invocation, not the one that
        %% wrote the checkpoint.  The fingerprint pins the scientific
        %% configuration, but SaveResults and ResultDirectory are not in
        %% it, so restoring them wholesale would make a resumed run ignore
        %% where it was told to write and whether to write at all.
        if exist('opt','var'); keep_opt = opt; else; keep_opt = []; end
        %% The checkpoint carries the provenance of the segment that wrote
        %% it.  Restoring it wholesale would make the final MAT report the
        %% first Slurm job and node as if the run had never been
        %% interrupted, so the segments are accumulated instead.
        if exist('run_provenance','var'); keep_provenance = run_provenance;
        else;                             keep_provenance = [];
        end
        load(ERW_ckpt_file);
        ERW_ckpt_file = keep_ckpt_file;
        ERW_ckpt_minutes = keep_ckpt_minutes;
        ERW_ckpt_fingerprint = keep_ckpt_fp;
        if ~isempty(keep_opt); opt = keep_opt; end
        %% A checkpoint written before a ledger existed restores a
        %% lateral_diag without it; start that ledger empty rather than
        %% fail at the first lateral step.
        for ledger_name = {'water_edge_abs_m3','water_net_m3'}
            if exist('lateral_diag','var') && isstruct(lateral_diag) && ...
                    ~isfield(lateral_diag,ledger_name{1})
                lateral_diag.(ledger_name{1}) = zeros(NNd,1);
            end
        end
        if ~isempty(keep_provenance)
            prior_segments = {};
            if exist('run_provenance','var') && isstruct(run_provenance)
                if isfield(run_provenance,'execution_segments')
                    prior_segments = run_provenance.execution_segments;
                    run_provenance = rmfield(run_provenance, ...
                        'execution_segments');
                end
                prior_segments{end+1} = run_provenance;
            end
            run_provenance = keep_provenance;
            run_provenance.execution_segments = prior_segments;
            run_provenance.resumed_from_iteration = i;
        end
        i_resume = i + 1;
        if ckpt_ok
            fprintf(['Resuming from checkpoint at iter %d / %d (%.1f%%)\n'], ...
                i, NN, 100*i/NN);
        else
            warning('ERW:CheckpointCodeMismatch', ...
                ['Resuming from a checkpoint written by DIFFERENT code ' ...
                 '(ERW_ckpt_force_resume = true). The trajectory before ' ...
                 'iter %d came from the previous version.'], i);
        end
        tic;
    else
        fprintf(['Checkpoint found but code/config changed; starting clean.\n']);
    end
end
ERW_ckpt_last_toc = toc;

%% ========================================================================
%% TIME LOOP
%% ========================================================================
for i = i_resume:NN
    if (mod(i,1000) == 0) || (i == 2)
        fprintf('Iter %d / %d   (%.1f%%)\n', i, NN, 100*i/NN);
        %% Pool progress, printed to the run log only. Read-only: it
        %% touches no state and cannot alter the simulation. Column 1 is
        %% enough to spot a pool heading to zero or running away.
        if OPT_SoilBiogeochemistry == 1 && j > 1
            fprintf(['  POOLS  SOMC %.4g  POC %.4g  MOC %.4g  DOC %.4g  ' ...
                     'Bact %.4g  Fung %.4g  AM %.4g  SOMN %.4g  ' ...
                     'NH4 %.4g  NO3 %.4g  Pmin %.4g  Ca %.4g\n'], ...
                P(j-1,9,1), P(j-1,10,1), P(j-1,11,1), ...
                P(j-1,12,1)+P(j-1,13,1), P(j-1,18,1), P(j-1,19,1), ...
                P(j-1,20,1), P(j-1,26,1), P(j-1,31,1), P(j-1,32,1), ...
                P(j-1,43,1), P(j-1,60,1));
        end
    end
    pdind = max(1,i-24):i-1;
    h_S = NaN;

    %% ------------------------------------------------------------------
    %% PER-CELL TOPOGRAPHIC RADIATION
    %%   slope/aspect project the direct beam (Essery & Marks 2007;
    %%   the T&C spatial-radiation form used by T&C-BG-2D); the sky-view
    %%   factor scales the diffuse and adds a terrain-reflected term.
    %%   Flat / no topography attached -> every cell keeps the station
    %%   series, so this is a no-op there.
    %% ------------------------------------------------------------------
    SAB1_c = repmat(SAB1(i),N_col,1); SAB2_c = repmat(SAB2(i),N_col,1);
    SAD1_c = repmat(SAD1(i),N_col,1); SAD2_c = repmat(SAD2(i),N_col,1);
    PARB_c = repmat(PARB(i),N_col,1); PARD_c = repmat(PARD(i),N_col,1);
    SvF_c  = arrayfun(@(c) PAR{c}.SvF, (1:N_col).');
    if topography_attached
        SvF_c = topo.cell_svf(:);
        %% One definition with the field spin-up's replay of a single cell
        %% (Lateral_Replay_Prepare), so a replayed cell sees this radiation.
        [SAB1_c,SAB2_c,SAD1_c,SAD2_c,PARB_c,PARD_c] = Topographic_Radiation_Projection( ...
            SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i),Datam(i,:), ...
            DeltaGMT,Lon,Lat,topo_t_bef,topo_t_aft, ...
            topo.cell_slope,topo.cell_aspect_deg,topo.cell_svf);
    end

    %% ------------------------------------------------------------------
    %% DAILY BLOCK
    %% ------------------------------------------------------------------
    if (Datam(i,4) == 1)
        j = j+1;
        [jDay_i] = JULIAN_DAY(Datam(i,:));
        [~,~,~,~,~,L_day_j] = SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        jDay(j) = jDay_i; L_day(j) = L_day_j; %#ok<SAGROW>
        if lateral_track_on
            lateral_track = Lateral_Track_Record('day',lateral_track,j,Date(i));
        end

        for c = 1:N_col
            p = PAR{c};
            Ccrown = Ccrown_c(c,:); Cbare = Cbare_c(c);
            ZR95_H = ZR95_H_c(c,:); ZR95_L = ZR95_L_c(c,:);
            RfH_Zs = RfH_Zs_c{c};   RfL_Zs = RfL_Zs_c{c};
            B_IO_c_now = B_IO_c{c};

            %% Rotation fertiliser rule for cropped cells
            if isempty(char(ColCfg(c).PFT_OVERRIDE))
                if ismember(Datam(i,1),B_IO_c_now.FertN_years)
                    B_IO_c_now.FertN(81:100) = 1.01; %% 202 kg N ha-1 (Kantola 2022; 2011/2014 were 168, outside 2016-2020)
                else
                    B_IO_c_now.FertN = zeros(1,366);
                end
                B_IO_c{c} = B_IO_c_now;
            end

            [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM,f_anaer_bio,WFPS_bio] = Biogeo_environment( ...
                Tdp(pdind,:,c),O(pdind,:,c),V(pdind,:,c),p.Soil_Param,p.Phy,p.SPAR,p.Bio_Zs,p.Zbio);
            Tdp_bio_c(c) = Tdp_bio; VSUM_c(c) = VSUM; VTSUM_c(c) = VTSUM;

            if OPT_SoilBiogeochemistry == 1
                IS    = Ccrown*squeeze(ISOIL_L(j-1,:,:,c)) + Ccrown*squeeze(ISOIL_H(j-1,:,:,c));
                REXMY = Ccrown*squeeze(Rexmy_L(j-1,:,:,c)) + Ccrown*squeeze(Rexmy_H(j-1,:,:,c));
                FireA = 1*((sum(ManIH(:,c)==-5) + sum(ManIL(:,c)==-5)) > 0);
                if sum(Ccrown) > 0
                    OPT_BfixN = p.ParEx_L(Ccrown>0).bfix;
                else
                    OPT_BfixN = 0;
                end
                if ~(isfinite(PH_c(c)) && PH_c(c) > 0 && PH_c(c) < 14)
                    error('BIOGEO:InvalidPorePH', ...
                        'Column %d pore-water pH must be finite and within (0,14).',c);
                end
                if ~(isfinite(PH_bulk_c(c)) && PH_bulk_c(c) > 0 && PH_bulk_c(c) < 14)
                    error('BIOGEO:InvalidBulkPH', ...
                        'Column %d microbial bulk-soil pH must be finite and within (0,14).',c);
                end
                PH_d(j,c) = PH_bulk_c(c);
                PH_pore_d(j,c) = PH_c(c);
                PH_bulk_d(j,c) = PH_bulk_c(c);
                %% The day's deposition: constant, or wet on the day's
                %% precipitation plus the year's dry (Deposition_Day).
                [B_IO_bgc_day,Dep_d(j,:,c)] = Deposition_Day(B_IO_c_now, ...
                    Datam(i,:),sum(Pr(pdind))*dth);

                [P(j,:,c),LEAK_NH4(j,c),LEAK_NO3(j,c),LEAK_P(j,c),LEAK_K(j,c), ...
                 LEAK_DOC(j,c),LEAK_DON(j,c),LEAK_DOP(j,c), ...
                 LEAK_Ca(j,c),LEAK_Mg(j,c),LEAK_Si(j,c), ...
                 Nuptake_H(j,:,c),NH4uptake_H(j,:,c),NO3uptake_H(j,:,c), ...
                 Puptake_H(j,:,c),Kuptake_H(j,:,c), ...
                 Nuptake_L(j,:,c),NH4uptake_L(j,:,c),NO3uptake_L(j,:,c), ...
                 Puptake_L(j,:,c),Kuptake_L(j,:,c), ...
                 Cauptake_H(j,:,c),Mguptake_H(j,:,c),Siuptake_H(j,:,c), ...
                 Cauptake_L(j,:,c),Mguptake_L(j,:,c),Siuptake_L(j,:,c), ...
                 RexmyI(j,:,c), ...
                 R_litter(j,c),R_microbe(j,c),R_litter_sur(j,c),R_ew(j,c), ...
                 VOL(j,c),N2flx(j,c),Min_N(j,c),Min_P(j,c), ...
                 R_bacteria(j,c),RmycAM(j,c),RmycEM(j,c),Prod_B(j,c),Prod_F(j,c), ...
                 BfixN(j,c),NavlI(j,:,c),LitFirEmi(j,:,c),N2Oflx(j,c), ...
                 FertN_d(j,c),DepN_d(j,c),NITflx(j,c),DENITflx(j,c),N2Oflx_nit(j,c)] = ...
                 BIOGEO_UNIT(P(j-1,:,c),IS,p.Zbio,sum(p.Bio_Zs.*p.rsd), ...
                    PH_bulk_c(c),Tdp_bio,mean(Ta(pdind)),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM, ...
                    Ccrown,p.Bio_Zs,RfH_Zs,RfL_Zs, ...
                    sum(Lk(pdind,1,c)),sum(T_H(pdind,:,c),1),sum(T_L(pdind,:,c),1), ...
                    B_H(j-1,:,3,c),B_L(j-1,:,3,c),LAI_H(j-1,:,c),LAI_L(j-1,:,c), ...
                    SupN_H(j-1,:,c),SupP_H(j-1,:,c),SupK_H(j-1,:,c), ...
                    SupN_L(j-1,:,c),SupP_L(j-1,:,c),SupK_L(j-1,:,c), ...
                    SupCa_H(j-1,:,c),SupMg_H(j-1,:,c),SupSi_H(j-1,:,c), ...
                    SupCa_L(j-1,:,c),SupMg_L(j-1,:,c),SupSi_L(j-1,:,c), ...
                    REXMY,RexmyI(j-1,:,c),p.ExEM,NavlI(j-1,:,c),p.Pcla,p.Psan, ...
                    B_IO_bgc_day,Datam(i,:),FireA,0, ...
                    OPT_ERW,OPT_BfixN,f_anaer_bio,WFPS_bio, ...
                    biogeo_exchange_scheme{c},biogeo_KNmaxd{c}, ...
                    biogeo_overrides{c});

                BLit(j,:,c) = 0.002*sum(P(j,1:5,c))*Ccrown';
                Bam_c(c) = P(j,20,c);
                Bem_c(c) = P(j,21,c);
            else
                BLit(j,:,c) = 0;
                NavlI(j,:,c) = 0; P(j,:,c) = 0;
                Bam_c(c) = 0; Bem_c(c) = 0;
            end

            %% ---- Vegetation, per PFT slot -----------------------------
            for k = 1:length(Ccrown)
                if ZR95_H(k) > 0
                    [LAI_H(j,k,c),B_H(j,k,:,c),NPP_H(j,k,c),ANPP_H(j,k,c),Rg_H(j,k,c), ...
                     RA_H(j,k,c),Rms_H(j,k,c),Rmr_H(j,k,c),Rmc_H(j,k,c),PHE_S_H(j,k,c), ...
                     dflo_H(j,k,c),AgeL_H(j,k,c),e_rel_H(j,k,c),e_relN_H(j,k,c), ...
                     LAIdead_H(j,k,c),NBLeaf_H(j,k,c), ...
                     Sr_H(j,k,c),Slf_H(j,k,c),Sfr_H(j,k,c),Swm_H(j,k,c),Sll_H(j,k,c), ...
                     Rexmy_H(j,k,:,c),Rrootl_H(j,k,c), ...
                     AgeDL_H(j,k,c),Bfac_dayH(j,k,c),Bfac_weekH(j,k,c),NPPI_H(j,k,c), ...
                     TdpI_H(j,k,c),NupI_H(j,k,:,c),PARI_H(j,k,:,c),NBLI_H(j,k,c), ...
                     RB_H(j,k,:,c),FNC_H(j,k,c), ...
                     Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                     Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                     rNc_H(j,k,c),rPc_H(j,k,c),rKc_H(j,k,c), ...
                     rCac_H(j,k,c),rMgc_H(j,k,c),rSic_H(j,k,c),ManIH(k,c),SeedInput_cc] = ...
                     VEGGIE_UNIT( ...
                        B_H(j-1,k,:,c),PHE_S_H(j-1,k,c),dflo_H(j-1,k,c), ...
                        AgeL_H(j-1,k,c),AgeDL_H(j-1,k,c), ...
                        Ta(pdind),Tdp_H(pdind,k,c),PARB(pdind)+PARD(pdind), ...
                        Psi_x_H(pdind,k,c),Psi_l_H(pdind,k,c), ...
                        An_H(pdind,k,c),Rdark_H(pdind,k,c),NPP_H(j-1,k,c),jDay(j),Datam(i,:), ...
                        NPPI_H(j-1,k,c),TdpI_H(j-1,k,c),Bfac_weekH(j-1,k,c), ...
                        NupI_H(j-1,k,:,c),NavlI(j,:,c),PARI_H(j-1,k,:,c), ...
                        NBLI_H(j-1,k,c),NBLeaf_H(j-1,k,c), ...
                        L_day(j),p.Lmax_day,p.VegH_Param_Dyn,k, ...
                        Nreserve_H(j-1,k,c),Preserve_H(j-1,k,c),Kreserve_H(j-1,k,c), ...
                        Careserve_H(j-1,k,c),Mgreserve_H(j-1,k,c),Sireserve_H(j-1,k,c), ...
                        Nuptake_H(j,k,c),Puptake_H(j,k,c),Kuptake_H(j,k,c), ...
                        Cauptake_H(j,k,c),Mguptake_H(j,k,c),Siuptake_H(j,k,c), ...
                        FNC_H(j-1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_H(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_H(k), ...
                        TBio_Ht(j-1,k,c),OPT_EnvLimitGrowth_c(c),OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [TexC_H(j,k,c),TexN_H(j,k,c),TexP_H(j,k,c),TexK_H(j,k,c), ...
                     TexCa_H(j,k,c),TexMg_H(j,k,c),TexSi_H(j,k,c), ...
                     TNIT_H(j,k,c),TPHO_H(j,k,c),TPOT_H(j,k,c), ...
                     TCAL_H(j,k,c),TMAG_H(j,k,c),TSIL_H(j,k,c), ...
                     NuLit_H(j,k,:,c),Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                     Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                     SupN_H(j,k,c),SupP_H(j,k,c),SupK_H(j,k,c), ...
                     SupCa_H(j,k,c),SupMg_H(j,k,c),SupSi_H(j,k,c), ...
                     ISOIL_H(j,k,:,c)] = Plant_Exports(B_H(j,k,:,c),B_H(j-1,k,:,c), ...
                        NuLit_H(j-1,k,:,c),Slf_H(j,k,c),Sfr_H(j,k,c),Swm_H(j,k,c), ...
                        Sll_H(j,k,c),Sr_H(j,k,c), ...
                        Rexmy_H(j,k,:,c),p.Stoich_H(k),p.Mpar_H(k),p.fab_H(k),p.fbe_H(k), ...
                        RB_H(j,k,:,c),Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                        Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                        rNc_H(j,k,c),rPc_H(j,k,c),rKc_H(j,k,c), ...
                        rCac_H(j,k,c),rMgc_H(j,k,c),rSic_H(j,k,c), ...
                        ManIH(k,c),OPT_SoilBiogeochemistry,id_location);

                    SAI_H(j,k,c) = SAI_H(j-1,k,c);
                    if p.aSE_H(k) == 2
                        hc_H(j,k,c) = GrassHeight(LAI_H(j,k,c),LAIdead_H(j,k,c));
                    elseif p.aSE_H(k) == 5
                        [hc_H(j,k,c),SAI_H(j,k,c),B_H(j,:,:,c),Ccrown,ZR95_H, ...
                         Nreserve_H(j,:,c),Preserve_H(j,:,c),Kreserve_H(j,:,c), ...
                         Careserve_H(j,:,c),Mgreserve_H(j,:,c),Sireserve_H(j,:,c), ...
                         AgrHarNut_cc,RfH_Zs] = CropHeightType( ...
                            LAI_H(j,k,c),LAIdead_H(j,k,c),k,ZR95_H,B_H(j,:,:,c), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            Nreserve_H(j,:,c),Preserve_H(j,:,c),Kreserve_H(j,:,c), ...
                            Careserve_H(j,:,c),Mgreserve_H(j,:,c),Sireserve_H(j,:,c), ...
                            ManIH(:,c),p.Mpar_H,p.VegH_Param_Dyn,OPT_SoilBiogeochemistry);
                        AgrHarNut(j,:,c) = Plant_Harvest_Reserve_Export( ...
                            AgrHarNut(j,:,c),AgrHarNut_cc,Ccrown(k));
                        SeedInput(j,:,c) = Plant_Seed_Input_Accumulator( ...
                            SeedInput(j,:,c),SeedInput_cc,Ccrown(k));
                    else
                        hc_H(j,k,c) = hc_H(j-1,k,c);
                    end
                else
                    LAI_H(j,k,c)=0; B_H(j,k,:,c)=0; NPP_H(j,k,c)=0; ANPP_H(j,k,c)=0;
                    Rg_H(j,k,c)=0; RA_H(j,k,c)=0; Rms_H(j,k,c)=0; Rmr_H(j,k,c)=0;
                    Rmc_H(j,k,c)=0; PHE_S_H(j,k,c)=0; dflo_H(j,k,c)=0; AgeL_H(j,k,c)=0;
                    e_rel_H(j,k,c)=0; e_relN_H(j,k,c)=0; SAI_H(j,k,c)=0; hc_H(j,k,c)=0;
                    LAIdead_H(j,k,c)=LAIdead_H(j-1,k,c);
                    Sr_H(j,k,c)=0; Slf_H(j,k,c)=0; Sfr_H(j,k,c)=0; Swm_H(j,k,c)=0;
                    Sll_H(j,k,c)=0; Rrootl_H(j,k,c)=0; Rexmy_H(j,k,:,c)=0;
                    AgeDL_H(j,k,c)=0; FNC_H(j,k,c)=1;
                    Nreserve_H(j,k,c)=0; Preserve_H(j,k,c)=0; Kreserve_H(j,k,c)=0;
                    Careserve_H(j,k,c)=0; Mgreserve_H(j,k,c)=0; Sireserve_H(j,k,c)=0;
                    rNc_H(j,k,c)=1; rPc_H(j,k,c)=1; rKc_H(j,k,c)=1;
                    rCac_H(j,k,c)=1; rMgc_H(j,k,c)=1; rSic_H(j,k,c)=1;
                    Bfac_dayH(j,k,c)=0; Bfac_weekH(j,k,c)=0; NPPI_H(j,k,c)=0;
                    TdpI_H(j,k,c)=0; NupI_H(j,k,:,c)=0; RB_H(j,k,:,c)=0;
                    TexC_H(j,k,c)=0; TexN_H(j,k,c)=0; TexP_H(j,k,c)=0; TexK_H(j,k,c)=0;
                    TexCa_H(j,k,c)=0; TexMg_H(j,k,c)=0; TexSi_H(j,k,c)=0;
                    TNIT_H(j,k,c)=0; TPHO_H(j,k,c)=0; TPOT_H(j,k,c)=0;
                    TCAL_H(j,k,c)=0; TMAG_H(j,k,c)=0; TSIL_H(j,k,c)=0;
                    ISOIL_H(j,k,:,c)=0; NuLit_H(j,k,:,c)=0;
                end

                if ZR95_L(k) > 0
                    [LAI_L(j,k,c),B_L(j,k,:,c),NPP_L(j,k,c),ANPP_L(j,k,c),Rg_L(j,k,c), ...
                     RA_L(j,k,c),Rms_L(j,k,c),Rmr_L(j,k,c),Rmc_L(j,k,c),PHE_S_L(j,k,c), ...
                     dflo_L(j,k,c),AgeL_L(j,k,c),e_rel_L(j,k,c),e_relN_L(j,k,c), ...
                     LAIdead_L(j,k,c),NBLeaf_L(j,k,c), ...
                     Sr_L(j,k,c),Slf_L(j,k,c),Sfr_L(j,k,c),Swm_L(j,k,c),Sll_L(j,k,c), ...
                     Rexmy_L(j,k,:,c),Rrootl_L(j,k,c), ...
                     AgeDL_L(j,k,c),Bfac_dayL(j,k,c),Bfac_weekL(j,k,c),NPPI_L(j,k,c), ...
                     TdpI_L(j,k,c),NupI_L(j,k,:,c),PARI_L(j,k,:,c),NBLI_L(j,k,c), ...
                     RB_L(j,k,:,c),FNC_L(j,k,c), ...
                     Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                     Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                     rNc_L(j,k,c),rPc_L(j,k,c),rKc_L(j,k,c), ...
                     rCac_L(j,k,c),rMgc_L(j,k,c),rSic_L(j,k,c),ManIL(k,c),SeedInput_cc] = ...
                     VEGGIE_UNIT( ...
                        B_L(j-1,k,:,c),PHE_S_L(j-1,k,c),dflo_L(j-1,k,c), ...
                        AgeL_L(j-1,k,c),AgeDL_L(j-1,k,c), ...
                        Ta(pdind),Tdp_L(pdind,k,c),PARB(pdind)+PARD(pdind), ...
                        Psi_x_L(pdind,k,c),Psi_l_L(pdind,k,c), ...
                        An_L(pdind,k,c),Rdark_L(pdind,k,c),NPP_L(j-1,k,c),jDay(j),Datam(i,:), ...
                        NPPI_L(j-1,k,c),TdpI_L(j-1,k,c),Bfac_weekL(j-1,k,c), ...
                        NupI_L(j-1,k,:,c),NavlI(j,:,c),PARI_L(j-1,k,:,c), ...
                        NBLI_L(j-1,k,c),NBLeaf_L(j-1,k,c), ...
                        L_day(j),p.Lmax_day,p.VegL_Param_Dyn,k, ...
                        Nreserve_L(j-1,k,c),Preserve_L(j-1,k,c),Kreserve_L(j-1,k,c), ...
                        Careserve_L(j-1,k,c),Mgreserve_L(j-1,k,c),Sireserve_L(j-1,k,c), ...
                        Nuptake_L(j,k,c),Puptake_L(j,k,c),Kuptake_L(j,k,c), ...
                        Cauptake_L(j,k,c),Mguptake_L(j,k,c),Siuptake_L(j,k,c), ...
                        FNC_L(j-1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_L(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_L(k), ...
                        TBio_Lt(j-1,k,c),OPT_EnvLimitGrowth_c(c),OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [TexC_L(j,k,c),TexN_L(j,k,c),TexP_L(j,k,c),TexK_L(j,k,c), ...
                     TexCa_L(j,k,c),TexMg_L(j,k,c),TexSi_L(j,k,c), ...
                     TNIT_L(j,k,c),TPHO_L(j,k,c),TPOT_L(j,k,c), ...
                     TCAL_L(j,k,c),TMAG_L(j,k,c),TSIL_L(j,k,c), ...
                     NuLit_L(j,k,:,c),Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                     Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                     SupN_L(j,k,c),SupP_L(j,k,c),SupK_L(j,k,c), ...
                     SupCa_L(j,k,c),SupMg_L(j,k,c),SupSi_L(j,k,c), ...
                     ISOIL_L(j,k,:,c)] = Plant_Exports(B_L(j,k,:,c),B_L(j-1,k,:,c), ...
                        NuLit_L(j-1,k,:,c),Slf_L(j,k,c),Sfr_L(j,k,c),Swm_L(j,k,c), ...
                        Sll_L(j,k,c),Sr_L(j,k,c), ...
                        Rexmy_L(j,k,:,c),p.Stoich_L(k),p.Mpar_L(k),p.fab_L(k),p.fbe_L(k), ...
                        RB_L(j,k,:,c),Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                        Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                        rNc_L(j,k,c),rPc_L(j,k,c),rKc_L(j,k,c), ...
                        rCac_L(j,k,c),rMgc_L(j,k,c),rSic_L(j,k,c), ...
                        ManIL(k,c),OPT_SoilBiogeochemistry,id_location);

                    SAI_L(j,k,c) = SAI_L(j-1,k,c);
                    if p.aSE_L(k) == 2
                        hc_L(j,k,c) = GrassHeight(LAI_L(j,k,c),LAIdead_L(j,k,c));
                    elseif p.aSE_L(k) == 5
                        [hc_L(j,k,c),SAI_L(j,k,c),B_L(j,:,:,c),Ccrown,ZR95_L, ...
                         Nreserve_L(j,:,c),Preserve_L(j,:,c),Kreserve_L(j,:,c), ...
                         Careserve_L(j,:,c),Mgreserve_L(j,:,c),Sireserve_L(j,:,c), ...
                         AgrHarNut_cc,RfL_Zs] = CropHeightType( ...
                            LAI_L(j,k,c),LAIdead_L(j,k,c),k,ZR95_L,B_L(j,:,:,c), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            Nreserve_L(j,:,c),Preserve_L(j,:,c),Kreserve_L(j,:,c), ...
                            Careserve_L(j,:,c),Mgreserve_L(j,:,c),Sireserve_L(j,:,c), ...
                            ManIL(:,c),p.Mpar_L,p.VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                        AgrHarNut(j,:,c) = Plant_Harvest_Reserve_Export( ...
                            AgrHarNut(j,:,c),AgrHarNut_cc,Ccrown(k));
                        SeedInput(j,:,c) = Plant_Seed_Input_Accumulator( ...
                            SeedInput(j,:,c),SeedInput_cc,Ccrown(k));
                    else
                        hc_L(j,k,c) = hc_L(j-1,k,c);
                    end
                else
                    LAI_L(j,k,c)=0; B_L(j,k,:,c)=0; NPP_L(j,k,c)=0; ANPP_L(j,k,c)=0;
                    Rg_L(j,k,c)=0; RA_L(j,k,c)=0; Rms_L(j,k,c)=0; Rmr_L(j,k,c)=0;
                    PHE_S_L(j,k,c)=0; dflo_L(j,k,c)=0; AgeL_L(j,k,c)=0;
                    e_rel_L(j,k,c)=0; e_relN_L(j,k,c)=0; SAI_L(j,k,c)=0; hc_L(j,k,c)=0;
                    LAIdead_L(j,k,c)=LAIdead_L(j-1,k,c);
                    Sr_L(j,k,c)=0; Slf_L(j,k,c)=0; Sfr_L(j,k,c)=0; Swm_L(j,k,c)=0;
                    Sll_L(j,k,c)=0; Rrootl_L(j,k,c)=0; Rexmy_L(j,k,:,c)=0;
                    AgeDL_L(j,k,c)=0; FNC_L(j,k,c)=1;
                    Nreserve_L(j,k,c)=0; Preserve_L(j,k,c)=0; Kreserve_L(j,k,c)=0;
                    Careserve_L(j,k,c)=0; Mgreserve_L(j,k,c)=0; Sireserve_L(j,k,c)=0;
                    rNc_L(j,k,c)=1; rPc_L(j,k,c)=1; rKc_L(j,k,c)=1;
                    rCac_L(j,k,c)=1; rMgc_L(j,k,c)=1; rSic_L(j,k,c)=1;
                    Bfac_dayL(j,k,c)=0; Bfac_weekL(j,k,c)=0; NPPI_L(j,k,c)=0;
                    TdpI_L(j,k,c)=0; NupI_L(j,k,:,c)=0; RB_L(j,k,:,c)=0;
                    TexC_L(j,k,c)=0; TexN_L(j,k,c)=0; TexP_L(j,k,c)=0; TexK_L(j,k,c)=0;
                    TexCa_L(j,k,c)=0; TexMg_L(j,k,c)=0; TexSi_L(j,k,c)=0;
                    TNIT_L(j,k,c)=0; TPHO_L(j,k,c)=0; TPOT_L(j,k,c)=0;
                    TCAL_L(j,k,c)=0; TMAG_L(j,k,c)=0; TSIL_L(j,k,c)=0;
                    ISOIL_L(j,k,:,c)=0; NuLit_L(j,k,:,c)=0;
                end
            end

            %% Crops may resize Ccrown / rooting depth: persist per column
            Ccrown_c(c,:) = Ccrown;
            ZR95_H_c(c,:) = ZR95_H;  ZR95_L_c(c,:) = ZR95_L;
            RfH_Zs_c{c} = RfH_Zs;    RfL_Zs_c{c} = RfL_Zs;
            if OPT_VCA == 0
                Ccrown_t(j,:,c) = Ccrown;
                Cbare_c(c) = 1 - sum(Ccrown) - p.Cwat - p.Curb - p.Crock;
                Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,Cbare_c(c),Ccrown);
            end
        end
    end

    %% ------------------------------------------------------------------
    %% HOURLY BLOCK  --  vertical physics, one column at a time
    %% ------------------------------------------------------------------
    for c = 1:N_col
        p = PAR{c};
        Ccrown = Ccrown_c(c,:);
        if N_col > 1
            %% The FV graph (Darcy_Lateral_Topology, total head H = psi + z)
            %% owns all subsurface lateral flow, so the column runs strictly
            %% 1-D: a non-zero Slo_top here would double-count via the ODE's
            %% own kinematic Qi_out (SOIL_MOISTURES_RICH_COMP, driven by
            %% sin(atan(Slo_top))).  See the WP4 note near the slope guard.
            Qi_in(i,:,c) = 0;
            Slo_top_cell = 0;
            Slo_pot_cell = zeros(1,ms);
            q_runon(i,1,c) = q_runon_next(c);   %% routed runoff, previous step
        else
            Qi_in(i,:,c) = Qi_out(i-1,:,c);
            Slo_top_cell = p.Slo_top;
            Slo_pot_cell = p.Slo_pot;
            q_runon(i,1,c) = 0;
        end

        [V(i,:,c),Vice(i,:,c),O(i,:,c),Oice(i,:,c),ZWT(i,c),OF(i,c),OS(i,c), ...
         OH(i,:,c),OL(i,:,c),Psi_s_H(i,:,c),Psi_s_L(i,:,c),Rd(i,c),Qi_out(i,:,c),WTR(i,:,c), ...
         Rh(i,c),Lk(i,c),f(i,c),WIS(i,c),Ts(i,c),Pr_sno(i,c),Pr_liq(i,c),Csno(i,c), ...
         Cice(i,c),NDVI(i,c),rb_H(i,:,c),rb_L(i,:,c),rs_sunH(i,:,c), ...
         rs_sunL(i,:,c),rs_shdH(i,:,c),rs_shdL(i,:,c),r_litter(i,:,c), ...
         An_L(i,:,c),An_H(i,:,c),Rdark_L(i,:,c),Rdark_H(i,:,c), ...
         Ci_sunH(i,:,c),Ci_sunL(i,:,c),Ci_shdH(i,:,c),Ci_shdL(i,:,c), ...
         rap_H(i,:,c),rap_L(i,:,c),r_soil(i,c),b_soil(i,c),alp_soil(i,c),ra(i,c),Rn(i,c), ...
         H(i,c),QE(i,c),Qv(i,c),Lpho(i,c),T_H(i,:,c),T_L(i,:,c),EIn_H(i,:,c),EIn_L(i,:,c), ...
         EG(i,c),ESN(i,c),ESN_In(i,c),ELitter(i,c),EWAT(i,c),EICE(i,c),EIn_urb(i,c), ...
         EIn_rock(i,c),dw_SNO(i,c), ...
         G(i,c),Gfin(i,c),Tdp(i,:,c),Tdeb(i,:,c),Tdamp(i,c),Tice(i,c), ...
         Tdp_H(i,:,c),Tdp_L(i,:,c),SWE(i,c),SND(i,c),ros(i,c),In_SWE(i,c),SP_wc(i,c), ...
         WR_SP(i,c),U_SWE(i,c),NIn_SWE(i,c),dQ(i,c),Qfm(i,c),t_sls(i,c),DQ(i,c),DT(i,c), ...
         WAT(i,c),ICE(i,c),ICE_D(i,c),IP_wc(i,c),WR_IP(i,c),NIce(i,c),Cicew(i,c), ...
         Csnow(i,c),FROCK(i,c),Imelt(i,c),Smelt(i,c), ...
         In_H(i,:,c),In_L(i,:,c),In_Litter(i,c),In_urb(i,c),In_rock(i,c), ...
         Dr_H(i,:,c),Dr_L(i,:,c),SE_rock(i,c),SE_urb(i,c),Lk_wat(i,c),Lk_rock(i,c),er(i,c), ...
         gsr_H(i,:,c),Psi_x_H(i,:,c),Psi_l_H(i,:,c),Jsx_H(i,:,c),Jxl_H(i,:,c), ...
         Kleaf_H(i,:,c),Kx_H(i,:,c),Vx_H(i,:,c),Vl_H(i,:,c), ...
         gsr_L(i,:,c),Psi_x_L(i,:,c),Psi_l_L(i,:,c),Jsx_L(i,:,c),Jxl_L(i,:,c), ...
         Kleaf_L(i,:,c),Kx_L(i,:,c),Vx_L(i,:,c),Vl_L(i,:,c), ...
         fapar_H(i,:,c),fapar_L(i,:,c),SIF_H(i,:,c),SIF_L(i,:,c), ...
         snow_alb_new,tau_sno(i,c),e_sno(i,c),Ws_under(i,c),dQVEG(i,c),HV(i,c), ...
         QEV(i,c),TsV(i,c),EK(i,c),POT(i,:,c)] = ...
         HYDROLOGIC_UNIT(V(i-1,:,c),Oice(i-1,:,c),p.aR,p.Zs, ...
            p.EvL_Zs,p.Inf_Zs,p.Zinf,RfH_Zs_c{c},RfL_Zs_c{c},p.dz,p.Dz,ms,p.Kbot, ...
            Pr(i),Ta(i),Ds(i),Ws(i),p.zatm,Ts(i-1,c),IrD(i),dt,dth,ea(i),N(i),Pre(i),Tstm0_c(c), ...
            LAI_H(j,:,c),SAI_H(j,:,c),LAI_L(j,:,c),SAI_L(j,:,c),LAIdead_H(j,:,c),LAIdead_L(j,:,c), ...
            Rrootl_H(j,:,c),Rrootl_L(j,:,c),BLit(j,:,c),p.Sllit,p.Kct, ...
            Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft,h_S, ...
            Ccrown,Cbare_c(c),p.Crock,p.Curb,p.Cwat, ...
            p.Soil_Param,p.Interc_Param,p.SnowIce_Param,p.VegH_Param,p.VegL_Param, ...
            p.Zs_deb,p.Deb_Par, ...
            ZR95_H_c(c,:),ZR95_L_c(c,:), ...
            SAB1_c(c),SAB2_c(c),SAD1_c(c),SAD2_c(c),PARB_c(c),PARD_c(c),SvF_c(c), ...
            SND(i-1,c),snow_alb_c{c},p.Color_Class,p.OM_H,p.OM_L, ...
            p.PFT_opt_H,p.PFT_opt_L,hc_H(j,:,c),hc_L(j,:,c),p.d_leaf_H,p.d_leaf_L, ...
            Ca(i),Oa,Ci_sunH(i-1,:,c),Ci_shdH(i-1,:,c),Ci_sunL(i-1,:,c),Ci_shdL(i-1,:,c), ...
            e_rel_H(j,:,c),e_relN_H(j,:,c),e_rel_L(j,:,c),e_relN_L(j,:,c), ...
            e_sno(i-1,c),In_H(i-1,:,c),In_L(i-1,:,c),In_Litter(i-1,c),In_urb(i-1,c), ...
            In_rock(i-1,c),SWE(i-1,c),In_SWE(i-1,c), ...
            Tdeb(i-1,:,c),Tdp(i-1,:,c),Tdamp(i-1,c),Tice(i-1,c), ...
            WAT(i-1,c),ICE(i-1,c),IP_wc(i-1,c),ICE_D(i-1,c),Cicew(i-1,c), ...
            Vx_H(i-1,:,c),Vl_H(i-1,:,c),Vx_L(i-1,:,c),Vl_L(i-1,:,c), ...
            Psi_x_H(i-1,:,c),Psi_l_H(i-1,:,c),Psi_x_L(i-1,:,c),Psi_l_L(i-1,:,c), ...
            FROCK(i-1,c),p.Krock,Ws_under(i-1,c), ...
            Tdew(i),t_sls(i-1,c),ros(i-1,c),SP_wc(i-1,c),p.fpr,Pr_sno(pdind,1,c), ...
            p.In_max_urb,p.In_max_rock,p.K_usle,tau_sno(i-1,c),Ta(pdind), ...
            Slo_top_cell,Slo_pot_cell,p.Asur,p.Ared,p.aTop,EK(i-1,c),q_runon(i,1,c),Qi_in(i,:,c), ...
            p.pow_dis,p.a_dis, ...
            p.SPAR,p.SN,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST, ...
            OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);

        snow_alb_c{c} = snow_alb_new;
        Tstm0_c(c) = 2*Ts(i,c) - Ts(i-1,c);
        STOT = SAB1_c(c)+SAB2_c(c)+SAD1_c(c)+SAD2_c(c);
        ALB(i,c) = SAB1_c(c)/STOT*snow_alb_new.dir_vis + SAD1_c(c)/STOT*snow_alb_new.dif_vis + ...
                   SAB2_c(c)/STOT*snow_alb_new.dir_nir + SAD2_c(c)/STOT*snow_alb_new.dif_nir;

        %% Vertical exit at the base of the column [m3]
        exit_water_vert_m3(c) = exit_water_vert_m3(c) + (Lk(i,c)/1000)*A_cell(c)*dth;
        if ~topography_attached
            %% Flat graph: unrouted runoff leaves each cell (unchanged).
            exit_water_surface_m3(c) = exit_water_surface_m3(c) + ...
                ((Rh(i,c)+Rd(i,c))/1000)*A_cell(c);
        end
        %% With topography attached, this step's ponded runoff is collected
        %% and routed after the lateral block (Surface_Runoff_Routing).
    end

    if lateral_track_on
        lateral_track = Lateral_Track_Record('runon',lateral_track,i, ...
            reshape(q_runon(i,1,:),N_col,1));
    end

    %% ------------------------------------------------------------------
    %% LATERAL COUPLING  --  operator-split, after every column has
    %% taken its vertical step.
    %% ------------------------------------------------------------------
    if N_col > 1
        theta_cell = reshape(O(i,:,:), ms, N_col).';        %% [N_col x ms]
        Oice_cell = reshape(Oice(i,:,:), ms, N_col).';
        T_cell = reshape(Tdp(i,:,:), ms, N_col).';

        %% Biogeo_Leakage defines which T&C pools are dissolved and mobile.
        %% The pools are vertically integrated [gX/m2]; the T&C leakage
        %% formulation treats them as well mixed over the column water.
        P_cells = reshape(P(j,:,:),79,N_col).';
        [~,mobile] = Biogeo_Lateral_Pools(P_cells);
        mobile = mobile .* lateral_fraction;
        C_aq = Biogeo_Lateral_Concentration( ...
            mobile,A_cell,theta_cell,Zlayer_m,bio_active_fraction);

        lambda_cell = zeros(N_col,ms);
        capp_cell = zeros(N_col,ms);
        for c = 1:N_col
            p = PAR{c};
            [lambda_cell(c,:),capp_cell(c,:)] = Soil_Lateral_Thermal_State( ...
                T_cell(c,:),theta_cell(c,:)+Oice_cell(c,:), ...
                Pre(i),p,OPT_FR_SOIL);
        end
        thermal = struct('T_cell',T_cell,'lambda_cell',lambda_cell, ...
            'capp_cell',capp_cell,'Oice_cell',Oice_cell);
        thermal.parameter_cell = PAR;
        thermal.Pre = Pre(i);
        thermal.OPT_FR_SOIL = OPT_FR_SOIL;
        [theta_new, dM_cell, ldiag, T_new, Oice_new] = Lateral_Routing_Step( ...
            topo,theta_cell,C_aq,soilpar_cell,Zlayer_m,dth,thermal, ...
            bio_active_fraction);
        if lateral_track_on
            lateral_track = Lateral_Track_Record('subsurface',lateral_track,i,j, ...
                struct('topo',topo,'P_cells',P_cells,'dtheta_water',ldiag.dtheta_water, ...
                'dM_edge',ldiag.dM_edge,'Q_solute',ldiag.Q_solute_edge, ...
                'theta_before',theta_cell,'A_cell',A_cell,'Zlayer_m',Zlayer_m, ...
                'active_fraction',bio_active_fraction,'dth',dth));
        end

        %% Write the redistributed moisture back, keeping V consistent
        %% with O through the column's own residual water content.
        for c = 1:N_col
            p = PAR{c};
            O(i,:,c) = theta_new(c,:);
            Oice(i,:,c) = Oice_new(c,:);
            V(i,:,c) = (theta_new(c,:) - p.Ohy) .* p.dz;
            Vice(i,:,c) = Oice_new(c,:) .* p.dz;
            Tdp(i,:,c) = T_new(c,:);
            [Tdp_H(i,:,c),Tdp_L(i,:,c)] = RootZone_Temp( ...
                T_new(c,:),RfH_Zs_c{c},RfL_Zs_c{c});
        end

        %% Return the advected mass to its T&C pools [gX/m2].
        dM_area = squeeze(sum(dM_cell,2)) ./ A_cell;
        P_cells = Biogeo_Lateral_Pools(P_cells,dM_area);
        P(j,:,:) = reshape(P_cells.',1,79,N_col);
        transported = unique([lateral_pool_idx{:}]);
        Pj = P(j,transported,:);
        if any(Pj(:) < -1e-10)
            error('MAIN_FRAME_BG_ERW_2D:NegativePool', ...
                'Conservative lateral transport produced a negative pool at i=%d.',i);
        end

        %% --------------------------------------------------------------
        %% SURFACE RUNOFF ROUTING -- this step's ponded runoff is routed
        %% down the DEM, carrying its dissolved load.  It leaves the
        %% domain only at the topographic outlet(s); otherwise it is
        %% offered to the next vertical step as q_runon (re-infiltration).
        %% --------------------------------------------------------------
        if topography_attached
            Rh_c = reshape(Rh(i,1,:),N_col,1);
            Rd_c = reshape(Rd(i,1,:),N_col,1);
            runoff_step_mm = max(Rh_c + Rd_c, 0);                 % [N_col x 1]
            depth_m = runoff_step_mm/1000;                        % [N_col x 1]
            [~,mobile_now] = Biogeo_Lateral_Pools(P_cells);
            mobile_now = mobile_now .* lateral_fraction;          % [N_col x N_sp]
            C_soil = Biogeo_Lateral_Concentration( ...
                mobile_now,A_cell,theta_new,Zlayer_m,bio_active_fraction);
            C_soil = reshape(C_soil(:,1,:),N_col,N_lateral_species);

            %% TWO SOURCES OF SURFACE SOLUTE.
            %%  Saturation-excess runoff (Rd) is soil water reaching the
            %%  surface, so it carries the active-zone concentration.
            %%  Infiltration-excess runoff (Rh) is ponded water: rain,
            %%  drip, melt (~0 for these pools; wet deposition is handled
            %%  apart) mixed with whatever run-on arrived last step and has
            %%  not infiltrated -- the pond_mass_area store.  Crediting
            %%  that store to the soil on arrival would strand it in the
            %%  first absorbing cell it reached.
            %%  WIS is exactly the water Rh is defined against
            %%  (Rh = (WIS - f)*dth in HYDROLOGIC_UNIT): throughfall +
            %%  drip + snowmelt + run-on, in mm on the surface, so the share
            %%  of the pond that leaves again is Rh/WIS, capped at 1.
            WIS_c    = reshape(WIS(i,1,:),N_col,1);               % [mm]
            runon_mm = reshape(q_runon(i,1,:),N_col,1)*dth;       % [mm]
            out_frac = min(1, Rh_c./max(WIS_c,eps));
            out_frac(WIS_c <= 0) = 0;
            pond_out_area = pond_mass_area .* out_frac;           % leaves with Rh
            pond_in_area  = pond_mass_area - pond_out_area;       % infiltrates here
            pond_in_water_mm = runon_mm.*(1-out_frac);            % [mm]

            %% Soil-derived load, capped at the mobile pool it comes from.
            soil_out_area = min((Rd_c/1000).*C_soil, mobile_now);
            soil_out_area(~isfinite(soil_out_area) | soil_out_area < 0) = 0;
            %% Concentration of the runoff that actually leaves [gX/m3].
            C_surf = (soil_out_area + pond_out_area) ./ max(depth_m,eps);
            C_surf(~isfinite(C_surf) | C_surf < 0) = 0;
            %% Recompute both removals from the concentration that is
            %% really routed, so the ledger cannot drift from it.
            removed_area = depth_m .* C_surf;                     % [gX/m2]
            soil_share = soil_out_area ./ max(soil_out_area+pond_out_area,eps);
            soil_removed_area = removed_area .* soil_share;

            %% FIELD-ORIGIN provenance as passive tracer columns: routed
            %% water mixes at confluences, so provenance cannot be
            %% recovered afterwards; carried as a tracer it mixes linearly
            %% and comes out exact.  Water: the Rd part is this cell's own
            %% (its overlap fraction); the Rh part is surface water whose
            %% field share is the pond's own provenance plus local rain.
            %% Solute likewise.  Without this the field ledger would weight
            %% the export by the OUTLET cell's overlap -- a clipped
            %% perimeter cell -- and a stop-over would re-label the load.
            surf_field_mm = pond_field_water_mm + ...
                max(WIS_c - runon_mm,0).*field_fraction;
            rh_field_frac = surf_field_mm ./ max(WIS_c,eps);
            water_tracer = (Rd_c.*field_fraction + Rh_c.*rh_field_frac) ./ ...
                max(runoff_step_mm,eps);
            pond_field_out_area = pond_field_mass_area .* out_frac;
            mass_tracer = (soil_removed_area.*field_fraction + pond_field_out_area) ./ ...
                max(depth_m,eps);
            water_tracer(~isfinite(water_tracer)) = 0;
            mass_tracer(~isfinite(mass_tracer)) = 0;
            C_route = [C_surf, water_tracer, mass_tracer];
            [q_deliver_mm,C_deliver,~,exit_tracer,~,exit_cell_m3,exit_cell_g] = ...
                Surface_Runoff_Routing(topo,runoff_step_mm,C_route,dth);
            delivered_area = (q_deliver_mm/1000) .* C_deliver(:,1:N_lateral_species);

            %% Only what infiltrated out of last step's run-on joins the
            %% soil pools; the routed load stays on the surface, with its
            %% provenance.
            if lateral_track_on
                lateral_track = Lateral_Track_Record('surface',lateral_track,j, ...
                    struct('P_cells',P_cells,'pond_in_area',pond_in_area, ...
                    'soil_removed_area',soil_removed_area,'Rd_mm',Rd_c, ...
                    'theta',theta_new,'Zlayer_m',Zlayer_m, ...
                    'active_fraction',bio_active_fraction));
            end
            P_cells = Biogeo_Lateral_Pools(P_cells, pond_in_area - soil_removed_area);
            P(j,:,:) = reshape(P_cells.',1,79,N_col);
            pond_mass_area       = delivered_area;
            pond_field_water_mm  = q_deliver_mm .* C_deliver(:,N_lateral_species+1);
            pond_field_mass_area = (q_deliver_mm/1000) .* ...
                C_deliver(:,N_lateral_species+2:end);

            %% ERW aqueous stores (Na, Al, residual anion An, DIC and hence
            %% alkalinity) ride the same water as the T&C cations in the
            %% subsurface step, so they must be exported overland too --
            %% otherwise water leaving the field carries Ca and Mg but no
            %% DIC and the exported solution is not charge balanced, while
            %% the CDR ledger misses the overland pathway entirely.  Same
            %% two-limb, pond-store treatment as the pools above.
            if OPT_ERW == 1 && all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
                [ERW,erw_surf] = ERW_Surface_Export(ERW,topo,runoff_step_mm,Rd_c, ...
                    out_frac,erw_pond_umol,erw_pond_field_umol,pond_in_water_mm, ...
                    field_fraction,dth);
                erw_pond_umol       = erw_surf.delivered_umol;
                erw_pond_field_umol = erw_surf.delivered_field_umol;
                exit_mass_surface_erw_umol = ...
                    exit_mass_surface_erw_umol + erw_surf.exit_umol;
                exit_mass_surface_erw_field_umol = ...
                    exit_mass_surface_erw_field_umol + erw_surf.exit_field_umol;
            end

            %% q_runon is a RATE [mm/h]; q_deliver_mm is a depth per step.
            q_runon_next          = q_deliver_mm/dth;
            exit_water_surface_m3 = exit_water_surface_m3 + exit_cell_m3;
            exit_mass_surface_g   = exit_mass_surface_g + ...
                sum(exit_cell_g(:,1:N_lateral_species),1).';
            exit_water_field_surface_m3 = exit_water_field_surface_m3 + ...
                exit_tracer(N_lateral_species+1);
            exit_mass_surface_field_g = exit_mass_surface_field_g + ...
                exit_tracer(N_lateral_species+2:end);
            lateral_diag.surface_water_out_m3(j) = ...
                lateral_diag.surface_water_out_m3(j) + sum(exit_cell_m3);
            lateral_diag.surface_mass_out_g(j,:) = ...
                lateral_diag.surface_mass_out_g(j,:) + ...
                sum(exit_cell_g(:,1:N_lateral_species),1);
            lateral_diag.surface_reinfilt_m3(j) = ...
                lateral_diag.surface_reinfilt_m3(j) + sum((q_deliver_mm/1000).*A_cell);

            if any(reshape(P(j,transported,:),1,[]) < -1e-9)
                error('MAIN_FRAME_BG_ERW_2D:NegativePool', ...
                    'Surface runoff transport produced a negative pool at i=%d.',i);
            end
        end

        lateral_diag.water_buf_to_main(j) = lateral_diag.water_buf_to_main(j) + ldiag.buf_to_main.water;
        lateral_diag.water_edge_abs_m3(j) = lateral_diag.water_edge_abs_m3(j) + sum(abs(ldiag.Q_edge(:)))*dth;
        Q_face = sum(ldiag.Q_edge,2);                                   % m3/h, a -> b
        Q_cell_net = accumarray(topo.edge_b(:),Q_face,[N_col 1]) - ...
            accumarray(topo.edge_a(:),Q_face,[N_col 1]);                % m3/h gained
        lateral_diag.water_net_m3(j) = lateral_diag.water_net_m3(j) + 0.5*sum(abs(Q_cell_net))*dth;
        lateral_diag.mass_buf_to_main(j,:) = lateral_diag.mass_buf_to_main(j,:) + ldiag.buf_to_main.net.';
        lateral_diag.NO3_buf_to_main(j)   = lateral_diag.NO3_buf_to_main(j) + ldiag.buf_to_main.net(i_lateral_NO3);
        lateral_diag.sub_steps_max(j)     = max(lateral_diag.sub_steps_max(j), ldiag.sub_steps);
        lateral_diag.heat_buf_to_main_J(j)= lateral_diag.heat_buf_to_main_J(j) + ldiag.buf_to_main.heat_J;
        lateral_diag.heat_residual_J(j)   = lateral_diag.heat_residual_J(j) + ldiag.heat.energy_residual_J;
        lateral_diag.heat_residual_relative(j) = max( ...
            lateral_diag.heat_residual_relative(j), ...
            ldiag.heat.energy_residual_relative);
        lateral_diag.heat_sub_steps_max(j)= max(lateral_diag.heat_sub_steps_max(j),ldiag.heat.sub_steps);
    end

    if N_col > 1 && OPT_ERW == 1 && ...
            all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
        [ERW,erw_ldiag] = ERW_Lateral_Transport( ...
            ERW,topo,ldiag.Q_edge,theta_cell,Zlayer_m,bio_active_fraction, ...
            dth);
        lateral_diag.erw_mass_buf_to_main_umol(j,:) = ...
            lateral_diag.erw_mass_buf_to_main_umol(j,:) + erw_ldiag.buf_to_main.aqueous_umol.';
        lateral_diag.erw_CO2_buf_to_main_umol(j) = ...
            lateral_diag.erw_CO2_buf_to_main_umol(j) + erw_ldiag.buf_to_main.CO2_umol;
        lateral_diag.erw_mass_residual_umol(j,:) = ...
            lateral_diag.erw_mass_residual_umol(j,:) + erw_ldiag.mass_residual_umol.';
        lateral_diag.erw_mass_abs_residual_umol(j,:) = ...
            lateral_diag.erw_mass_abs_residual_umol(j,:) + abs(erw_ldiag.mass_residual_umol).';
        erw_throughput_umol = sum(abs(erw_ldiag.aqueous_edge_umol),1);
        erw_throughput_umol(4) = erw_throughput_umol(4) + ...
            sum(abs(erw_ldiag.CO2_edge_umol));
        lateral_diag.erw_mass_throughput_umol(j,:) = ...
            lateral_diag.erw_mass_throughput_umol(j,:) + erw_throughput_umol;
    end

    %% ------------------------------------------------------------------
    %% ERW  --  once per day, after the current lateral transport
    %% ------------------------------------------------------------------
    if OPT_ERW == 1 && ((Datam(i,4) == 1) || i == NN)
        ERW_next = cell(N_col,1);
        parfor (c = 1:N_col,parallel_workers)
            ERW_PARAM_cell = ERW_PARAM_input;
            if isfield(ColCfg,'ERW_PARAM')
                ERW_PARAM_cell = ColCfg(c).ERW_PARAM;
            end
            ERW_next{c} = ERW_Daily_Step_2D(ERW{c}, ColCfg(c), PAR{c}, ...
                i, j, NN, NNd, Datam, Date, Ca, ...
                O(:,:,c), Lk(:,1,c), f(:,1,c), T_H(:,:,c), T_L(:,:,c), ...
                Tdp_bio_c(c), VSUM_c(c), VTSUM_c(c), ...
                R_microbe(:,1,c), R_litter(:,1,c), R_ew(:,1,c), R_litter_sur(:,1,c), ...
                Rmr_H(:,:,c), Rmr_L(:,:,c), Ccrown_c(c,:), ...
                P(:,:,c), PAR{c}.MOD_PARAM_input, ERW_PARAM_cell, ...
                conv_mol, conv_Al, N_col);
        end
        for c = 1:N_col
            ERW{c} = ERW_next{c};
            if j == 2 && isfield(ERW{c},'P_initial')
                P(j-1,:,c) = ERW{c}.P_initial;
            end
            %% Cold-start backfill of the pH records (see the seed comment
            %% above the PH_c/PH_bulk_c initialisation). Rows 1:j were
            %% written before this first ERW call from the bulk observation,
            %% so the pore record carried the bulk value. bg row 1 is the
            %% t = 0 ERW initial condition: pH(1) the pore-water pH that
            %% Invert_Bulk_pH recovered from the bulk observation, pH_bulk(1)
            %% the bulk pH Porewater_to_Bulk_pH returns for that state (the
            %% observation itself, to solver tolerance). The daily loop in
            %% ERW_Daily_Step_2D starts at row 2, so row 1 is never
            %% overwritten. Same pairing as SCEPTER: one pore-water H state
            %% solved at t = 0, soil pH derived from it. Records only;
            %% PH_c/PH_bulk_c are updated from the daily step below.
            if j == 2 && isfield(ERW{c},'bg') && isfield(ERW{c}.bg,'pH') && ...
                    isfield(ERW{c}.bg,'pH_bulk') && ...
                    isfinite(ERW{c}.bg.pH(1,1)) && isfinite(ERW{c}.bg.pH_bulk(1,1))
                PH_pore_d(1:j,c) = ERW{c}.bg.pH(1,1);
                PH_bulk_d(1:j,c) = ERW{c}.bg.pH_bulk(1,1);
                PH_d(1:j,c)      = ERW{c}.bg.pH_bulk(1,1);
            end
            if ERW{c}.has_update
                PH_c(c) = ERW{c}.PH;      %% pore-water pH for ERW chemistry
                if isfield(ERW{c},'PH_pore') && ~isempty(ERW{c}.PH_pore)
                    PH_c(c) = ERW{c}.PH_pore;
                end
                if isfield(ERW{c},'PH_bulk') && ~isempty(ERW{c}.PH_bulk)
                    PH_bulk_c(c) = ERW{c}.PH_bulk;
                end
                P(j,60,c) = ERW{c}.P60;  P(j,61,c) = ERW{c}.P61;
                P(j,68,c) = ERW{c}.P68;  P(j,69,c) = ERW{c}.P69;
                P(j,52,c) = ERW{c}.P52;  P(j,53,c) = ERW{c}.P53;
                P(j,76,c) = ERW{c}.P76;
                %% The full-pool path replaces P43; the dissolved-fraction
                %% path returns only new mineral weathering P additively.
                if isfield(ERW{c},'P43')
                    P(j,43,c) = ERW{c}.P43;
                end
            end
        end
        clear ERW_next
    end

    %% Rolling checkpoint after a completed hourly iteration. Chemistry
    %% arrays retain the last completed daily state when saved mid-day.
    if ~isempty(ERW_ckpt_file) && ...
            (toc - ERW_ckpt_last_toc) >= ERW_ckpt_minutes*60
        ERW_ckpt_tmp = [tempname(fileparts(ERW_ckpt_file)) '.mat'];
        try
            save(ERW_ckpt_tmp, '-regexp', '^(?!cleanup$).*', '-v7.3');
            movefile(ERW_ckpt_tmp, ERW_ckpt_file, 'f');
            fprintf(['  [checkpoint] iter %d / %d saved\n'], i, NN);
        catch ME_ckpt
            if isfile(ERW_ckpt_tmp); delete(ERW_ckpt_tmp); end
            warning('ERW:CheckpointFailed', ...
                'Checkpoint at iter %d failed: %s', i, ME_ckpt.message);
        end
        ERW_ckpt_last_toc = toc;
    end
end

%% The particle-class histories were kept in daily blocks during the run
%% (ERW_Daily_Step_2D); the results carry the rows each column kept:
%% those set at initialisation and the last row of every ERW day.
if OPT_ERW == 1
    for c = 1:N_col
        ERW{c} = ERW_Particle_History('assemble',ERW{c});
    end
end

Computational_Time = toc;
fprintf('\nCOMPUTATIONAL TIME [h] %.3f   (%.1f ms per column-step)\n', ...
    Computational_Time/3600, 1000*Computational_Time/(NN*N_col));

%% ========================================================================
%% DOMAIN AGGREGATION  --  polygon-overlap weighted where available
%% ========================================================================
%% A_field_total is the whole field, including any field land a majority-
%% rule grid leaves in its dropped squares (carried at the field mean).
[A_field,w_col,~,A_field_total] = Field_Aggregation_Areas(topo);
A_simulated_total = sum(A_cell);
Domain.area_m2 = A_field_total;
Domain.simulated_area_m2 = A_simulated_total;
Domain.grid_outside_area_m2 = A_simulated_total-sum(A_field);
if isfield(topo,'field_missed_area_m2')
    Domain.field_missed_area_m2 = A_field_total-sum(A_field);
end
Domain.area_main_m2 = sum(A_field(strcmp(topo.cell_region,'main')));
Domain.area_buffer_m2 = sum(A_field(~strcmp(topo.cell_region,'main')));
[Domain.N2Oflx,Domain.N2O_total_gN] = Field_Aggregate_Areal(squeeze(N2Oflx),topo);
[Domain.N2flx,Domain.N2_total_gN] = Field_Aggregate_Areal(squeeze(N2flx),topo);
[Domain.VOL,Domain.VOL_total_gN] = Field_Aggregate_Areal(squeeze(VOL),topo);
[Domain.BfixN,Domain.BfixN_total_gN] = Field_Aggregate_Areal(squeeze(BfixN),topo);
[Domain.FertN,Domain.FertN_total_gN] = Field_Aggregate_Areal(squeeze(FertN_d),topo);
[Domain.DepN,Domain.DepN_total_gN] = Field_Aggregate_Areal(squeeze(DepN_d),topo);
[Domain.NITflx,Domain.NIT_total_gN] = Field_Aggregate_Areal(squeeze(NITflx),topo);
[Domain.DENITflx,Domain.DENIT_total_gN] = Field_Aggregate_Areal(squeeze(DENITflx),topo);
[Domain.N2Oflx_nit,Domain.N2O_nit_total_gN] = ...
    Field_Aggregate_Areal(squeeze(N2Oflx_nit),topo);
[Domain.LEAK_NO3,Domain.LEAK_NO3_total_gN] = ...
    Field_Aggregate_Areal(squeeze(LEAK_NO3),topo);
[Domain.LEAK_NH4,Domain.LEAK_NH4_total_gN] = ...
    Field_Aggregate_Areal(squeeze(LEAK_NH4),topo);
Domain.PH_bulk_cell_d = squeeze(PH_bulk_d);
Domain.PH_pore_cell_d = squeeze(PH_pore_d);
Domain.PH_bulk_d      = pH_Area_Weighted(Domain.PH_bulk_cell_d,w_col);
Domain.PH_pore_d      = pH_Area_Weighted(Domain.PH_pore_cell_d,w_col);
Domain.PH_bulk_aggregation = 'area_weighted_hydrogen_activity_index_not_pooled_extract';
Domain.PH_d           = Domain.PH_bulk_d;
[Domain.NO3_pool,Domain.NO3_pool_total_gN] = ...
    Field_Aggregate_Areal(squeeze(P(:,32,:)),topo);
Domain.lateral       = lateral_diag;
Domain.exit_water_vert_m3 = exit_water_vert_m3;
Domain.exit_water_surface_m3 = exit_water_surface_m3;
Domain.exit_water_total_m3 = exit_water_vert_m3+exit_water_surface_m3;
Domain.exit_water_field_vert_m3 = sum(exit_water_vert_m3.*field_fraction);
if topography_attached
    %% Routed export is attributed to the OUTLET cell, which on a
    %% polygon-clipped grid is a partial perimeter cell, so weighting it
    %% by that cell's own overlap fraction would discard most of the
    %% export.  The field-origin volume is carried through the router as a
    %% passive tracer instead, which mixes correctly at confluences.
    Domain.exit_water_field_surface_m3 = exit_water_field_surface_m3;
    Domain.exit_mass_surface_field_g = exit_mass_surface_field_g;
else
    Domain.exit_water_field_surface_m3 = sum(exit_water_surface_m3.*field_fraction);
end
Domain.exit_water_field_total_m3 = Domain.exit_water_field_vert_m3+ ...
    Domain.exit_water_field_surface_m3;
%% Solute exported at the topographic outlet(s) with overland runoff
%% [g element], one entry per lateral species (Domain.lateral.species_names).
Domain.exit_mass_surface_g = exit_mass_surface_g;
Domain.exit_mass_surface_species = lateral_species;
%% ERW aqueous species exported overland [umol], ordered as
%% Domain.lateral.erw_species_names (Na, Al, An, DIC).
Domain.exit_mass_surface_erw_umol = exit_mass_surface_erw_umol;
Domain.exit_mass_surface_erw_field_umol = exit_mass_surface_erw_field_umol;
Domain.cell_region   = topo.cell_region;
Domain.cell_area_m2  = A_cell;
Domain.cell_overlap_area_m2 = A_field;
Domain.cell_overlap_fraction = field_fraction;
Domain.field_aggregation = 'polygon_overlap_area';

fprintf('Domain: %.0f m2 field inside %.0f m2 simulated grid\n', ...
    Domain.area_m2,Domain.simulated_area_m2);
fprintf('Cumulative N2O   : %.4f gN/m2\n', sum(Domain.N2Oflx,'omitnan'));
fprintf('Cumulative NO3   : %.4f gN/m2\n', sum(Domain.LEAK_NO3,'omitnan'));
fprintf('Net lateral water margin->field: %.4f m3\n', ...
    sum(lateral_diag.water_buf_to_main,'omitnan'));
fprintf('Net lateral NO3   margin->field: %.4f g\n', ...
    sum(lateral_diag.NO3_buf_to_main,'omitnan'));

%% ========================================================================
function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end

function ctx = Cell_BG_State_Context(ctx,state)
%% The pool state v0.9's parameter scripts read, plus the exchanger seed
%% (Load_Cell_Parameters -> bulk_pH_params.initial_f_CEC / initial_Al_w).
%% Respiration, the Na-seed water loss and the plant reserves start the way
%% v0.9 starts them, from the run itself and the parameter file; v1.0: a
%% perennial stand's end state (Veg_out) goes on (Perennial_Vegetation_Apply).
names = {'B_out','O_out','RexmyI_out','NavlI_out','NupI_L_out', ...
    'spinup_erw_handoff_diagnostic','Veg_out'};
for k = 1:numel(names)
    if isfield(state,names{k})
        ctx.(names{k}) = state.(names{k});
    end
end
end
