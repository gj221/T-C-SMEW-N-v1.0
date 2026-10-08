function chem = Representative_Soil_Chemistry(lab_file,target_pH,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's representative exchanger from laboratory-characterised pedons
%% (NRCS Kellogg Soil Survey Laboratory layers: exchangeable Ca, Mg, K and
%% Na by NH4OAc, KCl-extractable Al, CEC at pH 7, pH in water), as the
%% model takes it: the CEC at pH 7 and each cation as a fraction of it.
%%
%% The model's Appelo exchanger has a fixed capacity and carries the soil's
%% pH-dependent charge as exchangeable H (ERW_Resolve_CEC_pH_Sensitivity),
%% so its capacity is a buffered CEC, as the Energy Farm's is (13 cmolc/kg
%% by summation of Mehlich-3 bases and Sikora-buffer acidity, Pihlap et al.
%% 2025), not the effective CEC (bases + KCl Al) at the soil's own pH.  The
%% sites the bases and the KCl Al leave are the exchanger's H.
%%
%% The exchanger's composition follows its pH -- exchangeable Al replaces
%% the bases steeply below pH(H2O) ~5.5 -- so only layers of the active
%% depth whose pH lies within 'pHWindow' of the site's pH (target_pH, the
%% area mean of its cells) are pooled, and only those with a measured KCl Al
%% (an effective CEC).  Pooled, not averaged: CEC is the mean of those
%% layers' CEC at pH 7, each fraction is the layers' summed cation over
%% their summed CEC, so the fractions are of the very CEC reported, and
%% layers_pH, their CEC-weighted pH, is the pH that exchanger was measured
%% at.  ECEC_cmolc_kg, the layers' mean effective CEC, is reported beside.
%%
%% Name/value
%%   'DepthCm'    active depth [top bottom] (default [0 30]); a layer counts
%%                when it starts above the bottom
%%   'pHWindow'   half-width of the pH window (default 0.3)
%%   'SoilGrids'  SoilGrids table at the site (property, top_cm, bottom_cm,
%%                Q50 ...): its CEC at pH 7, thickness-weighted over the
%%                active depth, is reported as the cross-check
%%   'Output'     CSV to write the one-row result to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('DepthCm',[0 30],'pHWindow',0.3,'SoilGrids','','Output','');
for k = 1:2:numel(varargin)
    if ~isfield(opt,char(varargin{k}))
        error('Representative_Soil_Chemistry:BadOptions','Unknown option %s.',char(varargin{k}));
    end
    opt.(char(varargin{k})) = varargin{k+1};
end
L = readtable(lab_file,'TextType','string');
use = L.hzn_top < opt.DepthCm(2) & abs(L.ph_h2o - target_pH) <= opt.pHWindow & ...
    isfinite(L.ecec_base_plus_aluminum) & L.ecec_base_plus_aluminum > 0 & ...
    isfinite(L.aluminum_kcl_extractable) & isfinite(L.cec_nh4_ph_7) & L.cec_nh4_ph_7 > 0 & ...
    isfinite(L.ca_nh4_ph_7) & isfinite(L.mg_nh4_ph_7) & isfinite(L.k_nh4_ph_7) & ...
    isfinite(L.na_nh4_ph_7);
if ~any(use)
    error('Representative_Soil_Chemistry:NoLayers', ...
        'No layer of %s within %g of pH %.2f carries a CEC at pH 7 and a KCl Al.', ...
        lab_file,opt.pHWindow,target_pH);
end
P = L(use,:);
cec = P.cec_nh4_ph_7;
f = [sum(P.ca_nh4_ph_7) sum(P.mg_nh4_ph_7) sum(P.k_nh4_ph_7) sum(P.na_nh4_ph_7) ...
    sum(P.aluminum_kcl_extractable)]/sum(cec);
if any(f(1:3) <= 0) || sum(f) >= 1
    error('Representative_Soil_Chemistry:BadFractions', ...
        'Pooled fractions %s (Ca Mg K Na Al) cannot seed an exchanger.',mat2str(f,4));
end
layers = strjoin(compose('%s %g-%g cm pH %.1f',P.pedlabsampnum,P.hzn_top,P.hzn_bot,P.ph_h2o),'; ');
chem = table(mean(cec),f(1),f(2),f(3),f(4),f(5),sum(cec.*P.ph_h2o)/sum(cec),target_pH, ...
    height(P),string(layers),mean(P.ecec_base_plus_aluminum),'VariableNames', ...
    {'CEC_cmolc_kg','bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction','bs_Na_fraction', ...
    'Al_fraction','layers_pH','target_pH','n_layers','layers','ECEC_cmolc_kg'});
chem.soilgrids_cec7_cmolc_kg = NaN;
if ~isempty(opt.SoilGrids)
    G = readtable(opt.SoilGrids,'TextType','string');
    G = G(G.property == "cec" & G.top_cm < opt.DepthCm(2),:);
    w = min(G.bottom_cm,opt.DepthCm(2)) - max(G.top_cm,opt.DepthCm(1));
    chem.soilgrids_cec7_cmolc_kg = sum(w.*G.Q50)/sum(w);
end
if ~isempty(opt.Output)
    writetable(chem,opt.Output);
end
end
