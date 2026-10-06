function chem = Representative_Soil_Chemistry(lab_file,target_pH,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's representative exchanger from laboratory-characterised pedons
%% (NRCS Kellogg Soil Survey Laboratory layers: exchangeable Ca, Mg and K by
%% NH4OAc, KCl-extractable Al, CEC at pH 7, pH in water), as the model takes
%% it: the CEC at the soil's own pH (effective CEC, bases + exchangeable Al)
%% and Ca, Mg and K as fractions of it.
%%
%% The exchanger's composition follows its pH -- exchangeable Al replaces
%% the bases steeply below pH(H2O) ~5.5 -- so only layers of the active
%% depth whose pH lies within 'pHWindow' of the site's pH (target_pH, the
%% area mean of its cells) are pooled, and only those with a measured
%% effective CEC.  Pooled, not averaged: CEC is the mean of those layers and
%% each fraction is the layers' summed base over their summed effective CEC,
%% so the fractions are of the very CEC reported.
%%
%% Name/value
%%   'DepthCm'    active depth [top bottom] (default [0 30]); a layer counts
%%                when it starts above the bottom
%%   'pHWindow'   half-width of the pH window (default 0.3)
%%   'SoilGrids'  SoilGrids table at the site (property, top_cm, bottom_cm,
%%                Q50 ...): its CEC at pH 7, thickness-weighted over the
%%                active depth and scaled by the pooled layers' effective /
%%                pH-7 CEC ratio, is reported as a cross-check
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
    isfinite(L.ca_nh4_ph_7) & isfinite(L.mg_nh4_ph_7) & isfinite(L.k_nh4_ph_7);
if ~any(use)
    error('Representative_Soil_Chemistry:NoLayers', ...
        'No layer of %s within %g of pH %.2f carries an effective CEC.',lab_file,opt.pHWindow,target_pH);
end
P = L(use,:);
ecec = P.ecec_base_plus_aluminum;
f = [sum(P.ca_nh4_ph_7) sum(P.mg_nh4_ph_7) sum(P.k_nh4_ph_7)]/sum(ecec);
if any(f <= 0) || sum(f) >= 1
    error('Representative_Soil_Chemistry:BadFractions', ...
        'Pooled base fractions %s cannot seed an exchanger.',mat2str(f,4));
end
ratio = sum(ecec)/sum(P.cec_nh4_ph_7);
layers = strjoin(compose('%s %g-%g cm pH %.1f',P.pedlabsampnum,P.hzn_top,P.hzn_bot,P.ph_h2o),'; ');
chem = table(mean(ecec),f(1),f(2),f(3),target_pH,height(P),string(layers),ratio, ...
    'VariableNames',{'CEC_cmolc_kg','bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction', ...
    'target_pH','n_layers','layers','ecec_to_cec7'});
chem.soilgrids_cec7_cmolc_kg = NaN; chem.soilgrids_cec_eff_cmolc_kg = NaN;
if ~isempty(opt.SoilGrids)
    G = readtable(opt.SoilGrids,'TextType','string');
    G = G(G.property == "cec" & G.top_cm < opt.DepthCm(2),:);
    w = min(G.bottom_cm,opt.DepthCm(2)) - max(G.top_cm,opt.DepthCm(1));
    chem.soilgrids_cec7_cmolc_kg = sum(w.*G.Q50)/sum(w);
    chem.soilgrids_cec_eff_cmolc_kg = ratio*chem.soilgrids_cec7_cmolc_kg;
end
if ~isempty(opt.Output)
    writetable(chem,opt.Output);
end
end
