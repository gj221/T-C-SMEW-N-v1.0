function T = Plant_Tissue_Nutrients(B,St,NuLit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plant_Tissue_Nutrients
%%
%% The nutrients held in one crown's tissues [g m-2 PFT], [N P K Ca Mg Si],
%% as Plant_Exports totals them each day (TNIT ... TSIL): each live pool's
%% carbon over its carbon:nutrient ratio -- leaves, sapwood, fine roots,
%% carbohydrate reserve (at the sapwood ratio), fruit, heartwood -- plus the
%% nutrients in leaf litter still in transit (NuLit).  The frames fill
%% TNIT ... TSIL from their second day, so this gives the stock of the
%% initial state the parameter file sets.
%%
%% B      carbon pools [gC m-2], at least the first six of the frame's eight
%% St     stoichiometry (Veg_Stoichiometric_Parameter)
%% NuLit  nutrients in transit [g m-2], [N P K Ca Mg Si]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = double(B(:).');
NuLit = double(NuLit(:).');
if numel(B) < 6 || numel(NuLit) ~= 6
    error('Plant_Tissue_Nutrients:BadInput', ...
        'B needs the six live carbon pools and NuLit six nutrients.');
end
ratio = [St.Nl    St.Ns    St.Nr    St.Ns    St.Nf    St.Nh; ...
         St.Phol  St.Phos  St.Phor  St.Phos  St.Phof  St.Phoh; ...
         St.Kpotl St.Kpots St.Kpotr St.Kpots St.Kpotf St.Kpoth; ...
         St.Call  St.Cals  St.Calr  St.Cals  St.Calf  St.Calh; ...
         St.Magl  St.Mags  St.Magr  St.Mags  St.Magf  St.Magh; ...
         St.Sill  St.Sils  St.Silr  St.Sils  St.Silf  St.Silh];
T = sum(B(1:6)./ratio,2).'+NuLit;
end
