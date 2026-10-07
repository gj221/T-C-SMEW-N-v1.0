function [names,has_depth] = Perennial_Vegetation_States()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perennial_Vegetation_States
%%
%% The daily vegetation states a perennial stand carries from one run of a
%% driver to the next (Perennial_Vegetation_Capture, Perennial_Vegetation_
%% Apply): every vegetation state the 10-minute and 2-D drivers read from
%% the previous day -- the (j-1) arguments of VEGGIE_UNIT and Plant_Exports,
%% and the litter and exudates BIOGEO_UNIT receives the day after they are
%% shed -- with the canopy the first day's hours see, for the trees (_H) and
%% the grass (_L).
%%
%% names      the variables, each held per day and crown (rows x cc)
%% has_depth  true where a variable also has a third dimension (pools,
%%            nutrients, litter fractions): rows x cc x k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer = {'B','LAI','LAIdead','hc','SAI','PHE_S','dflo','AgeL','AgeDL', ...
    'e_rel','e_relN','Rrootl','NPP','NPPI','TdpI','NupI','PARI','NBLI','NBLeaf','FNC', ...
    'Nreserve','Preserve','Kreserve','Careserve','Mgreserve','Sireserve','NuLit', ...
    'SupN','SupP','SupK','SupCa','SupMg','SupSi','ISOIL','Rexmy'};
depth = ismember(layer,{'B','NupI','PARI','NuLit','ISOIL','Rexmy'});
names = [strcat(layer,'_H') strcat(layer,'_L') {'Bfac_weekH','Bfac_weekL','TBio_Ht','TBio_Lt'}];
has_depth = [depth depth false(1,4)];
end
