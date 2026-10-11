function q = Biogeo_Layer_Drainage(layer,V_prev,Vice_prev,V_now,Vice_now, ...
    f,Rd,EG,Jsx_H,Jsx_L,Qi_in,Qi_out,EvL_Zs,RfH_Zs,RfL_Zs,Asur,Ared,dth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One hydrology step's drainage out of the bottom of the biogeochemical
%% layer [mm/h, positive downward]: the water the layer's dissolved pools
%% leach with under leaching_water 'layer' (Leaching_Water_Basis).
%%
%% T&C returns the leakage out of the bottom of the whole column (Lk) but
%% not the fluxes between its layers, so this is the step's water balance
%% of the layers the biogeochemistry spans (layer: Bio_Zs, nonzero where
%% it counts): what infiltrated less what left those layers by every other
%% path and what they kept,
%%
%%   q dth = f dth Asur Ared - Rd - (EG sum_top(EvL_Zs) + J_top) dth
%%           + sum_top(Qi_in - Qi_out) dth - Asur Ared dS_top
%%   J_top = sum_crowns(Jsx_H RfH_top + Jsx_L RfL_top),  S = sum_top(V + Vice)
%%
%% term by term the part of T&C's own step balance (CK1 in the drivers)
%% that falls in those layers: the step's whole balance is this plus the
%% same sum over the layers below, with Lk leaving their bottom.  Root
%% uptake and soil evaporation are split over the layers as the solver
%% splits them (SVAT_UNIT: Jsx*RfH_Zs with each crown's profile rescaled
%% to sum to its Jsx, EG*EvL_Zs).  A negative value is
%% water rising into the layer from below.  Exact when the biogeochemical
%% depth ends on a layer boundary (the default mesh and Zbio = 300 mm do);
%% otherwise it is the drainage out of the deepest layer the depth reaches.
%%
%% V_prev and Vice_prev are the state the step started from; V_now and
%% Vice_now the state the vertical step returned, before any lateral
%% exchange (the 2-D routing, a lateral replay) changes it, since lateral
%% transfer is not drainage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in = double(layer(:).' > 0);
dS = sum((V_now - V_prev + Vice_now - Vice_prev).*in);
J_top = root_uptake_in(Jsx_H,RfH_Zs,in) + root_uptake_in(Jsx_L,RfL_Zs,in);
q = (f*dth*Asur*Ared - Rd - (EG*sum(EvL_Zs.*in) + J_top)*dth ...
    + sum((Qi_in - Qi_out).*in)*dth - Asur*Ared*dS)/dth;
end

function J = root_uptake_in(Jsx,Rf,in)
%% The crowns' root uptake from the layers in 'in', distributed as
%% SVAT_UNIT distributes it: Jsx*Rf, each crown's row rescaled to sum to its
%% Jsx where that sum is positive.
share = Rf*in.';
total = sum(Rf,2);
Jsx = Jsx(:);
rescale = Jsx.*total > 0;
share(rescale) = share(rescale)./total(rescale);
J = sum(Jsx.*share);
end
