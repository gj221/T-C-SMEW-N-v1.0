function q = ERW_Mean_Water_Loss(T_H, T_L, Lk, Bio_Zs, RfH_Zs, RfL_Zs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The mean water loss of the ERW's biogeochemical layer over a run of the
%% driver [L m-2 d-1]: the transpiration its roots draw from the layer plus
%% the leakage, built as the drivers build the ERW's Tw and Lw
%% (MAIN_FRAME_BG_ERW_10mins, ERW_Daily_Step_2D), from the hourly T_H, T_L
%% [mm/h, one column per crown] and Lk [mm/h].  Row 1, the initial
%% condition, is left out.  The root profiles are the run's last: exact
%% for a perennial stand, while a crop's seasonal change in the share of
%% its roots in the layer, which the drivers' daily ERW step follows, is
%% not.
%%
%% The spin-up hands its ring's mean on as the long-term water loss the
%% rain Na seed divides by (seed_water_loss_l_m2_d, ERW_Seed_Water_Loss).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_layer = double(Bio_Zs(:)' > 0);
cfH = (RfH_Zs*in_layer')';
cfL = (RfL_Zs*in_layer')';
w = sum(cfH.*T_H(2:end,:) + cfL.*T_L(2:end,:), 2) + Lk(2:end);
q = 24*mean(w);
end
