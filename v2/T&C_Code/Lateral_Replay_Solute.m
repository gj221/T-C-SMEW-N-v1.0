function [Pj,rep] = Lateral_Replay_Solute(rep,j,Pj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral_Replay_Solute
%%
%% Apply day j's recorded lateral exchange of the mobile pools to a
%% replayed column, right after its daily biogeochemistry step:
%%     T <- retain .* T + input
%% (Lateral_Track_Record).  In the 2-D run the same exchange acted on the
%% day's pools hour by hour between this step and the next one, and only
%% the next step reads them, so applying the composed map at once gives
%% the next step the same pools.  The cell's own export scales with its own
%% pool; what arrived from neighbours is imposed as recorded.
%%
%% rep.solute_in_g / solute_out_g (g m-2) book the gross exchange for the
%% spin-up's nutrient budgets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = rep.pool_idx;
T = Pj(idx);
a = rep.retain(j,:);
b = rep.input(j,:);
rep.solute_out_g(j,:) = (1 - a).*T;
rep.solute_in_g(j,:) = b;
Pj(idx) = a.*T + b;
end
