function [O_i,V_i,rep] = Lateral_Replay_Water(rep,i,O_i,Oice_i,Ohy,Osat,dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral_Replay_Water
%%
%% Apply step i's recorded lateral (Darcy) water exchange to a replayed
%% column, after its vertical step, as the 2-D driver applied it: the
%% layer water contents move by the recorded change within the solver's
%% own bounds [Ohy, Osat - Oice], and V is rebuilt from them exactly as the
%% 2-D driver writes the lateral result back.  A change the bounds refuse
%% (the replayed column is drier or wetter than the field cell was) is not
%% applied and is booked in rep.water_clip_mm, so the replay's water
%% ledger stays exact.
%%
%%   O_i, Oice_i, Ohy, Osat, dz   [1 x ms]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_try = O_i + rep.dtheta(i,:);
O_new = min(max(O_try,Ohy),Osat - Oice_i);
rep.water_mm(i) = sum((O_new - O_i).*dz);
rep.water_clip_mm(i) = sum((O_try - O_new).*dz);
O_i = O_new;
V_i = (O_i - Ohy).*dz;
end
