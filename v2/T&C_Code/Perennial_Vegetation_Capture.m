%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perennial_Vegetation_Capture  (script)
%%
%% Run after the 10-minute driver: the stand at the end of the run, row j of
%% every state Perennial_Vegetation_States lists, as Veg_out, which a
%% perennial stand's vegetation file applies at the start of the next run
%% (Perennial_Vegetation_Apply).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
veg_names = Perennial_Vegetation_States();
Veg_out = struct();
for veg_k = 1:numel(veg_names)
    veg_x = eval(veg_names{veg_k});
    Veg_out.(veg_names{veg_k}) = veg_x(j,:,:);
end
clear veg_names veg_k veg_x
