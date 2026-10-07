function txt = Perennial_Stand_Summary(LAI_H,LAI_L,Veg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perennial_Stand_Summary
%%
%% One log line on a perennial stand after a year: the year's peak leaf
%% area of the trees and the grass in each crown, and the carbon the next
%% year starts from [gC m-2] -- the trees' fine roots, sapwood and reserve,
%% the grass's roots and reserve.  Leaf mass cannot exceed LtR times the
%% fine roots, and the spring flush draws on the reserve
%% (VEGETATION_DYNAMIC), so these show whether a canopy can grow back.
%%
%% LAI_H, LAI_L  the year's daily leaf area (days x crowns)
%% Veg           the stand handed on (Perennial_Vegetation_Capture)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txt = sprintf(['    perennial stand: peak LAI trees %s, grass %s; trees: fine roots %s, ' ...
    'sapwood %s, reserve %s; grass: roots %s, reserve %s gC m-2 at the year''s end'], ...
    mat2str(max(LAI_H,[],1),3),mat2str(max(LAI_L,[],1),3), ...
    mat2str(round(Veg.B_H(1,:,3))),mat2str(round(Veg.B_H(1,:,2))),mat2str(round(Veg.B_H(1,:,4))), ...
    mat2str(round(Veg.B_L(1,:,3))),mat2str(round(Veg.B_L(1,:,4))));
end
