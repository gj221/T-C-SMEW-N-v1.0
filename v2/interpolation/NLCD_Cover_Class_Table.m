function T = NLCD_Cover_Class_Table()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NLCD_Cover_Class_Table
%%
%% How each USGS NLCD land-cover class divides among the model's cover
%% components (Cell_Land_Cover):
%%   cover_tree   crown 1: tree canopy, with the low vegetation beneath it
%%   cover_grass  crown 2: open low vegetation
%%   cover_urban  impervious surface (Curb)
%%   cover_water  open water (Cwat)
%%   cover_bare   bare soil (Cbare)
%%   cover_rock   exposed rock (Crock)
%% Each row sums to 1.  The developed classes split at the middle of their
%% NLCD impervious-cover ranges (open space < 20 %, low 20-49 %, medium
%% 50-79 %, high 80-100 %), the rest being lawn.  Every forest and woody
%% wetland class is tree crown, evergreen included, so a site with
%% evergreen or shrub cover needs its own tree parameters; cultivated
%% crops count as low vegetation.  Edit a copy (writetable / readtable) to
%% change the mapping for a site.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows = { ...
    11, "Open water",                    0, 0,    0,    1, 0, 0
    12, "Perennial ice/snow",            0, 0,    0,    1, 0, 0
    21, "Developed, open space",         0, 0.90, 0.10, 0, 0, 0
    22, "Developed, low intensity",      0, 0.65, 0.35, 0, 0, 0
    23, "Developed, medium intensity",   0, 0.35, 0.65, 0, 0, 0
    24, "Developed, high intensity",     0, 0.10, 0.90, 0, 0, 0
    31, "Barren land",                   0, 0,    0,    0, 1, 0
    41, "Deciduous forest",              1, 0,    0,    0, 0, 0
    42, "Evergreen forest",              1, 0,    0,    0, 0, 0
    43, "Mixed forest",                  1, 0,    0,    0, 0, 0
    52, "Shrub/scrub",                   1, 0,    0,    0, 0, 0
    71, "Grassland/herbaceous",          0, 1,    0,    0, 0, 0
    81, "Pasture/hay",                   0, 1,    0,    0, 0, 0
    82, "Cultivated crops",              0, 1,    0,    0, 0, 0
    90, "Woody wetlands",                1, 0,    0,    0, 0, 0
    95, "Emergent herbaceous wetlands",  0, 1,    0,    0, 0, 0};
T = cell2table(rows,'VariableNames',{'class','name','cover_tree','cover_grass', ...
    'cover_urban','cover_water','cover_bare','cover_rock'});
end
