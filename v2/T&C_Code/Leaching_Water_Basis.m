function basis = Leaching_Water_Basis(MOD_PARAM_input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The water the dissolved pools of the biogeochemical layer leach with,
%% from a site's MOD_PARAM row (column leaching_water).
%%
%%   'column'  T&C-BG's leaching, and v0.9's: the day's drainage out of the
%%             bottom of the soil column (Lk) against the water of the
%%             whole column (VSUM), so a pool B leaches B*Lk/VSUM.
%%   'layer'   v1.0: the day's drainage out of the bottom of the
%%             biogeochemical layer (Biogeo_Layer_Drainage) against the
%%             layer's own water, Zbio*theta_bio, the volume the ERW's
%%             pore-water concentrations (cov_moll) and the 2-D lateral
%%             transport (Biogeo_Lateral_Concentration) already use.
%%
%% The pools are held in the layer, so under 'column' they leave at their
%% concentration diluted by the whole column's water (about 1000 mm on
%% the default 5 m column against 50-80 mm in a 300 mm layer), carried by
%% the drainage of the column's bottom, which lags the layer's by weeks
%% to months.  Uptake, in both cases, sees the layer's mobile water (VT).
%%
%% A row without the column, or with an empty, missing or 'column' value,
%% gives 'column', so it runs as v0.9 did.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basis = 'column';
raw = [];
if istable(MOD_PARAM_input)
    if ismember('leaching_water',MOD_PARAM_input.Properties.VariableNames)
        raw = MOD_PARAM_input.leaching_water;
    end
elseif isstruct(MOD_PARAM_input) && isfield(MOD_PARAM_input,'leaching_water')
    raw = MOD_PARAM_input.leaching_water;
end
if iscell(raw)
    if isempty(raw); return; end
    raw = string(raw);
end
if isempty(raw)
    return
end
%% A workbook column left empty reads as NaN.
if isnumeric(raw) && all(isnan(raw(:)))
    return
end
if ~(ischar(raw) || isstring(raw))
    error('Leaching_Water_Basis:Unknown', ...
        'leaching_water must be ''column'' or ''layer''.');
end
raw = lower(strtrim(string(raw)));
raw(ismissing(raw)) = "";
raw(raw == "") = "column";
if any(raw ~= raw(1))
    error('Leaching_Water_Basis:Ambiguous', ...
        'The parameter rows disagree on leaching_water (%s).',strjoin(unique(raw),', '));
end
switch raw(1)
    case "column"
    case "layer"
        basis = 'layer';
    otherwise
        error('Leaching_Water_Basis:Unknown', ...
            'leaching_water must be ''column'' or ''layer'', not ''%s''.',raw(1));
end
end
