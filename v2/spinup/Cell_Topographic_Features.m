function [X, names] = Cell_Topographic_Features(topo, resolved)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell_Topographic_Features
%%
%% The covariate matrix the field spin-up's random forest is trained on
%% (Lian et al. 2025): elevation, slope, aspect as a unit vector,
%% profile and plan curvature, TWI, plus whichever resolved soil columns
%% the agronomic table carries.  One definition, shared by the
%% reconstruction (Reconstruct_Cell_Pools_RF) and its validation
%% (Run_RF_Validation_Experiment), so the two can never disagree on the
%% feature space.
%%
%%   X      [N_cells x n_feat]
%%   names  {1 x n_feat}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needed = {'cell_z_m','cell_slope','cell_aspect_deg','cell_curv_prof', ...
    'cell_curv_plan','cell_twi'};
missing = needed(~isfield(topo,needed));
if ~isempty(missing)
    error('Cell_Topographic_Features:NoTopography', ...
        ['The topology carries no %s; attach a DEM with ' ...
         'Attach_Topography_To_Topology first.'], strjoin(missing,', '));
end
X = [topo.cell_z_m(:), topo.cell_slope(:), cosd(topo.cell_aspect_deg(:)), ...
    sind(topo.cell_aspect_deg(:)), topo.cell_curv_prof(:), ...
    topo.cell_curv_plan(:), topo.cell_twi(:)];
names = {'elevation_m','slope','aspect_cos','aspect_sin', ...
    'curv_profile','curv_plan','twi'};
if ~isempty(resolved) && istable(resolved)
    for col = {'sand_fraction','clay_fraction','OM_fraction','bulk_pH', ...
            'CEC_cmolc_kg'}
        hit = find(strcmpi(resolved.Properties.VariableNames,col{1}),1);
        if ~isempty(hit)
            X = [X, double(resolved.(resolved.Properties.VariableNames{hit}))]; %#ok<AGROW>
            names{end+1} = col{1}; %#ok<AGROW>
        end
    end
end
end
