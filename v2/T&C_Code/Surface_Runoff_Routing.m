function [q_deliver_mm, C_deliver, exit_water_m3, exit_mass, diag, ...
    exit_by_cell_m3, exit_mass_by_cell] = ...
    Surface_Runoff_Routing(topo, runoff_mm, C_runoff, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surface_Runoff_Routing
%%
%% Route ponded overland runoff one cell downslope on the DEM topology,
%% carrying its dissolved load.  Returns (a) the depth + concentration
%% delivered to each cell as its NEXT vertical step's re-infiltration
%% opportunity (q_runon) and (b) the water + mass that leaves the domain
%% at the topographic outlet(s).
%%
%% MODEL
%%   Each cell's ponded runoff this step is partitioned among its
%%   downslope neighbours by the Quinn (1991) multiple-flow-direction
%%   weights (from Attach_Topography_To_Topology) and delivered ONE cell
%%   hop.  A share that lands on an outlet cell leaves the domain; a share
%%   that lands on an interior cell becomes that cell's next-step q_runon,
%%   where the Richards / infiltration solver decides how much
%%   re-infiltrates and how much becomes runoff again (routed on the
%%   following step).  A pit cell keeps its own runoff.  Solute moves with
%%   the water, upwind from the contributing cell.
%%
%%   One hop per step deliberately under-propagates fast overland flow
%%   (real travel times are sub-hourly), which lags -- but conserves --
%%   streamflow timing while giving every downslope cell a re-infiltration
%%   chance.  It is mass-exact.
%%
%%   No topography attached -> no weights, no outlet: every cell keeps its
%%   own runoff (q_deliver == runoff), the pre-topography behaviour.
%%
%% INPUTS
%%   topo        topology; uses cell_z_m, edge_overland_weight_ab/ba,
%%               cell_is_outlet, cell_area_m2, edge_a/edge_b, N_edges
%%   runoff_mm   [N_cells x 1] ponded runoff depth this step [mm]
%%   C_runoff    [N_cells x N_species] its aqueous concentration
%%               [mass / m3 water].  [] => no solute tracking.
%%   dt_h        hydrologic timestep [h] (interface symmetry / future use)
%%
%% OUTPUTS
%%   q_deliver_mm      [N_cells x 1] depth delivered to each cell [mm]
%%   C_deliver         [N_cells x N_species] its concentration [mass/m3]
%%   exit_water_m3     scalar water leaving at outlets [m3]
%%   exit_mass         [N_species x 1] solute mass leaving at outlets
%%   diag              .edge_flux_m3 [N_edges x 1] overland water per edge
%%                     .pit_cells    interior cells with nowhere to drain
%%   exit_by_cell_m3   [N_cells x 1] outlet water attributed per cell
%%   exit_mass_by_cell [N_cells x N_species] outlet mass per cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
runoff_mm = runoff_mm(:);
if numel(runoff_mm) ~= N_cells || any(~isfinite(runoff_mm)) || any(runoff_mm < -1e-9)
    error('Surface_Runoff_Routing:BadRunoff', ...
        'runoff_mm must hold one finite nonnegative depth per cell.');
end
runoff_mm = max(0,runoff_mm);
track_solute = nargin >= 3 && ~isempty(C_runoff);
if track_solute
    N_sp = size(C_runoff,2);
    if size(C_runoff,1) ~= N_cells || any(~isfinite(C_runoff(:))) || any(C_runoff(:) < 0)
        error('Surface_Runoff_Routing:BadConc', ...
            'C_runoff must be [N_cells x N_species], finite and nonnegative.');
    end
else
    N_sp = 0;
end
if nargin < 4; dt_h = 1; end %#ok<NASGU>

A = topo.cell_area_m2(:);
has_topo = isfield(topo,'edge_overland_weight_ab') && ...
    isfield(topo,'cell_z_m') && ~isempty(topo.cell_z_m);

diag = struct('edge_flux_m3',zeros(topo.N_edges,1),'pit_cells',[]);
exit_water_m3      = 0;
exit_by_cell_m3    = zeros(N_cells,1);
exit_mass_by_cell  = zeros(N_cells,max(N_sp,1));

if ~has_topo || all(runoff_mm == 0)
    q_deliver_mm = runoff_mm;
    if track_solute
        C_deliver = C_runoff; exit_mass = zeros(N_sp,1);
        exit_mass_by_cell = exit_mass_by_cell(:,1:N_sp);
    else
        C_deliver = zeros(N_cells,0); exit_mass = zeros(0,1);
        exit_mass_by_cell = zeros(N_cells,0);
    end
    return
end

is_outlet = false(N_cells,1);
if isfield(topo,'cell_is_outlet'); is_outlet = logical(topo.cell_is_outlet(:)); end
w_ab = topo.edge_overland_weight_ab(:);
if isfield(topo,'edge_overland_weight_ba')
    w_ba = topo.edge_overland_weight_ba(:);
else
    w_ba = zeros(topo.N_edges,1);
end

pond_m3 = (runoff_mm/1000).*A;                       % this step's own runoff
mass    = zeros(N_cells,max(N_sp,1));
if track_solute; mass(:,1:N_sp) = pond_m3 .* C_runoff; end

delivered_m3   = zeros(N_cells,1);
delivered_mass = zeros(N_cells,max(N_sp,1));
pit = false(N_cells,1);

for c = 1:N_cells
    Vc = pond_m3(c);
    if Vc <= 0; continue; end
    if is_outlet(c)
        exit_water_m3      = exit_water_m3 + Vc;
        exit_by_cell_m3(c) = exit_by_cell_m3(c) + Vc;
        if track_solute
            exit_mass_by_cell(c,1:N_sp) = exit_mass_by_cell(c,1:N_sp) + mass(c,1:N_sp);
        end
        continue
    end
    ea_out = find(topo.edge_a == c & w_ab > 0);
    eb_out = find(topo.edge_b == c & w_ba > 0);
    targets = [topo.edge_b(ea_out); topo.edge_a(eb_out)];
    weights = [w_ab(ea_out);        w_ba(eb_out)];
    edges   = [ea_out;              eb_out];
    esign   = [ones(numel(ea_out),1); -ones(numel(eb_out),1)];   % a->b positive
    wsum = sum(weights);
    if wsum <= 0
        pit(c) = true;
        delivered_m3(c) = delivered_m3(c) + Vc;
        if track_solute
            delivered_mass(c,1:N_sp) = delivered_mass(c,1:N_sp) + mass(c,1:N_sp);
        end
        continue
    end
    weights = weights/wsum;
    if track_solute; Cc = mass(c,1:N_sp)/max(Vc,eps); end
    for t = 1:numel(targets)
        nt = targets(t);
        dV = Vc*weights(t);
        diag.edge_flux_m3(edges(t)) = diag.edge_flux_m3(edges(t)) + esign(t)*dV;
        if is_outlet(nt)
            exit_water_m3       = exit_water_m3 + dV;
            exit_by_cell_m3(nt) = exit_by_cell_m3(nt) + dV;
            if track_solute
                exit_mass_by_cell(nt,1:N_sp) = exit_mass_by_cell(nt,1:N_sp) + dV*Cc;
            end
        else
            delivered_m3(nt) = delivered_m3(nt) + dV;
            if track_solute
                delivered_mass(nt,1:N_sp) = delivered_mass(nt,1:N_sp) + dV*Cc;
            end
        end
    end
end

q_deliver_mm = 1000*delivered_m3./A;
if track_solute
    C_deliver = zeros(N_cells,N_sp);
    nz = delivered_m3 > 0;
    C_deliver(nz,:) = delivered_mass(nz,1:N_sp)./delivered_m3(nz);
    exit_mass         = sum(exit_mass_by_cell(:,1:N_sp),1).';
    exit_mass_by_cell = exit_mass_by_cell(:,1:N_sp);
else
    C_deliver = zeros(N_cells,0);
    exit_mass = zeros(0,1);
    exit_mass_by_cell = zeros(N_cells,0);
end
diag.pit_cells = find(pit);

%% Mass-balance guard.
in_water  = sum((runoff_mm/1000).*A);
out_water = sum(delivered_m3) + exit_water_m3;
if abs(in_water-out_water) > 1e-8*max(1,in_water)
    error('Surface_Runoff_Routing:WaterImbalance', ...
        'Routed water (%.6g m3) != input (%.6g m3).',out_water,in_water);
end
if track_solute
    in_mass  = sum(pond_m3.*C_runoff,1).';
    out_mass = sum(delivered_mass(:,1:N_sp),1).' + exit_mass;
    bad = abs(in_mass-out_mass) > 1e-7*max(1,max(abs(in_mass)));
    if any(bad)
        error('Surface_Runoff_Routing:MassImbalance', ...
            'Routed solute species %d does not conserve.',find(bad,1));
    end
end
end
