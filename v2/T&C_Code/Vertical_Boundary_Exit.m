function exit = Vertical_Boundary_Exit(topo, exit, Q_edge_lat, ...
    C_aq_cell, leak_bot_cell, C_aq_bot_cell, dt_h, day_idx, N_days, Q_outer_lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical_Boundary_Exit
%%
%% Mass-balance accumulator for water and solute leaving the 24-cell
%% domain through (i) the perimeter-drain faces of the buffer cells
%% (non-zero only when opt_bc='open_outflow') and (ii) the bottom of the
%% soil column via Leakage_Bottom.
%%
%% This is bookkeeping only -- no process model.  Called once per
%% hydrologic timestep after Darcy_Lateral_Topology and
%% Solute_Lateral_Transport.
%%
%% Inputs
%%   topo            24-cell topology struct
%%   exit            running accumulator struct (see below)
%%   Q_edge_lat      [N_edges x ms]  per-edge lateral water flux from
%%                                    Darcy_Lateral_Topology [m3/h]
%%   C_aq_cell       [N_cells x ms x N_species] aqueous conc [mass/m3]
%%   leak_bot_cell   [N_cells x 1]   vertical leakage at the bottom of
%%                                    the soil column for each cell
%%                                    (from Leakage_Bottom.m) [mm/h]
%%   C_aq_bot_cell   [N_cells x N_species] aqueous conc in the bottom
%%                                    layer at each cell [mass/m3]
%%   dt_h            timestep [h]
%%   day_idx         (optional) current day index 1..N_days.  If supplied
%%                   together with N_days, the accumulator also fills
%%                   exit.NO3_leach_ts_gN_m2 (daily time series, gN/m2/d).
%%                   If omitted, only the cumulative total is updated.
%%   N_days          (optional) total number of simulated days (for
%%                   preallocating the daily time series on first call).
%%   Q_outer_lat     (optional) [N_cells x ms] perimeter water outflow
%%                   rate [m3/h], positive leaving the domain. If omitted,
%%                   lateral exit is assumed zero.
%%
%% Outputs
%%   exit  struct with cumulative fields:
%%      .water_lat_main   [m3]
%%      .water_lat_buffer [m3]
%%      .water_vert_main   [m3]
%%      .water_vert_buffer [m3]
%%      .solute_lat_main   [N_species x 1] [mass]
%%      .solute_lat_buffer [N_species x 1]
%%      .solute_vert_main  [N_species x 1]
%%      .solute_vert_buffer[N_species x 1]
%%      .NO3_leach_gN      [gN] cumulative NO3-N leaving via all paths
%%                         (species index 7 in the 12-species transport list:
%%                          Ca=1, Mg=2, K=3, Na=4, Al=5, Si=6, NO3=7, NH4=8,
%%                          DOC=9, DIC=10, PO4=11, An=12)
%%      .NO3_leach_ts_gN_m2 [N_days x 1] daily time series of NO3-N
%%                         leaving the domain normalised by total domain
%%                         area [gN m-2 d-1].  Only populated when the
%%                         caller supplies day_idx + N_days; otherwise
%%                         left absent to avoid preallocating arrays the
%%                         caller will never query.
%%
%% Note
%%   Under the strict zero-flux Neumann ('closed_perimeter') lateral BC,
%%   lateral exit fluxes are zero by construction. When the caller uses
%%   'open_outflow', this accumulator additionally records the perimeter
%%   water and solute export passed in via Q_outer_lat. Vertical exit via
%%   Leakage_Bottom is always accounted for here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells   = topo.N_cells;
N_species = size(C_aq_cell, 3);

if nargin < 8; day_idx = []; end
if nargin < 9; N_days  = []; end
if nargin < 10 || isempty(Q_outer_lat)
    Q_outer_lat = zeros(N_cells, size(C_aq_cell, 2));
end
track_ts = ~isempty(day_idx) && ~isempty(N_days);

% Initialise accumulator on first call
if isempty(exit) || ~isstruct(exit)
    exit.water_lat_main    = 0;
    exit.water_lat_buffer  = 0;
    exit.water_vert_main   = 0;
    exit.water_vert_buffer = 0;
    exit.solute_lat_main    = zeros(N_species,1);
    exit.solute_lat_buffer  = zeros(N_species,1);
    exit.solute_vert_main   = zeros(N_species,1);
    exit.solute_vert_buffer = zeros(N_species,1);
    exit.NO3_leach_gN       = 0;   %% cumulative NO3-N [gN] leaving domain
    exit.used_grid_routing = 0;  %% 0 = single-column, 1 = multi-cell grid active
    if track_ts
        exit.NO3_leach_ts_gN_m2 = zeros(N_days,1);   %% [gN m-2 d-1]
    end
end
% If this is a non-first call and caller newly starts tracking the daily
% time series, lazily allocate the field now.
if track_ts && ~isfield(exit, 'NO3_leach_ts_gN_m2')
    exit.NO3_leach_ts_gN_m2 = zeros(N_days,1);
end

NO3_IDX = 7;   %% NO3 position in 12-species list (Ca,Mg,K,Na,Al,Si,NO3,...)
NO3_G_PER_UMOL = 14e-6;  %% routed NO3 is stored as umol N; summaries use gN

is_main = strcmp(topo.cell_region, 'main');

%% Domain area for normalising the daily time series (NO3_leach_ts is
%% expressed per m2 of total domain so that buffer and main contributions
%% are area-weighted correctly).
domain_area_m2 = sum(topo.cell_area_m2);   %% [m2]

%% --- 1. Lateral water + solute exit through outer faces -------------
% Under 'closed_perimeter' (zero-flux Neumann) Q_outer_lat is zero by
% construction. When the caller uses 'open_outflow', the same accumulator
% fields capture the perimeter export without changing the downstream
% mass-balance interface.
%#ok<*INUSD>
for c = 1:N_cells
    for k = 1:size(C_aq_cell, 2)
        Q_out_h = Q_outer_lat(c,k);                % m3/h
        if Q_out_h <= 0
            continue;
        end
        Vw = Q_out_h * dt_h;                       % m3 over dt_h
        if is_main(c)
            exit.water_lat_main = exit.water_lat_main + Vw;
        else
            exit.water_lat_buffer = exit.water_lat_buffer + Vw;
        end
        for s = 1:N_species
            m = Vw * C_aq_cell(c,k,s);
            if is_main(c)
                exit.solute_lat_main(s) = exit.solute_lat_main(s) + m;
            else
                exit.solute_lat_buffer(s) = exit.solute_lat_buffer(s) + m;
            end
            if s == NO3_IDX
                m_gN = m * NO3_G_PER_UMOL;
                exit.NO3_leach_gN = exit.NO3_leach_gN + m_gN;
                if track_ts
                    exit.NO3_leach_ts_gN_m2(day_idx) = ...
                        exit.NO3_leach_ts_gN_m2(day_idx) + m_gN / domain_area_m2;
                end
            end
        end
    end
end

%% --- 2. Vertical water + solute exit via Leakage_Bottom -------------
% leak_bot_cell is in mm/h, cell areas in m2; convert mm/h * m2 -> L/h,
% then divide by 1000 -> m3/h.
for c = 1:N_cells
    A_c    = topo.cell_area_m2(c);                  % m2
    qv_mmh = leak_bot_cell(c);                      % mm/h
    if qv_mmh <= 0; continue; end
    Vw_h   = (qv_mmh/1000) * A_c;                   % m3/h
    Vw     = Vw_h * dt_h;                            % m3 over dt_h
    if is_main(c)
        exit.water_vert_main = exit.water_vert_main + Vw;
    else
        exit.water_vert_buffer = exit.water_vert_buffer + Vw;
    end
    for s = 1:N_species
        m = Vw * C_aq_bot_cell(c,s);                % mass
        if is_main(c)
            exit.solute_vert_main(s) = exit.solute_vert_main(s) + m;
        else
            exit.solute_vert_buffer(s) = exit.solute_vert_buffer(s) + m;
        end
        if s == NO3_IDX
            m_gN = m * NO3_G_PER_UMOL;
            exit.NO3_leach_gN = exit.NO3_leach_gN + m_gN;   %% [gN]
            if track_ts
                %% Per-area daily flux [gN m-2 d-1].  m_gN is the step mass
                %% [gN] over dt_h; dividing by dt_h*domain_area_m2 gives
                %% an instantaneous rate, but since we want a daily total
                %% normalised by area we just add the per-m2 mass into
                %% the day bucket (dt sums to 24 h -> day).
                exit.NO3_leach_ts_gN_m2(day_idx) = ...
                    exit.NO3_leach_ts_gN_m2(day_idx) + m_gN / domain_area_m2;
            end
        end
    end
end

end
