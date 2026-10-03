function tr = Lateral_Track_Record(mode, tr, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral_Track_Record
%%
%% Flux tracking for the field spin-up (Lian et al. 2025, Fig. 3c).  While
%% the 2-D field run exchanges water and dissolved pools between cells,
%% record for a subset of cells what entered and left each of them, so the
%% cell can later be run as a plot-scale column that receives the same
%% lateral exchange (Lateral_Replay_Prepare / _Water / _Solute).
%%
%% Recorded per tracked cell
%%   q_runon [NN x 1]    overland run-on offered to the column at step i
%%                       [mm/h] (the routed runoff of step i-1)
%%   dtheta  [NN x ms]   net lateral Darcy change of liquid water content
%%                       per layer at step i [m3 m-3]
%%   retain, input [NNd x N_species]
%%                       the day's lateral exchange of each mobile pool as
%%                       one affine map, composed exactly from the hourly
%%                       upwind steps of that day:
%%                           T_end = retain .* T_start + input
%%                       retain is the share of the cell's own pool that
%%                       did not leave it (subsurface export and soil
%%                       solute carried off with saturation-excess runoff),
%%                       input the mass that arrived from neighbours
%%                       [g m-2] (subsurface inflow and re-infiltrated
%%                       run-on solute).  Replayed this way the cell's
%%                       export stays proportional to its own, evolving
%%                       pool -- the upwind scheme with the neighbours'
%%                       concentrations frozen at their recorded values.
%%   day_dn  [NNd x 1]   calendar day (datenum) of each daily block; row 1
%%                       is the initial state and carries no exchange.
%%
%% Modes
%%   tr = Lateral_Track_Record('init',[],cells,N_col,NN,NNd,ms, ...
%%            species,pool_idx,lateral_fraction,Date)
%%   tr = Lateral_Track_Record('day',tr,j,date_dn)
%%   tr = Lateral_Track_Record('runon',tr,i,q_runon_col)
%%   tr = Lateral_Track_Record('subsurface',tr,i,j,S)
%%        S: topo, P_cells (pools before the transport), dtheta_water,
%%           dM_edge, Q_solute, theta_before, A_cell, Zlayer_m,
%%           active_fraction, dth
%%   tr = Lateral_Track_Record('surface',tr,j,S)
%%        S: P_cells (pools after the subsurface step), pond_in_area,
%%           soil_removed_area, Rd_mm, theta, A_cell, Zlayer_m,
%%           active_fraction
%%   rec = Lateral_Track_Record('cell',tr,k)   one tracked cell's record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(mode)
    case 'init'
        [cells,N_col,NN,NNd,ms,species,pool_idx,lateral_fraction,Date] = varargin{:};
        cells = cells(:);
        if isempty(cells) || any(~isfinite(cells)) || any(cells ~= fix(cells)) || ...
                any(cells < 1) || any(cells > N_col) || numel(unique(cells)) ~= numel(cells)
            error('Lateral_Track_Record:BadCells', ...
                'Tracked cells must be distinct indices within 1..%d.',N_col);
        end
        if ~all(cellfun(@isscalar,pool_idx))
            error('Lateral_Track_Record:BadPools', ...
                'Every laterally mobile species must map to one pool.');
        end
        n = numel(cells);
        n_sp = numel(species);
        tr = struct();
        tr.cells = cells;
        tr.species = species;
        tr.pool_idx = cellfun(@(x) x,pool_idx);
        tr.lateral_fraction = lateral_fraction(cells,:);       %% [n x n_sp]
        tr.hour_dn = Date(:);
        tr.q_runon = zeros(NN,n);
        tr.dtheta = zeros(NN,ms,n);
        tr.day_dn = nan(NNd,1);
        tr.retain = ones(NNd,n_sp,n);
        tr.input = zeros(NNd,n_sp,n);

    case 'day'
        [j,date_dn] = varargin{:};
        tr.day_dn(j) = floor(date_dn);

    case 'runon'
        [i,q_runon_col] = varargin{:};
        tr.q_runon(i,:) = reshape(q_runon_col(tr.cells),1,[]);

    case 'subsurface'
        [i,j,S] = varargin{:};
        c = tr.cells;
        tr.dtheta(i,:,:) = reshape(S.dtheta_water(c,:).',[1 size(S.dtheta_water,2) numel(c)]);
        %% Per-edge masses [g] over the step, signed a -> b.  The two
        %% directions are split layer by layer before summing: one face
        %% can carry water downslope in one layer and upslope in another.
        N_col = size(S.P_cells,1);
        n_edge = numel(S.topo.edge_a);
        a_to_b = reshape(sum(max(S.dM_edge,0),2),n_edge,[]);   %% [N_edges x n_sp]
        b_to_a = reshape(sum(max(-S.dM_edge,0),2),n_edge,[]);
        Ia = sparse(S.topo.edge_a(:),(1:n_edge).',1,N_col,n_edge);
        Ib = sparse(S.topo.edge_b(:),(1:n_edge).',1,N_col,n_edge);
        out_g = full(Ia(c,:)*a_to_b + Ib(c,:)*b_to_a);         %% [n x n_sp]
        in_g  = full(Ib(c,:)*a_to_b + Ia(c,:)*b_to_a);
        area = S.A_cell(c); area = area(:);
        %% The share of the active-zone water that left each tracked cell:
        %% the upwind export of every species is this share of its mobile
        %% part (the donor's layer outflow capped at the layer's water,
        %% exactly as Solute_Lateral_Transport limits it).  It stands in
        %% where a pool is empty (see retention below).
        Vw = S.theta_before .* (S.A_cell(:)*S.Zlayer_m(:).') .* S.active_fraction;
        Qa = abs(S.Q_solute)*S.dth;                            %% [N_edges x ms] m3
        out_w = Ia*(Qa.*(S.Q_solute > 0)) + Ib*(Qa.*(S.Q_solute < 0));
        out_w = full(out_w(c,:));                              %% [n x ms] m3
        f_water = sum(min(out_w,Vw(c,:)),2)./max(sum(Vw(c,:),2),eps);
        f_water(sum(Vw(c,:),2) <= 0) = 0;
        T = S.P_cells(c,tr.pool_idx);                          %% [n x n_sp] g m-2
        tr = compose(tr,j,retention(out_g./area,T,f_water.*tr.lateral_fraction), ...
            in_g./area);

    case 'surface'
        [j,S] = varargin{:};
        c = tr.cells;
        T = S.P_cells(c,tr.pool_idx);
        %% fallback share for an empty pool: the saturation-excess runoff
        %% drains this share of the active-zone water
        W_act = sum(S.theta(c,:).*S.active_fraction(c,:).*S.Zlayer_m(:).',2);  %% [m]
        f_water = min(1,(S.Rd_mm(c)/1000)./max(W_act,eps));
        f_water(W_act <= 0) = 0;
        tr = compose(tr,j,retention(S.soil_removed_area(c,:),T, ...
            f_water(:).*tr.lateral_fraction),S.pond_in_area(c,:));

    case 'cell'
        k = varargin{1};
        rec = struct();
        rec.hour_dn = tr.hour_dn;
        rec.q_runon = tr.q_runon(:,k);
        rec.dtheta = tr.dtheta(:,:,k);
        rec.day_dn = tr.day_dn;
        rec.retain = tr.retain(:,:,k);
        rec.input = tr.input(:,:,k);
        rec.species = tr.species;
        rec.pool_idx = tr.pool_idx;
        tr = rec;

    otherwise
        error('Lateral_Track_Record:BadMode','Unknown mode %s.',mode);
end
end

%% =======================================================================
function a = retention(out_area,T,f_empty)
%% The retained share of each pool: one minus what left over what was
%% there.  An empty pool loses nothing in the field run, but a replayed
%% cell whose pool is not empty at that hour must still lose the share
%% the water carried, so the water-based share stands in.
a = 1 - f_empty;
has = T > 0;
a(has) = 1 - out_area(has)./T(has);
a = min(1,max(0,a));
end

function tr = compose(tr,j,a,b)
%% Append one hourly map to day j's: T -> a.*T + b after the day so far.
n = numel(tr.cells);
A = reshape(tr.retain(j,:,:),[],n).';
B = reshape(tr.input(j,:,:),[],n).';
A = a.*A;
B = a.*B + b;
tr.retain(j,:,:) = reshape(A.',[1 size(A,2) n]);
tr.input(j,:,:)  = reshape(B.',[1 size(B,2) n]);
end
