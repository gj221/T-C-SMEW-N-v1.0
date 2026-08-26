function bg_state = build_bg_state(P, O, RexmyI, NavlI, NupI_L, ...
    Nreserve_L, Preserve_L, Kreserve_L, Careserve_L, Mgreserve_L, Sireserve_L, ...
    Nreserve_H, Preserve_H, Kreserve_H, Careserve_H, Mgreserve_H, Sireserve_H, id_location)
daily_idx = last_nonzero_row(P, 'P');
hydro_idx = last_finite_row(O, 'O');
plant_idx = min(daily_idx, size(NavlI,1));
bg_state = struct();
bg_state.B_out = P(daily_idx,:);
bg_state.O_out = O(hydro_idx,:);
bg_state.RexmyI_out = RexmyI(min(daily_idx,size(RexmyI,1)),:);
bg_state.NavlI_out = NavlI(plant_idx,:);
bg_state.NupI_L_out = NupI_L(min(daily_idx,size(NupI_L,1)),:,:);
bg_state.Nreserve_L_out = Nreserve_L(plant_idx,:);
bg_state.Preserve_L_out = Preserve_L(plant_idx,:);
bg_state.Kreserve_L_out = Kreserve_L(plant_idx,:);
bg_state.Careserve_L_out = Careserve_L(plant_idx,:);
bg_state.Mgreserve_L_out = Mgreserve_L(plant_idx,:);
bg_state.Sireserve_L_out = Sireserve_L(plant_idx,:);
bg_state.Nreserve_H_out = Nreserve_H(plant_idx,:);
bg_state.Preserve_H_out = Preserve_H(plant_idx,:);
bg_state.Kreserve_H_out = Kreserve_H(plant_idx,:);
bg_state.Careserve_H_out = Careserve_H(plant_idx,:);
bg_state.Mgreserve_H_out = Mgreserve_H(plant_idx,:);
bg_state.Sireserve_H_out = Sireserve_H(plant_idx,:);
bg_state.created_on = datestr(now, 'yyyy-mm-dd HH:MM:SS');

function idx = last_nonzero_row(A, name)
if isempty(A) || ~isnumeric(A)
    error('last_nonzero_row:InvalidBGState', ...
        'Cannot build BG state: %s is empty or non-numeric.', name);
end
A2 = reshape(A, size(A,1), []);
ok = any(isfinite(A2) & abs(A2) > 0, 2);
idx = find(ok, 1, 'last');
if isempty(idx)
    error('last_nonzero_row:InvalidBGState', ...
        'Cannot build BG state: %s has no finite non-zero rows.', name);
end
end

function idx = last_finite_row(A, name)
if isempty(A) || ~isnumeric(A)
    error('last_finite_row:InvalidBGState', ...
        'Cannot build BG state: %s is empty or non-numeric.', name);
end
A2 = reshape(A, size(A,1), []);
idx = find(all(isfinite(A2), 2), 1, 'last');
if isempty(idx)
    error('last_finite_row:InvalidBGState', 'Cannot build BG state: %s has no finite rows.', name);
end
end
bg_state.source = ['TC_SMEW_N_run spinup ' char(string(id_location))];
end
