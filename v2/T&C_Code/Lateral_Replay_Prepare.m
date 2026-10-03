function [rep,SAB1,SAB2,SAD1,SAD2,PARB,PARD,SvF] = Lateral_Replay_Prepare( ...
    R,Date,Datam,NN,NNd,ms,SAB1,SAB2,SAD1,SAD2,PARB,PARD,SvF,DeltaGMT,Lon,Lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral_Replay_Prepare
%%
%% Map a cell's lateral-exchange record from the 2-D field run
%% (Lateral_Track_Record, assembled over the spin-up years by
%% Run_Field_Lateral_Tracking) onto the hours and days of one plot-scale
%% run of that cell, and give the column the topographic radiation the
%% cell received in the field run (Topographic_Radiation_Projection).
%% This is the plot-scale half of the flux-tracking spin-up of Lian et al.
%% (2025, Fig. 3d); MAIN_FRAME_BG_ERW_10mins calls it once, after the
%% parameter file, when the caller defines LATERAL_REPLAY.
%%
%% R (the record)
%%   hour_dn [Nh x 1]  q_runon [Nh x 1]  dtheta [Nh x ms]
%%   day_dn  [Nd x 1]  retain [Nd x Nsp] input [Nd x Nsp]  pool_idx [1 x Nsp]
%%   required          true: every step and day of this run must be in the
%%                     record (a spin-up year); false: hours outside it
%%                     carry no exchange (the two-day ERW seed check)
%%   slope, aspect_deg, svf, solar_window   (optional) the cell's terrain;
%%                     absent for a flat field run
%%
%% rep: the per-step series the driver hooks read, and the ledgers they
%% fill: water_mm / water_clip_mm [NN x 1], solute_in_g / solute_out_g
%% [NNd x Nsp] (g m-2 entering / leaving each pool per day).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

need = {'hour_dn','q_runon','dtheta','day_dn','retain','input','pool_idx','required'};
missing = need(~isfield(R,need));
if ~isempty(missing)
    error('Lateral_Replay:BadRecord','The lateral record lacks %s.',strjoin(missing,', '));
end
Nh = numel(R.hour_dn); Nd = numel(R.day_dn); n_sp = numel(R.pool_idx);
if numel(R.q_runon) ~= Nh || ~isequal(size(R.dtheta),[Nh ms]) || ...
        ~isequal(size(R.retain),[Nd n_sp]) || ~isequal(size(R.input),[Nd n_sp])
    error('Lateral_Replay:BadRecord', ...
        'The lateral record is inconsistent with itself or with %d soil layers.',ms);
end
if any(~isfinite(R.q_runon)) || any(R.q_runon < 0) || any(~isfinite(R.dtheta(:))) || ...
        any(~isfinite(R.retain(:))) || any(R.retain(:) < 0) || any(R.retain(:) > 1) || ...
        any(~isfinite(R.input(:))) || any(R.input(:) < 0)
    error('Lateral_Replay:BadRecord', ...
        'Run-on and inputs must be finite and non-negative, retained shares within [0,1].');
end
if any(R.pool_idx < 1) || any(R.pool_idx > 79) || any(R.pool_idx ~= fix(R.pool_idx))
    error('Lateral_Replay:BadRecord','Pool indices must address the 79-pool state.');
end
required = logical(R.required);

rep = struct();
rep.required = required;
rep.pool_idx = R.pool_idx(:).';
rep.q_runon = zeros(NN,1);
rep.dtheta = zeros(NN,ms);
rep.water_mm = zeros(NN,1);
rep.water_clip_mm = zeros(NN,1);
rep.solute_in_g = zeros(NNd,n_sp);
rep.solute_out_g = zeros(NNd,n_sp);

%% Steps: matched on the whole hour, so a datenum rounding difference
%% between the two runs cannot drop an hour.
[hit,loc] = ismember(round(Date(:)*24),round(R.hour_dn(:)*24));
steps = 2:NN;
if required && ~all(hit(steps))
    error('Lateral_Replay:NotCovered', ...
        'The lateral record does not cover %d of the %d steps from %s.', ...
        nnz(~hit(steps)),numel(steps),datestr(Date(2),'yyyy-mm-dd HH:MM'));
end
use = steps(hit(steps));
rep.q_runon(use) = R.q_runon(loc(use));
rep.dtheta(use,:) = R.dtheta(loc(use),:);

%% Days: the driver opens day j at the step whose hour is 1, exactly as
%% the 2-D run did, so its j counter is rebuilt here the same way.
rep.retain = ones(NNd,n_sp);
rep.input = zeros(NNd,n_sp);
day_key = round(R.day_dn(:));
j = 1; n_missing = 0;
for i = steps
    if Datam(i,4) == 1
        j = j + 1;
        r = find(day_key == floor(Date(i)),1);
        if isempty(r)
            n_missing = n_missing + 1;
        else
            rep.retain(j,:) = R.retain(r,:);
            rep.input(j,:) = R.input(r,:);
        end
    end
end
if required && n_missing > 0
    error('Lateral_Replay:NotCovered', ...
        'The lateral record does not cover %d of the days from %s.', ...
        n_missing,datestr(Date(2),'yyyy-mm-dd'));
end

%% The radiation the cell received in the field run.
if all(isfield(R,{'slope','aspect_deg','svf','solar_window'}))
    for i = steps
        [SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i)] = ...
            Topographic_Radiation_Projection(SAB1(i),SAB2(i),SAD1(i),SAD2(i), ...
            PARB(i),PARD(i),Datam(i,:),DeltaGMT,Lon,Lat, ...
            R.solar_window(1),R.solar_window(2),R.slope,R.aspect_deg,R.svf);
    end
    SvF = R.svf;
end
end
