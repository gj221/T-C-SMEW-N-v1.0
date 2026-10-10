function varargout = Field_Transient_Carry(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Transient_Carry
%%
%% A long field transient runs in segments (Run_Field_Transient). The 2-D
%% driver holds every column's hourly state for the whole window it is
%% given, and each column's ERW chemistry a ten-minute record of it: for
%% Hamden's 936 columns over 2020-2026 that is several hundred GB. Each
%% segment therefore goes on from the end state of the one before it,
%% exactly as a single run would:
%%   - a segment starts at 00:00 of the last but one day of the one before
%%     (Date(1) on a day boundary, so the drivers' calendar days and daily
%%     stamps fall where a single run puts them), and those two days' 48
%%     hours, the overlap, are copied into its first rows. Every hourly
%%     look-back then reads what one run would read: the daily step's 24
%%     hours, and the ERW day's water sums, which also read the hour before
%%     the day (ERW_Daily_Step_2D sums [Tw(1); Tw], shifted by one hour);
%%   - the last three daily rows (the day before the overlap and the two
%%     overlap days) become its rows 1 to 3, its day counter j is 3, and the
%%     time loop resumes at step 49, the first hour after the overlap;
%%   - the columns' carried state (cover, rooting, pH, ponded run-on and
%%     its solute, management indices, the runoff and export ledgers, the
%%     biogeochemical inputs record) is handed on as it stands;
%%   - each ERW column is put on the new segment's ten-minute and daily
%%     axes: its last two days' rows become the first rows and every later row
%%     takes the value the array had at its end, which a segment with no
%%     final flush never wrote, so it is the array's preallocation. The
%%     ten-minute anchor the next ERW day starts from is therefore the
%%     state one run would hold there. An array the ERW days grow row by
%%     row keeps its rows to the anchor, and the CO2 record the first ERW
%%     day lays on its whole axis is laid on the new one.
%% The longest day, which the parameter files find over the days a run is
%% given, is the whole window's in every segment (Run_Field_Transient,
%% MAIN_FRAME_BG_ERW_2D's FIELD_TRANSIENT_LMAX_DAY).
%% The segment before must hold at least two days and leave its last day to
%% the next: no final ERW flush (MAIN_FRAME_BG_ERW_2D,
%% FIELD_TRANSIENT_CONTINUES). A column with
%% particle classes (rock) is refused: its application schedule counts
%% from the start of the run.
%%
%%   D     = Field_Transient_Carry('overlap_days')
%%       the days a segment repeats of the one before (2)
%%   names = Field_Transient_Carry('names')
%%       the driver arrays carried: .hourly, .daily, .column
%%   carry = Field_Transient_Carry('capture',W)
%%       from the driver's workspace values W (the names above, ERW,
%%       lateral_diag, jDay, L_day, Date, NN, NNd, j, and W.whos, the
%%       workspace listing, to check that every time-indexed array is
%%       carried)
%%   snap  = Field_Transient_Carry('snapshot',W)
%%       the same, for a segment that ran its final flush (the end of the
%%       run): what 'digest' identifies
%%   E     = Field_Transient_Carry('erw',E,NN,NNd,Ca)
%%       a carried ERW column on the axes of a segment of NN hours with CO2
%%       series Ca (Field_Transient_Restore)
%%   Field_Transient_Carry('check',carry,Date,N_col)
%%       refuses a carry that does not line up with the segment
%%   hex   = Field_Transient_Carry('digest',x)
%%       SHA-256 of any carry, snapshot or record set
%%   parts = Field_Transient_Carry('parts',x)
%%       {name, SHA-256} of each part of a carry or snapshot, an ERW field
%%       across the columns: the parts two runs disagree on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch action
    case 'overlap_days'
        varargout{1} = overlap_days();
    case 'names'
        varargout{1} = carried_names();
    case 'capture'
        varargout{1} = capture(varargin{1},true);
    case 'snapshot'
        varargout{1} = capture(varargin{1},false);
    case 'erw'
        varargout{1} = erw_on_axes(varargin{:});
    case 'check'
        check(varargin{:});
    case 'digest'
        varargout{1} = digest(varargin{1});
    case 'parts'
        varargout{1} = part_digests(varargin{1});
    otherwise
        error('Field_Transient_Carry:UnknownAction','Unknown action %s.',action);
end
end

%% =======================================================================
function D = overlap_days()
%% Two: the ERW day's water sums read the hour before the day.
D = 2;
end

function N = carried_names()
%% Every state array MAIN_FRAME_BG_ERW_2D allocates by hour (first
%% dimension NN) or by day (NNd), and the per-column state it carries from
%% one step to the next. 'capture' refuses a workspace holding a
%% time-indexed driver array that is not listed here.
N.hourly = [ ...
    {'V','O','Qi_out','WTR','POT','Tdp','Vice','Oice','Qi_in','Tdeb'}, ...
    {'OF','OS','ZWT','Pr_sno','Pr_liq','er','Ts','ra','r_soil','b_soil', ...
     'Tice','Rn','H','QE','Qv','Lpho','EG','ESN','ESN_In','EWAT', ...
     'EIn_urb','EIn_rock','dw_SNO','G','SWE','SND','ros','In_SWE','SP_wc', ...
     'WR_SP','U_SWE','NIn_SWE','dQ','Qfm','t_sls','DQ','DT','In_urb', ...
     'In_rock','SE_rock','SE_urb','Csno','WIS','Lk','f','Rh','Rd','NDVI', ...
     'alp_soil','tau_sno','e_sno','ALB','EK','dQVEG','TsV','HV','QEV', ...
     'Cice','Lk_wat','Lk_rock','EICE','WAT','ICE','ICE_D','WR_IP','NIce', ...
     'Cicew','IP_wc','Csnow','FROCK','Imelt','Smelt','Tdamp','Gfin', ...
     'In_Litter','ELitter','Ws_under','q_runon'}, ...
    {'r_litter','In_H','In_L','OH','OL','Psi_s_H','Psi_s_L','Tdp_H','Tdp_L', ...
     'rap_H','rap_L','Dr_H','Dr_L','T_H','T_L','EIn_H','EIn_L','rb_H','rb_L', ...
     'An_H','Rdark_H','An_L','Rdark_L','Ci_sunH','Ci_sunL','Ci_shdH', ...
     'Ci_shdL','rs_sunH','rs_sunL','rs_shdH','rs_shdL','Vx_H','Vl_H','Vx_L', ...
     'Vl_L','Psi_x_H','Psi_l_H','Psi_x_L','Psi_l_L','gsr_H','Jsx_H','Jxl_H', ...
     'Kleaf_H','Kx_H','gsr_L','Jsx_L','Jxl_L','Kleaf_L','Kx_L','fapar_H', ...
     'fapar_L','SIF_H','SIF_L'}];
veg = {'LAI','B','NPP','Rg','RA','Rmc','ANPP','Rms','Rmr','PHE_S','dflo', ...
    'AgeL','e_rel','SAI','hc','e_relN','BA','Rrootl','Bfac_day','Bfac_week', ...
    'NPPI','TdpI','Sr','Slf','Sfr','Swm','Sll','LAIdead','Rexmy','AgeDL', ...
    'Nreserve','Preserve','Kreserve','Careserve','Mgreserve','Sireserve', ...
    'rNc','rPc','rKc','rCac','rMgc','rSic','NBLeaf','PARI','NBLI','NupI', ...
    'RB','NuLit','Nuptake','NH4uptake','NO3uptake','Puptake','Kuptake', ...
    'Cauptake','Mguptake','Siuptake','FNC','TexC','TexN','TexP','TexK', ...
    'TNIT','TPHO','TPOT','SupN','SupP','SupK','TexCa','TexMg','TexSi', ...
    'TCAL','TMAG','TSIL','SupCa','SupMg','SupSi','ISOIL','AgePl','Tden'};
%% The driver spells three of them without the crown's underscore.
glued = {'Bfac_day','Bfac_week'};
veg_names = cell(1,2*numel(veg));
for k = 1:numel(veg)
    sep = '_';
    if any(strcmp(veg{k},glued)); sep = ''; end
    veg_names{2*k-1} = [veg{k} sep 'H'];
    veg_names{2*k}   = [veg{k} sep 'L'];
end
N.daily = [veg_names, {'TBio_Ht','TBio_Lt','Ccrown_t','BLit','P', ...
    'R_litter','R_microbe','R_litter_sur','R_ew','VOL','N2flx','N2Oflx', ...
    'NITflx','DENITflx','N2Oflx_nit','FertN_d','DepN_d','PH_d','PH_pore_d', ...
    'PH_bulk_d','Dep_d','Dep_erw_d','Min_N','Min_P','R_bacteria','RmycAM', ...
    'RmycEM','Prod_B','Prod_F','BfixN','LitFirEmi','LEAK_NH4','LEAK_NO3', ...
    'LEAK_P','LEAK_K','LEAK_DOC','LEAK_DON','LEAK_DOP','LEAK_Ca','LEAK_Mg', ...
    'LEAK_Si','RexmyI','NavlI','AgrHarNut','SeedInput'}];
N.column = {'Ccrown_c','Cbare_c','ZR95_H_c','ZR95_L_c','RfH_Zs_c','RfL_Zs_c', ...
    'PH_c','PH_bulk_c','Tstm0_c','Bam_c','Bem_c','snow_alb_c','B_IO_c', ...
    'TBio_H_c','TBio_L_c','CcrownFIX','ManIH','ManIL','Tdp_bio_c','VSUM_c', ...
    'VTSUM_c','q_runon_next','pond_mass_area','pond_field_water_mm', ...
    'pond_field_mass_area','erw_pond_umol','erw_pond_field_umol', ...
    'exit_water_vert_m3','exit_water_surface_m3','exit_mass_surface_g', ...
    'exit_mass_surface_erw_umol','exit_mass_surface_erw_field_umol', ...
    'exit_water_field_surface_m3','exit_mass_surface_field_g'};
end

%% =======================================================================
function C = capture(W,strict)
%% The overlap and the state the next segment starts from. strict: the
%% segment left its last day to the next (no final flush).
required = {'Date','NN','NNd','j','ERW','lateral_diag','jDay','L_day'};
missing = required(~isfield(W,required));
if ~isempty(missing)
    error('Field_Transient_Carry:MissingWorkspace', ...
        'The capture needs %s from the driver''s workspace.',strjoin(missing,', '));
end
NN = W.NN; NNd = W.NNd; Date = W.Date(:);
D = overlap_days(); L = 24*D;
if mod(NN,24) ~= 0 || NN < L || numel(Date) ~= NN || ~at_midnight(Date(1))
    error('Field_Transient_Carry:NotWholeDays', ...
        ['A segment that another follows runs at least %d whole days from 00:00 ' ...
         '(%d hours from %s).'],D,NN,datestr(Date(1),'yyyy-mm-dd HH:MM'));
end
if W.j ~= NNd
    error('Field_Transient_Carry:IncompleteDays', ...
        'The segment ended on day %d of %d.',W.j,NNd);
end
names = carried_names();
if isfield(W,'whos'); guard(W.whos,names,NN,NNd,numel(W.ERW)); end
C = struct();
C.overlap_days = D;
C.overlap_hours = L;
C.hour_dn = Date(NN-L+1:NN);
C.N_col = numel(W.ERW);
C.hourly = struct(); C.daily = struct(); C.column = struct();
for k = 1:numel(names.hourly)
    x = field_of(W,names.hourly{k});
    if size(x,1) ~= NN
        error('Field_Transient_Carry:HourlySize','%s has %d rows, not %d hours.', ...
            names.hourly{k},size(x,1),NN);
    end
    C.hourly.(names.hourly{k}) = x(NN-L+1:NN,:);
end
for k = 1:numel(names.daily)
    x = field_of(W,names.daily{k});
    if size(x,1) ~= NNd
        error('Field_Transient_Carry:DailySize','%s has %d rows, not %d days.', ...
            names.daily{k},size(x,1),NNd);
    end
    C.daily.(names.daily{k}) = x(NNd-D:NNd,:);
end
for k = 1:numel(names.column)
    C.column.(names.column{k}) = field_of(W,names.column{k});
end
if numel(W.jDay) < NNd || numel(W.L_day) < NNd
    error('Field_Transient_Carry:IncompleteDays','jDay and L_day do not reach day %d.',NNd);
end
C.jDay = W.jDay(NNd-D:NNd);
C.L_day = W.L_day(NNd-D:NNd);
C.lateral_diag = struct();
for name = fieldnames(W.lateral_diag).'
    x = W.lateral_diag.(name{1});
    if (isnumeric(x) || islogical(x)) && size(x,1) == NNd
        C.lateral_diag.(name{1}) = x(NNd-D:NNd,:);
    end
end
C.ERW = cell(C.N_col,1);
for c = 1:C.N_col
    C.ERW{c} = erw_trim(W.ERW{c},NN,NNd,strict,c);
end
end

function x = field_of(W,name)
if ~isfield(W,name)
    error('Field_Transient_Carry:MissingWorkspace', ...
        'The driver''s workspace holds no %s.',name);
end
x = W.(name);
end

function guard(info,names,NN,NNd,N_col)
%% A time-indexed driver array (first dimension NN or NNd, last N_col) that
%% the lists do not carry would restart from its parameter-file value.
%% Per-column arrays cannot be told apart from them when N_col is NN or NNd.
if any(N_col == [NN NNd]); return; end
known = [names.hourly names.daily names.column {'Date','Pr','Ta','Ws','ea', ...
    'SAD1','SAD2','SAB1','SAB2','Pre','N','Tdew','esat','PARB','PARD','Ds', ...
    'Ca','IrD','Datam','jDay','L_day','field_fraction','A_cell'}];
for k = 1:numel(info)
    sz = info(k).size;
    numeric = any(strcmp(info(k).class,{'double','single','logical', ...
        'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}));
    if ~numeric || numel(sz) < 2 || ~any(sz(1) == [NN NNd]) || sz(end) ~= N_col
        continue
    end
    if ~any(strcmp(info(k).name,known))
        error('Field_Transient_Carry:UncarriedArray', ...
            ['%s [%s] looks like a time-indexed driver array that the segment ' ...
             'carry does not hand on: add it to Field_Transient_Carry.'], ...
            info(k).name,num2str(sz));
    end
end
end

%% =======================================================================
function E = erw_trim(E,NN,NNd,strict,c)
%% One ERW column cut to its last two days: the ten-minute rows of the
%% overlap (the last of them never written without a final flush) and the
%% last three daily rows.
if ~(isstruct(E) && isfield(E,'init') && E.init && isfield(E,'init_done_j2') && E.init_done_j2)
    error('Field_Transient_Carry:ERWNotStarted', ...
        'ERW column %d has not been initialised; a carried segment needs its first ERW day.',c);
end
if isfield(E,'hist') && isfield(E.hist,'names') && ~isempty(E.hist.names)
    error('Field_Transient_Carry:RockColumn', ...
        ['ERW column %d holds particle classes (rock): its application schedule ' ...
         'counts from the start of the run and cannot be carried.'],c);
end
res = E.res; n_t = E.n_t;
D = overlap_days();
keep = 24*D*res;
if n_t ~= NN*res
    error('Field_Transient_Carry:ERWAxis','ERW column %d has %d rows for %d hours.',c,n_t,NN);
end
if strict && E.dm_end ~= (NN-23)*res
    error('Field_Transient_Carry:FinalFlush', ...
        ['ERW column %d stopped at row %d, not at 01:00 of the last day (%d): ' ...
         'the segment ran its final flush and cannot be continued.'], ...
        c,E.dm_end,(NN-23)*res);
end
shift = n_t - keep;
rows = shift+1:n_t;
M = struct('rows',keep,'ten_minute',{{}},'grown',{{}},'daily',{{}},'hourly',{{}});
for part = {'bg','weath','const','plant','TC_const'}
    if ~isfield(E,part{1}); continue; end
    S = E.(part{1});
    for name = fieldnames(S).'
        x = S.(name{1});
        if ~(isnumeric(x) || islogical(x)); continue; end
        if size(x,1) == n_t
            S.(name{1}) = take_rows(x,rows);
            M.ten_minute(end+1,:) = {part{1},name{1}};
        elseif size(x,1) == E.dm_end
            %% written row by row as the days go, not preallocated
            %% (bg.K_CEC): its rows to the anchor, as one run holds them
            S.(name{1}) = take_rows(x,shift+1:E.dm_end);
            M.grown(end+1,:) = {part{1},name{1}};
        end
    end
    E.(part{1}) = S;
end
for name = fieldnames(E.TC).'
    x = E.TC.(name{1});
    if (isnumeric(x) || islogical(x)) && size(x,1) == NNd
        E.TC.(name{1}) = take_rows(x,NNd-D:NNd);
        M.daily{end+1} = name{1};
    elseif (isnumeric(x) || islogical(x)) && size(x,1) == NN
        %% the day's hourly water and CO2, which every ERW day rebuilds
        %% before reading them
        E.TC.(name{1}) = [];
        M.hourly{end+1} = name{1};
    end
end
%% The hourly CO2 stamp at 00:00 of the overlap: the first stamp of the
%% next segment's axis, which holds no hour before it.
M.CO2_h_first = E.CO2_h(NN-24*D+1);
E.dm_end = E.dm_end - shift;
E.field_transient_carry = M;
%% The segment's axes, which the next rebuilds on its own length.
E = rmfield(E,intersect(fieldnames(E),{'t','tm','th','td','CO2_h','t_end','n_t'}));
end

function E = erw_on_axes(E,NN,NNd,Ca)
%% The carried column on a segment of NN hours: its axes rebuilt as
%% ERW_Daily_Step_2D builds them, its ten-minute arrays padded with the
%% value each held at the end, its daily series with zeros.
M = E.field_transient_carry;
E.t_end = NN/24;
E.t = E.dtw:E.dtw:E.t_end;
E.n_t = round((NN/24)/(1/(24*E.res)));
E.th = (0:NN)';
E.td = E.th(1:24:end);
E.tm = ((0:E.n_t-1)')/E.res;
E.CO2_h = [M.CO2_h_first; Ca(:)/22.41];
for k = 1:size(M.ten_minute,1)
    part = M.ten_minute{k,1}; name = M.ten_minute{k,2};
    E.(part).(name) = pad_rows(E.(part).(name),E.n_t,'last');
end
for k = 1:numel(M.daily)
    E.TC.(M.daily{k}) = pad_rows(E.TC.(M.daily{k}),NNd,'zero');
end
for k = 1:numel(M.hourly)
    E.TC.(M.hourly{k}) = zeros(NN,1);
end
%% The CO2 record the first ERW day lays on its whole axis (bg.CO2_atm;
%% no later day reads it, each interpolating its own), laid on this one's
%% as ERW_Daily_Step_2D lays it.
if isfield(E,'bg') && isfield(E.bg,'CO2_atm')
    TT = synchronize(timetable(hours(E.th),E.CO2_h,'VariableNames',{'CO2_atm'}), ...
        timetable(hours(E.td),E.TC.s,'VariableNames',{'s'}),'regular','linear', ...
        'TimeStep',minutes(60/E.res));
    E.bg.CO2_atm = TT.CO2_atm(1:E.n_t);
end
E = rmfield(E,'field_transient_carry');
end

function y = take_rows(x,rows)
sz = size(x);
y = reshape(x(rows,:),[numel(rows) sz(2:end)]);
end

function y = pad_rows(x,n,fill)
sz = size(x);
m = sz(1);
if n < m
    error('Field_Transient_Carry:ShortSegment', ...
        'A segment of %d rows cannot hold the %d carried.',n,m);
end
flat = x(:,:);
if strcmp(fill,'last')
    tail = repmat(flat(m,:),n-m,1);
else
    tail = zeros(n-m,size(flat,2),'like',flat);
end
y = reshape([flat; tail],[n sz(2:end)]);
end

%% =======================================================================
function check(C,Date,N_col)
L = C.overlap_hours;
Date = Date(:);
if numel(Date) < L + 24 || mod(numel(Date),24) ~= 0
    error('Field_Transient_Carry:ShortSegment', ...
        'A carried segment needs the %d-hour overlap and at least one day after it.',L);
end
if any(abs(Date(1:L) - C.hour_dn(:)) > 1/86400)
    error('Field_Transient_Carry:Misaligned', ...
        ['The segment starts at %s; the state carried into it ends with the day ' ...
         'of %s.'],datestr(Date(1),'yyyy-mm-dd HH:MM'),datestr(C.hour_dn(1),'yyyy-mm-dd HH:MM'));
end
if C.N_col ~= N_col
    error('Field_Transient_Carry:ColumnCount','The carry holds %d columns, the segment %d.', ...
        C.N_col,N_col);
end
end

%% =======================================================================
function P = part_digests(x)
%% Every field of a carry or snapshot, one level into its structs; the ERW
%% columns' fields, two levels in, gathered across the columns.
P = cell(0,2);
for name = fieldnames(x).'
    v = x.(name{1});
    if strcmp(name{1},'ERW') && iscell(v) && ~isempty(v)
        for f = fieldnames(v{1}).'
            vals = cellfun(@(e) field_or_empty(e,f{1}),v,'UniformOutput',false);
            if isstruct(vals{1}) && isscalar(vals{1})
                for g = fieldnames(vals{1}).'
                    sub = cellfun(@(s) field_or_empty(s,g{1}),vals,'UniformOutput',false);
                    P(end+1,:) = {sprintf('ERW.%s.%s',f{1},g{1}),digest(sub)}; %#ok<AGROW>
                end
            else
                P(end+1,:) = {['ERW.' f{1}],digest(vals)}; %#ok<AGROW>
            end
        end
    elseif isstruct(v) && isscalar(v)
        for f = fieldnames(v).'
            P(end+1,:) = {[name{1} '.' f{1}],digest(v.(f{1}))}; %#ok<AGROW>
        end
    else
        P(end+1,:) = {name{1},digest(v)}; %#ok<AGROW>
    end
end
end

function v = field_or_empty(s,name)
v = [];
if isstruct(s) && isfield(s,name); v = s.(name); end
end

function hex = digest(x)
bytes = getByteStreamFromArray(x);
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(bytes(:),'int8'));
raw = typecast(md.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end

function tf = at_midnight(dn)
h = round(dn*24);
tf = abs(dn*24 - h) < 1e-6 && mod(h,24) == 0;
end
