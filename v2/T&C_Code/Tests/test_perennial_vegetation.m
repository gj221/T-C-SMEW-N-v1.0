function test_perennial_vegetation()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A perennial stand carried from one run of a driver to the next
%% (Perennial_Vegetation_States, Perennial_Vegetation_Capture, _Apply).
%%   1  the list is every vegetation state the drivers read from the
%%      previous day, with the litter and exudates the soil takes the next
%%      day and the canopy the first day's hours see; the 10-minute driver
%%      creates each before the parameter files run, so their first row
%%      stands
%%   2  capture takes the run's last day; apply sets the next run's first
%%      row and nothing else, creates a state the parameter files have not
%%      (the 2-D driver's per-cell workspace), refuses an incomplete or
%%      non-finite state, and leaves a cold start alone
%%   3  the spin-up carries the stand from Phase 3 on (Phase 1 has no
%%      nutrient uptake to refill its reserves), through every checkpoint
%%      and the state it hands off, and the loaders and field drivers pass
%%      it on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_directory = fileparts(mfilename('fullpath'));
code = fileparts(test_directory);
v2 = fileparts(code);
addpath(code);
[names,has_depth] = Perennial_Vegetation_States();
assert(numel(names) == numel(unique(names)) && numel(has_depth) == numel(names), ...
    'Each carried state is listed once, with its depth.');

%% --- 1 --------------------------------------------------------------------
frame_1d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_10mins.m'));
frame_2d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_2D.m'));
read_prev = union(previous_day_reads(frame_1d),previous_day_reads(frame_2d));
%% Soil, exchange and ERW states (carried as the soil's state, B_out ...)
%% and the stand-structure states of the optional demography (OPT_VCA).
not_vegetation = {'Ca','K','Mg','NH4','NO3','P','PO4','Si','NavlI','RexmyI', ...
    'Ccrown_t','AgePl_H','AgePl_L','Tden_H','Tden_L'};
missing = setdiff(read_prev,[names not_vegetation]);
assert(isempty(missing),'Read from the previous day but not carried: %s.',strjoin(missing,', '));
first_day = {'e_rel_H','e_rel_L','e_relN_H','e_relN_L','Rrootl_H','Rrootl_L'};
extra = setdiff(names,[read_prev first_day]);
assert(isempty(extra),'Carried but read by neither driver: %s.',strjoin(extra,', '));
param_at = strfind(frame_1d,'run(PARAM_IC);');
for k = 1:numel(names)
    made = regexp(frame_1d,['(^|[\s;,])' names{k} '\s*=\s*(zeros|ones)\(NNd,cc(,\w+)?\)'], ...
        'tokens','once','lineanchors');
    at = regexp(frame_1d,['(^|[\s;,])' names{k} '\s*=\s*(zeros|ones)\(NNd,cc'],'once','lineanchors');
    assert(~isempty(at) && at < param_at(1), ...
        '%s must exist before the parameter files run.',names{k});
    assert(has_depth(k) == ~isempty(made{3}), ...
        'The depth of %s must match the driver''s array.',names{k});
end

%% --- 2 --------------------------------------------------------------------
NNd = 6; cc = 2; j = NNd;
depth = struct('B',8,'NupI',6,'PARI',3,'NuLit',6,'ISOIL',27,'Rexmy',3);
rng(7);
W = struct(); W0 = struct();
for k = 1:numel(names)
    sz = [NNd cc];
    if has_depth(k); sz = [NNd cc depth.(regexprep(names{k},'_[HL]$',''))]; end
    W.(names{k}) = rand(sz);
    W0.(names{k}) = zeros(sz);
end
Veg = capture_from(W,j);
ok = true;
for k = 1:numel(names)
    ok = ok && isequal(Veg.(names{k}),W.(names{k})(j,:,:)) && size(Veg.(names{k}),1) == 1;
end
assert(ok && isequal(sort(fieldnames(Veg)),sort(names(:))), ...
    'Capture must take every state''s last day.');
A = apply_to(W0,Veg);
ok = true;
for k = 1:numel(names)
    ok = ok && isequal(A.(names{k})(1,:,:),Veg.(names{k})) && ...
        ~any(A.(names{k})(2:end,:,:),'all');
end
assert(ok,'Apply must set the first row and nothing else.');
assert(isequal(apply_to(W,[]),W),'Without a carried stand nothing changes.');
E = apply_to(struct(),Veg);
ok = true;
for k = 1:numel(names); ok = ok && isequal(E.(names{k}),Veg.(names{k})); end
assert(ok,'Apply must create a state the parameter files have not.');
%% The 2-D driver's per-cell parameter workspace preallocates nothing:
%% x(1,:) = 0 there makes a scalar, which the carried row replaces.
R = struct('TdpI_H',0,'FNC_L',[0 0],'B_H',zeros(1,2,8));
R = apply_to(R,Veg);
assert(isequal(R.TdpI_H,Veg.TdpI_H) && isequal(R.FNC_L,Veg.FNC_L) && isequal(R.B_H,Veg.B_H), ...
    'A state set as its first row alone must take the carried row.');
bad = rmfield(Veg,'NuLit_L');
expect_error(@() apply_to(W0,bad),'Perennial_Vegetation_Apply:IncompleteState');
bad = Veg; bad.B_H(1,1,3) = NaN;
expect_error(@() apply_to(W0,bad),'Perennial_Vegetation_Apply:InvalidState');

%% --- 3 --------------------------------------------------------------------
spinup = fileread(fullfile(code,'Run_BG_Spinup.m'));
loader = fileread(fullfile(v2,'spinup','Run_Automatic_Site_Spinup.m'));
tracking = fileread(fullfile(v2,'spinup','Run_Field_Lateral_Tracking.m'));
field = fileread(fullfile(v2,'spinup','Run_Field_Cell_Spinup.m'));
pft = fileread(fullfile(code,'Parameters','PFT_Hamden_Deciduous_Grass.m'));
%% Phase 1 has no biogeochemistry, so no uptake refills a carried stand's
%% nutrient reserves: the stand is carried from Phase 3 on (each Phase-3
%% year and the handoff extension), through its checkpoint, a warm start, a
%% handoff-only run and the state handed off.
assert(numel(strfind(spinup,'Perennial_Vegetation_Capture;')) == 2 && ...
    ~contains(spinup,'p1_avgs.Veg_out') && ~contains(spinup,'p2_bg_state.Veg_out = Veg_out;') && ...
    contains(spinup,'A perennial stand is not carried through Phase 1') && ...
    contains(spinup,'''rsd''; ''Veg_out'';') && ...
    contains(spinup,'Veg_out = p2_bg_state.Veg_out;') && ...
    contains(spinup,'Veg_out = tmp.Veg_out;') && ...
    contains(spinup,'Veg_out = p3_handoff_source.Veg_out;') && ...
    contains(spinup,'bg_state_no_erw.Veg_out = Veg_out;'), ...
    ['The spin-up must carry a perennial stand from Phase 3 on, through its ' ...
     'checkpoint, the handoff extension and the state it hands off.']);
assert(numel(strfind(loader,'Veg_out = bg_state.Veg_out;')) == 2, ...
    'A seed or a transient must start a perennial stand from the state handed in.');
assert(contains(frame_2d,'''spinup_erw_handoff_diagnostic'',''Veg_out''};') && ...
    contains(frame_2d,'ctx.Veg_out = Veg_out;'), ...
    'The 2-D driver must hand each cell''s stand to its vegetation file.');
assert(contains(tracking,'ColCfg(c).BG_STATE.Veg_out = state(c).Veg_out;') && ...
    contains(tracking,'V.(veg_names{k}) = x(j,:,:,c);') && ...
    contains(tracking,'V.(veg_names{k}) = x(j,:,c);') && ...
    contains(tracking,'vegetation_vector(seeds(c))') && ...
    numel(strfind(field,'''Veg_out'',{stand_of(')) == 2, ...
    'The field tracking must start each cell from its seed''s stand and carry it between segments.');
at_apply = strfind(pft,'Perennial_Vegetation_Apply;');
at_zatm = strfind(pft,'zatm = ');
assert(contains(pft,'PFT_perennial = true;') && ~isempty(at_apply) && at_apply(1) < at_zatm(1), ...
    'The Hamden stand is perennial and applies a carried state before its reference height.');

fprintf('test_perennial_vegetation: all checks passed\n');
end

function names = previous_day_reads(src)
t = regexp(src,'([A-Za-z_]\w*)\(j-1','tokens');
names = unique(cellfun(@(x) x{1},t,'UniformOutput',false));
end

function Veg_out = capture_from(W,j) %#ok<INUSD>
for f = fieldnames(W).'
    eval([f{1} ' = W.(f{1});']);
end
Perennial_Vegetation_Capture;
end

function W = apply_to(W,Veg)
if ~isempty(Veg)
    Veg_out = Veg; %#ok<NASGU>
end
for f = fieldnames(W).'
    eval([f{1} ' = W.(f{1});']);
end
Perennial_Vegetation_Apply;
names = Perennial_Vegetation_States();
for k = 1:numel(names)
    W.(names{k}) = eval(names{k});
end
end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier),'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
