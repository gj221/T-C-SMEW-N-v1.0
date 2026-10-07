function test_event_deposition()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Event deposition (Deposition_Day, Build_Site_Deposition).
%%   1  without an event record the day's deposition is the constant rate
%%      and B_IO is returned untouched (every site but Hamden)
%%   2  with one: dry every day, wet on the day's precipitation, the year's
%%      row; a span's mean rate; a year the record lacks and a negative
%%      precipitation are errors, never a silent zero; the day's record
%%      stays double beside the global maps' single-precision P and Si
%%   3  the Hamden record: one row per year 2014-2026; complete NTN years
%%      keep their own concentration, the rest pool the three nearest
%%      complete years; years beyond TDep's grids take the mean of the
%%      three nearest; its long-term rates are the record applied to the
%%      forcing's precipitation, and Phase 2's ring mean is the same sum
%%   4  the frames apply it: the plot-scale and 2-D frames hand
%%      BIOGEO_UNIT the day's deposition, Phase 2 (running the row's
%%      vegetation override, as the frames do) the ring's mean and the
%%      Phase 3 budgets, from the first computed day, what the frame applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_directory = fileparts(mfilename('fullpath'));
code = fileparts(test_directory);
addpath(code,fullfile(code,'Inputs','Preprocess'));

%% --- 1 --------------------------------------------------------------------
B = struct('DepN',2e-3,'DepP',1e-5,'DepK',1e-4,'DepCa',2e-4,'DepMg',1e-4,'DepSi',3e-6);
[B_day,dep] = Deposition_Day(B,[2020 6 1 1 0 0],25);
assert(isequal(B_day,B) && isequal(dep,[2e-3 1e-5 1e-4 2e-4 1e-4 3e-6]), ...
    'Without an event record the constant rates must stand.');

%% --- 2 --------------------------------------------------------------------
E = B;
E.Dep_events = struct('year',[2019;2020], ...
    'wet_g_m2_per_mm',[2e-4 2e-5 5e-5 3e-5; 3e-4 1e-5 4e-5 2e-5], ...
    'dry_g_m2_d',[1e-3 2e-6 7e-6 3e-6; 2e-3 3e-6 8e-6 4e-6], ...
    'elements',{{'N','K','Ca','Mg'}});
[~,dry_day] = Deposition_Day(E,[2020 3 4 1 0 0],0);
[D_wet,wet_day] = Deposition_Day(E,[2020 3 5 1 0 0],12.5);
assert(max(abs(dry_day([1 3 4 5])-[2e-3 3e-6 8e-6 4e-6])) < 1e-15 && ...
    max(abs(wet_day([1 3 4 5])-([2e-3 3e-6 8e-6 4e-6]+12.5*[3e-4 1e-5 4e-5 2e-5]))) < 1e-15, ...
    'Dry deposition every day, wet on the day''s precipitation.');
assert(isequal(dry_day([2 6]),[1e-5 3e-6]) && D_wet.DepN == wet_day(1) && ...
    isequal(D_wet.Dep_events,E.Dep_events),'P and Si keep their rates.');
[~,year_2019] = Deposition_Day(E,[2019 12 31 1 0 0],10);
assert(abs(year_2019(1)-(1e-3+10*2e-4)) < 1e-15,'The date''s year selects the row.');
[~,span] = Deposition_Day(E,[2020 1 1 0 0 0],366*3.2,366);
assert(abs(span(1)-(2e-3+3.2*3e-4)) < 1e-15,'A span gives its mean daily rate.');
E_single = E; E_single.DepP = single(1e-5); E_single.DepSi = single(3e-6);
[~,mixed] = Deposition_Day(E_single,[2020 3 5 1 0 0],12.5);
assert(isa(mixed,'double') && isequal(mixed([1 3 4 5]),wet_day([1 3 4 5])), ...
    'A single-precision P or Si rate must not round the day''s record.');
expect_error(@() Deposition_Day(E,[2021 1 1 1 0 0],0),'Deposition_Day:NoYear');
expect_error(@() Deposition_Day(E,[2020 1 1 1 0 0],-1),'Deposition_Day:BadPrecipitation');
expect_error(@() Deposition_Day(E,[2020 1 1 1 0 0],NaN),'Deposition_Day:BadPrecipitation');

%% --- 3 --------------------------------------------------------------------
S = load(fullfile(code,'Inputs','Deposition_Hamden.mat'),'Deposition');
D = S.Deposition;
ev = D.events;
assert(isequal(ev.year(:).',2014:2026) && all(ev.wet_g_m2_per_mm(:) > 0) && ...
    all(ev.dry_g_m2_d(:) > 0) && isequal(ev.elements,{'N','K','Ca','Mg'}), ...
    'The Hamden record holds every forcing year with positive rates.');
W = D.wet_mg_L;
own = W.complete;
assert(all(W.source(own) == "NTN "+string(W.year(own))) && ...
    W.source(W.year == 2026) == "NTN 2023+2024+2025 pooled" && ...
    W.source(W.year == 2014) == "NTN 2016+2017+2018 pooled", ...
    'Complete years keep their own concentration, the others pool the nearest three.');
assert(all(W.sample_span(own) >= 0.9 & W.valid_precip_fraction(own) >= 0.75) && ...
    all(~own(W.year == 2020)),'NADP''s completeness criteria decide.');
Y = D.dry_kg_ha_yr;
assert(Y.source(Y.year == 2025) == "TDep 2022+2023+2024 mean" && ...
    abs(Y.N(Y.year == 2025)-mean(Y.N(ismember(Y.year,2022:2024)))) < 1e-12, ...
    'A year beyond the grids takes the mean of the three nearest.');
F = load(fullfile(code,'Inputs','Data_Hamden.mat'),'Date','Pr');
forcing_year = year(datetime(F.Date,'ConvertFrom','datenum'));
total = zeros(1,4); days = 0;
for y = 2015:2025
    n = datenum(y+1,1,1)-datenum(y,1,1);
    [~,rate] = Deposition_Day(struct('DepN',0,'DepP',0,'DepK',0,'DepCa',0, ...
        'DepMg',0,'DepSi',0,'Dep_events',ev),[y 1 1 0 0 0], ...
        sum(F.Pr(forcing_year == y)),n);
    total = total+n*rate([1 3 4 5]);
    days = days+n;
end
assert(max(abs(total/days-D.mean_g_m2_d)./D.mean_g_m2_d) < 1e-12, ...
    'The long-term rates are the record applied to the forcing''s precipitation.');
A = D.annual;
assert(abs(mean(A.total_N_kg_ha)-mean(A.wet_N_kg_ha)-mean(A.dry_N_kg_ha)) < 1e-12 && ...
    mean(A.dry_N_kg_ha) > mean(A.wet_N_kg_ha),'Dry N outweighs wet N at the parcel.');

%% --- 4 --------------------------------------------------------------------
frame_1d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_10mins.m'));
frame_2d = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_2D.m'));
spinup = fileread(fullfile(code,'Run_BG_Spinup.m'));
assert(contains(frame_1d,'[B_IO_bgc_day,Dep_d(j,:)] = Deposition_Day(B_IO,Datam(i,:),sum(Pr(pdind))*dth);') && ...
    contains(frame_1d,'B_IO_bgc_day,Datam(i,:),FireA,0,...'), ...
    'The plot-scale frame must hand BIOGEO_UNIT the day''s deposition.');
assert(contains(frame_2d,'[B_IO_bgc_day,Dep_d(j,:,c)] = Deposition_Day(B_IO_c_now, ...') && ...
    contains(frame_2d,'B_IO_bgc_day,Datam(i,:),FireA,0, ...'), ...
    'The 2-D frame must hand BIOGEO_UNIT each cell''s day''s deposition.');
assert(contains(spinup,'B_IO_p2 = rmfield(B_IO_p2,''Dep_events'');') && ...
    contains(spinup,'p3_dep_sum = sum(Dep_d(p3_rows,[1 3 4 5 6]),1);') && ...
    contains(spinup,'p3_terms.deposition = p3_dep_sum(p3_x);'), ...
    'Phase 2 must take the ring''s mean and the Phase 3 budgets what was applied.');
assert(contains(spinup,"ismember('pft_override',MOD_PARAM_input.Properties.VariableNames)") && ...
    contains(spinup,'p3_r0 = 2;') && contains(spinup,'p3_rows = p3_r0+1:j;') && ...
    contains(spinup,'sum(P(p3_r0,p3_soil_idx{p3_x}))'), ...
    ['Phase 2 must run the row''s vegetation override (else it never sees the ' ...
     'event record), and the Phase 3 budgets start on the first computed day.']);

fprintf('test_event_deposition: all checks passed\n');
end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
