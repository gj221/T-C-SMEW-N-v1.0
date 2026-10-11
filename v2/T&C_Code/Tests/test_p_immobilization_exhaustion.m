function test_p_immobilization_exhaustion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Microbes that take all the mineral P must leave the pool at zero, not
%% below it. Hamden's class seed under the layer's leaching stopped in
%% Phase 3 loop 11 (23 July 2017) with pool 43 at -2.8e-9 g m^-2: fungi
%% regrowing on a 100 kg C m^-2 POC pile wanted 0.064 g P m^-2 d^-1 from
%% a pool and its day's inputs holding 4e-5. The deposition map holds P in
%% single precision, so the split of that demand between mineral P and SOM
%% was single, and its rounding (4e-8 of the demand) overdrew the pool the
%% demand had emptied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
code_dir = fileparts(here);
addpath(code_dir);
F = load(fullfile(here,'Fixtures','p_immobilization_hamden_20170723.mat'), ...
    'fixture');
a = F.fixture.args;
soil = a{3}*a{4}*1000; %% g soil m^-2, as BIOGEO_UNIT passes it on
tol = 16*eps(1)*soil;  %% Biogeo_Pool_Update's rounding band at an empty pool
assert(isa(a{16},'single') && a{2}(43) < 1e-9 && ...
    a{2}(19)/a{2}(40) > a{50}.R_CP_fun, ['The fixture must hold single ' ...
    'deposition P, an empty mineral pool and P-poor fungi.']);

%% The day itself: the pool ends at zero and the update accepts it.
[dB,~,~,~,~,~,~,~,Min_P] = BIOGEOCHEMISTRY_DYNAMIC3(a{:});
P = Biogeo_Pool_Update(a{2},dB,soil);
assert(P(43) == 0, ...
    'Microbes that take all the mineral P must leave the pool at zero.');
assert(Min_P < 0 && -Min_P > a{2}(43), ...
    'The microbes take the pool and its day''s inputs.');
P_in = sum(a{5}(13:15))+a{52}(5)+a{15}+double(a{16})+a{42};
P_out = a{27}+a{35}+a{38};
residual = sum(dB(35:47))-(P_in-P_out);
assert(abs(residual) < 1e-10, ...
    'The emptied pool''s day must conserve phosphorus (%.3g g m^-2).',residual);

%% With the deposition in double the same day never overdrew.
b = a;
b{16} = double(b{16});
dB_double = BIOGEOCHEMISTRY_DYNAMIC3(b{:});
P_double = Biogeo_Pool_Update(b{2},dB_double,soil);
assert(P_double(43) == 0,'Double deposition P must also empty the pool exactly.');

%% Any demand that empties the pool, and any pool size it empties: the
%% single-precision split overdrew 133 of these 301 fungal states and 87
%% of these 201 pools, by up to 1.6e-6 g m^-2.
for fungi = linspace(150,900,301)
    b = a;
    b{2}(19) = fungi;
    dB_k = BIOGEOCHEMISTRY_DYNAMIC3(b{:});
    assert(b{2}(43)+dB_k(43) >= -tol, ...
        'Fungi of %.4g g C m^-2 overdraw the mineral P pool by %.3g g m^-2.', ...
        fungi,-(b{2}(43)+dB_k(43)));
end
for pool = logspace(-6,-1,201)
    b = a;
    b{2}(43) = pool; b{2}(19) = 2000; b{2}(40) = 1; b{27} = 0;
    dB_k = BIOGEOCHEMISTRY_DYNAMIC3(b{:});
    P_k = b{2}+dB_k;
    assert(P_k(43) >= -16*eps(max(1,abs(P_k(43))/soil))*soil, ...
        'Emptying %.3g g m^-2 of mineral P overdraws it by %.3g g m^-2.', ...
        pool,-P_k(43));
end

%% A pool the demand does not empty keeps the arithmetic it always had
%% (v0.9's Kantola holds 180 g P m^-2 and never comes close).
b = a;
b{2}(43) = 10;
[~,~,~,~,~,~,~,~,Min_P_full] = BIOGEOCHEMISTRY_DYNAMIC3(b{:});
assert(isa(Min_P_full,'single') && Min_P_full < 0, ...
    'A mineral pool the microbes do not empty must keep its single-precision split.');

fprintf('PASS: microbes that empty the mineral P leave it at zero.\n');
end
