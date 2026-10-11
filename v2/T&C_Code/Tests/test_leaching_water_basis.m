function test_leaching_water_basis
%% The leaching_water option (Leaching_Water_Basis): the default leaches
%% the biogeochemical layer's dissolved pools as v0.9 does, with the
%% column's bottom drainage over the column's water; 'layer' leaches them
%% with the drainage out of the bottom of the layer (Biogeo_Layer_Drainage)
%% over the layer's own water, the same way in both drivers, the ERW and
%% the spin-up's stand-alone phase.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2 = fileparts(code_dir);
addpath(code_dir);

%% 1. The row's column, default v0.9.
T = table(5,'VariableNames',{'Number'});
assert(strcmp(Leaching_Water_Basis(T),'column'),'A row without the column must leach as v0.9.');
assert(strcmp(Leaching_Water_Basis([]),'column'),'No row must leach as v0.9.');
for v = {"",string(missing),"column"," Column ",NaN,{''}}
    T.leaching_water = v{1};
    assert(strcmp(Leaching_Water_Basis(T),'column'), ...
        'An empty, missing or ''column'' value must leach as v0.9.');
end
T.leaching_water = " Layer ";
assert(strcmp(Leaching_Water_Basis(T),'layer'),'''layer'' must select the layer.');
assert(strcmp(Leaching_Water_Basis(struct('leaching_water','layer')),'layer'), ...
    'A struct row must be read the same way.');
T.leaching_water = "profile";
assert_throws(@() Leaching_Water_Basis(T),'Leaching_Water_Basis:Unknown');
T.leaching_water = 3;
assert_throws(@() Leaching_Water_Basis(T),'Leaching_Water_Basis:Unknown');
T2 = table([1;2],["layer";"column"],'VariableNames',{'Number','leaching_water'});
assert_throws(@() Leaching_Water_Basis(T2),'Leaching_Water_Basis:Ambiguous');

%% 2. The layer's drainage is the step's water balance of its layers.
%% Three layers, the biogeochemistry spanning the top two.
bio = [0.4 0.6 0];
dz = [100 150 250];
V_prev = [20 30 60];  Vice_prev = [1 0 0];
V_now  = [22 31 61];  Vice_now  = [0 0 0];
f = 9; Rd = 0.5; EG = 0.4; EvL = [0.8 0.2 0];
Jsx_H = 0.6; RfH = [0.5 0.3 0.2];
Jsx_L = [0.3 0.1]; RfL = [0.6 0.4 0; 0.2 0.2 0.5];   % second row sums to 0.9
Qi_in = [0.2 0 0]; Qi_out = [0.1 0.1 0];
q = Biogeo_Layer_Drainage(bio,V_prev,Vice_prev,V_now,Vice_now,f,Rd,EG, ...
    Jsx_H,Jsx_L,Qi_in,Qi_out,EvL,RfH,RfL,1,1,1);
J_top = 0.6*0.8 + 0.3*1.0 + 0.1*(0.4/0.9);   % each crown's profile rescaled to its Jsx
dS = (22+31+0) - (20+30+1);
expected = f - Rd - (EG*1.0 + J_top) + (0.2-0.2) - dS;
assert(abs(q-expected) < 1e-12,'Layer drainage must close the layer''s step balance.');
%% Over every layer it is T&C's own column balance: the bottom leakage
%% plus the step's closure error (CK1 in the drivers).
q_all = Biogeo_Layer_Drainage(ones(1,3),V_prev,Vice_prev,V_now,Vice_now,f,Rd,EG, ...
    Jsx_H,Jsx_L,Qi_in,Qi_out,EvL,RfH,RfL,1,1,1);
Lk = 0.7;
CK1 = f + sum(V_prev-V_now) + sum(Vice_prev-Vice_now) - EG - Lk - sum(Qi_out) - Rd ...
    - sum(Jsx_L) - sum(Jsx_H) + sum(Qi_in);
assert(abs(q_all-(Lk+CK1)) < 1e-12,'Over the whole column it must be Lk plus T&C''s closure error.');
%% Slope and rock scale the terms as T&C's balance does, and a per-hour
%% rate comes back whatever the step.
q_s = Biogeo_Layer_Drainage(bio,V_prev,Vice_prev,V_now,Vice_now,f,Rd,EG, ...
    Jsx_H,Jsx_L,Qi_in,Qi_out,EvL,RfH,RfL,1.02,0.9,1);
expected_s = f*1.02*0.9 - Rd - (EG + J_top) - 1.02*0.9*dS;
assert(abs(q_s-expected_s) < 1e-12,'Asur and Ared must enter as in T&C''s balance.');
q_2h = Biogeo_Layer_Drainage(bio,V_prev,Vice_prev,V_now,Vice_now,f,Rd,EG, ...
    Jsx_H,Jsx_L,Qi_in,Qi_out,EvL,RfH,RfL,1,1,2);
assert(abs(q_2h-(2*f - Rd - 2*(EG + J_top) - dS)/2) < 1e-12, ...
    'The drainage must come back as a rate per hour.');
%% Water drawn up into a drying layer is negative drainage.
q_up = Biogeo_Layer_Drainage(bio,V_prev,Vice_prev,V_prev,Vice_prev,0,0,0, ...
    1,[0 0],[0 0 0],[0 0 0],EvL,RfH,RfL,1,1,1);
assert(abs(q_up+0.8) < 1e-12,'Uptake from a layer with constant storage is water rising into it.');

%% 3. At the layer basis a pool leaches at its concentration in the layer.
par = Biogeochemistry_Parameter(1);
pool = zeros(1,79); pool(32) = 3; pool(60) = 2;
W = 300*0.25;
[~,leak_no3,~,~,~,~,~,leak_ca] = Biogeo_Leakage(pool,6,W,par);
assert(abs(leak_no3-3*6/W) < 1e-12 && abs(leak_ca-2*6/W) < 1e-12, ...
    'Leaching must carry the layer''s concentration times its drainage.');

%% 4. Both drivers, the ERW step and the spin-up route the option the same way.
one = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
two = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
step = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
spin = fileread(fullfile(code_dir,'Run_BG_Spinup.m'));
sub = fileread(fullfile(v2,'run','Run_One_Subfield.m'));
field = fileread(fullfile(v2,'spinup','Run_Field_Cell_Spinup.m'));
wb = fileread(fullfile(code_dir,'Read_Parameter_Workbook.m'));
assert(contains(one,"strcmp(Leaching_Water_Basis(MOD_PARAM_input),'layer')") && ...
    contains(two,"strcmp(Leaching_Water_Basis(PAR{c}.MOD_PARAM_input),'layer')"), ...
    'Both drivers must read the row''s leaching_water.');
assert(contains(one,'Se_fc,V_leach,VTSUM,') && contains(one,'Lk_leach,sum(T_H(pdind,:),1)') && ...
    contains(two,'Se_fc,V_leach,VTSUM,') && contains(two,'Lk_leach,sum(T_H(pdind,:,c),1)'), ...
    'Both drivers must hand BIOGEO_UNIT the leaching water and its volume.');
assert(contains(one,'Lk_leach = sum(Lk(pdind));') && contains(one,'V_leach = VSUM;') && ...
    contains(two,'Lk_leach = sum(Lk(pdind,1,c));') && contains(two,'V_leach = VSUM;'), ...
    'The default must be v0.9''s column drainage over the column''s water.');
assert(contains(one,'V_leach = Zbio*sum(Bio_Zs.*mean(O(pdind,:),1));') && ...
    contains(two,'V_leach = p.Zbio*sum(p.Bio_Zs.*mean(O(pdind,:,c),1));'), ...
    'The layer basis must use the layer''s own water, the ERW''s cov_moll.');
%% The drainage comes from the vertical step, before lateral exchange.
i1 = strfind(one,'Lk_bio(i) = Biogeo_Layer_Drainage(');
i2 = strfind(one,'Lateral_Replay_Water(lat_rep,i,');
assert(isscalar(i1) && isscalar(i2) && i1 < i2, ...
    'The 1-D layer drainage must precede the lateral replay.');
j1 = strfind(two,'Lk_bion(1,1,c) = Biogeo_Layer_Drainage(');
j2 = strfind(two,'LATERAL COUPLING');
assert(isscalar(j1) && ~isempty(j2) && j1 < j2(1) && ...
    contains(two,'Vn(1,:,c),Vicen(1,:,c),fn(1,1,c),Rdn(1,1,c),EGn(1,1,c)') && ...
    contains(two,'if leach_layer; Lk_bio(i,:,:) = Lk_bion; end'), ...
    'The 2-D layer drainage must read the hour''s column buffers and precede the lateral routing.');
assert(contains(one,'ERW_TC.Lw = Lk_bio/1000;') && contains(one,'ERW_TC.f_Lw_bio(j) = Lk_tot/V_leach;') && ...
    contains(two,'Lk_erw = Lk_bio; V_erw_c = V_leach_c;') && ...
    contains(two,'O(:,:,c), Lk_erw(:,1,c), f(:,1,c)') && ...
    contains(two,'Tdp_bio_c(c), V_erw_c(c), VTSUM_c(c)'), ...
    'The ERW must leach with the same water as the biogeochemistry.');
assert(contains(step,"if Lk_tot < 0 && strcmp(Leaching_Water_Basis(MOD_PARAM_input),'layer')"), ...
    'The 2-D ERW step must stop at a day of upward drainage.');
assert(contains(spin,'p1_VSUM_acc = p1_VSUM_acc + mean(V_leach_out(jrows),''omitnan'');') && ...
    contains(spin,'p1_Lk_acc = p1_Lk_acc + mean(Lk_leach_out(jrows),''omitnan'');') && ...
    contains(spin,'ERW_Mean_Water_Loss(T_H, T_L, Lk_bio, Bio_Zs, RfH_Zs, RfL_Zs)'), ...
    'The spin-up''s Phase 2 and seed water loss must follow the layer basis.');
assert(contains(sub,'''LeachingWater'',[]') && contains(sub,'__leachlayer') && ...
    contains(sub,'MOD_PARAM_input.leaching_water = repmat("layer",height(MOD_PARAM_input),1);'), ...
    'Run_One_Subfield must carry LeachingWater into the row and tag its results.');
assert(contains(field,'''co2_file'',''leaching_water''}') && contains(wb,'''leaching_water''}'), ...
    'A site''s MOD columns and the workbook must accept leaching_water.');


%% 5. DOC mobility: the column's path keeps v0.9's 0.04 everywhere; under
%% 'layer' BIOGEO_UNIT leaches DOC at 0.004 and nothing else moves.
assert(par.aDOC == 0.04,'The default DOC mobility must stay at v0.9''s 0.04.');
pool(12) = 30; pool(13) = 50;
[~,~,~,~,leak_doc] = Biogeo_Leakage(pool,6,W,par);
assert(abs(leak_doc-0.04*80*6/W) < 1e-12, ...
    'The column''s DOC leaching must be aDOC*DOC*Lk/V with aDOC 0.04.');
par_layer = par; par_layer.aDOC = 0.004;
[~,leak_no3_l,~,~,leak_doc_l] = Biogeo_Leakage(pool,6,W,par_layer);
assert(abs(leak_doc_l-0.004*80*6/W) < 1e-12 && leak_no3_l == leak_no3, ...
    'Under the layer''s mobility DOC leaches at 0.004 and the other pools as before.');
unit = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
par_src = fileread(fullfile(code_dir,'Biogeochemistry_Parameter.m'));
assert(contains(par_src,'BiogeoPar.aDOC = 0.04;'), ...
    'The parameter script must keep aDOC at 0.04.');
k_if = strfind(unit,'if leach_layer');
k_set = strfind(unit,'BiogeoPar.aDOC = 0.004;');
k_leak = strfind(unit,'Biogeo_Leakage(Ptm1,Lk,V,BiogeoPar)');
assert(isscalar(k_if) && isscalar(k_set) && isscalar(k_leak) && k_if < k_set && k_set < k_leak && ...
    contains(unit,'biogeo_overrides,leach_layer)') && contains(unit,'leach_layer = false;'), ...
    'BIOGEO_UNIT must lower aDOC only under leach_layer, before the leaching, and default to the column.');
assert(contains(one,'biogeo_exchange_scheme,[],[],leach_layer);') && ...
    contains(two,'biogeo_overrides{c},leach_layer);') && ...
    contains(spin,'[], [], [], leach_layer);') && ...
    contains(spin,"leach_layer = strcmp(Leaching_Water_Basis(MOD_PARAM_input),'layer');"), ...
    'Both drivers and the spin-up''s Phase 2 must hand BIOGEO_UNIT the leaching water basis.');

fprintf('PASS: leaching_water selects the column''s or the layer''s drainage and water.\n');
end

function assert_throws(fun,identifier)
threw = false;
try
    fun();
catch ME
    threw = strcmp(ME.identifier,identifier);
end
assert(threw,'Expected error %s.',identifier);
end
