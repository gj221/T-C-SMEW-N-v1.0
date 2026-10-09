function test_dead_crown_reserve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A crown whose living tissue has died away and left a reserve behind.
%% The case is Hamden's open-grass crown at the end of its class seed: the
%% cell has no open grass (Ccrown(2) = 0), so the grass has no share of the
%% soil water, starved within four years, and its leaf, root and
%% carbohydrate carbon then decayed to ~1e-174, 1e-259 and 1e-30 gC/m2
%% while 0.47 gN/m2 stayed in its N reserve.  The relative N concentration
%% of that crown diverges and Nutrients_Available's e_relN > 100 guard
%% trips every day.  The guard reports and nothing more: every output is
%% the one the equations give, the crown takes nothing up (SupN = 1), and
%% the message comes at the 1st, 10th, 100th, ... occurrence instead of
%% once a day for the rest of the run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);
St = Veg_Stoichiometric_Parameter(23);     % the grass of PFT_Hamden_Deciduous_Grass
B = [1.18904e-174 0 8.51393e-260 1.59953e-30 1.35374e-321 0 3.33256 0];
res = [0.470607 0.0336148 0.0770398 0.176736 0.0565162 0.0951838];

clear Nutrients_Available                  % start the occurrence count afresh
restore = onCleanup(@() clear('Nutrients_Available')); %#ok<NASGU>
out = cell(1,16);
txt = evalc(['[out{:}] = Nutrients_Available(B,1,St,res(1),res(2),res(3),' ...
    'res(4),res(5),res(6),1);']);
assert(all(cellfun(@(x) isscalar(x) && isfinite(x),out)), ...
    'test_dead_crown_reserve:Unassigned', ...
    'Every output must be assigned and finite when the guard trips.');
[FNC,e_relN,rNc,~,~,~,~,~,~,rNcR,Navail] = out{1:11};

%% The equations, unchanged by the guard.
TNITns = B(1)/St.Nl + B(3)/St.Nr + B(4)/St.Ns + B(5)/St.Nf;
Nsto = St.FiS*B(4)*(1/St.Nl - 1/St.Ns);
rNc_eq = 1 + (res(1) - Nsto)/TNITns;
assert(e_relN > 100, 'test_dead_crown_reserve:NotTripped', ...
    'The dead crown must trip the e_relN > 100 guard (e_relN = %.3g).',e_relN);
assert(abs(rNc/rNc_eq - 1) < 1e-12 && abs(e_relN/(0.5*(rNc_eq-1)+1) - 1) < 1e-12 && ...
    rNcR == 1.85 && FNC == 1 && abs(Navail - (res(1) + 0.35*TNITns)) < 1e-15, ...
    'test_dead_crown_reserve:Changed', ...
    ['The guard must leave the outputs as the equations give them ' ...
     '(rNc %.6g vs %.6g, e_relN %.6g, rNcR %.3g, FNC %.3g).'],rNc,rNc_eq,e_relN,rNcR,FNC);
%% Its reserve is full, so the crown takes nothing up (BIOGEO_UNIT scales
%% each crown's uptake by 1 - SupN).
SupN = Plant_Nutrient_Suppression(res(1),Nsto,TNITns,'N');
assert(SupN == 1, 'test_dead_crown_reserve:Uptake', ...
    'A dead crown with a reserve must take no N up (SupN = %.3g).',SupN);

%% The message: the 1st and the 10th of 25 occurrences.
for k = 2:25
    txt = [txt evalc(['[out{:}] = Nutrients_Available(B,1,St,res(1),res(2),res(3),' ...
        'res(4),res(5),res(6),1);'])]; %#ok<AGROW>
end
n_msg = numel(strfind(txt,'unrealistic Photosynthesis N efficiency'));
assert(n_msg == 2 && contains(txt,'occurrence 1 ') && contains(txt,'occurrence 10 '), ...
    'test_dead_crown_reserve:Messages', ...
    '25 occurrences must print the 1st and the 10th only; printed %d:\n%s',n_msg,txt);

%% A living crown is left alone: no message.
B_live = [40 0 120 80 0 0 2 0];
TN_live = B_live(1)/St.Nl + B_live(3)/St.Nr + B_live(4)/St.Ns + B_live(5)/St.Nf;
Nsto_live = St.FiS*B_live(4)*(1/St.Nl - 1/St.Ns);
txt_live = evalc(['[out{:}] = Nutrients_Available(B_live,1,St,Nsto_live+0.1*TN_live,' ...
    '0,0,0,0,0,1);']);
assert(~contains(txt_live,'unrealistic') && abs(out{3} - 1.1) < 1e-12, ...
    'test_dead_crown_reserve:Living', ...
    'A living crown 10%% above its N storage must read rNc = 1.1 without a message (rNc %.6g).',out{3});

fprintf(['PASS: a dead crown''s reserve trips the e_relN guard (e_relN %.3g) with its ' ...
    'outputs unchanged and SupN = 1; 25 occurrences print 2 messages.\n'],e_relN);
end
