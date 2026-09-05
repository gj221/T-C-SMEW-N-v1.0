function test_plant_reserve_fraction()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_plant_reserve_fraction
%%
%% The reserve-fill fraction and the uptake suppression used by
%% Plant_Exports:
%%   1  reserve within capacity -> the ratio
%%   2  reserve above capacity  -> saturated at 1
%%   3  a NEGATIVE reserve is a broken state under the physical-inventory
%%      design (Plant_Reserve_Growth_Limit never overdraws it; VEG_DYN_RES
%%      clears round-off), reported with its value
%%   4  zero capacity keeps the saturated limit: 1 if there is reserve,
%%      0 otherwise
%%   5  non-finite reserve, non-finite or negative capacity are errors
%%      that report their values
%%   6  the suppression is clipped to [0,1] and equals the reserve
%%      fraction when its denominator vanishes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
fprintf('\n=== test_plant_reserve_fraction ===\n');
n_pass = 0; n_fail = 0;
    function check(ok,msg)
        if ok; n_pass = n_pass+1; fprintf('  PASS  %s\n',msg);
        else;  n_fail = n_fail+1; fprintf('  FAIL  %s\n',msg); end
    end
check(abs(Plant_Reserve_Fraction(0.25,1,'N') - 0.25) < 1e-15,'reserve within capacity gives the ratio');
check(Plant_Reserve_Fraction(3,1,'N') == 1,'reserve above capacity saturates at 1');
ok = false;
try; Plant_Reserve_Fraction(-3.7e-11,0.005,'N');
catch ME; ok = strcmp(ME.identifier,'Plant_Reserve_Fraction:InvalidState') && contains(ME.message,'-3.7e-11'); end
check(ok,'a negative reserve is an error that reports its value');
check(Plant_Reserve_Fraction(0.2,0,'N') == 1 && Plant_Reserve_Fraction(0,0,'N') == 0, ...
    'zero capacity keeps the saturated limit');
bad = {{NaN,1},{1,NaN},{1,-1},{Inf,1},{-1,1}};
ok = true;
for k = 1:numel(bad)
    try
        Plant_Reserve_Fraction(bad{k}{1},bad{k}{2},'K'); ok = false;
    catch ME
        ok = ok && strcmp(ME.identifier,'Plant_Reserve_Fraction:InvalidState') && ...
            contains(ME.message,'capacity');
    end
end
check(ok,'non-finite reserve and non-finite or negative capacity are errors reporting their values');
check(Plant_Nutrient_Suppression(5,1,1,'N') == 1 && ...
    abs(Plant_Nutrient_Suppression(0.9,1,1,'N') - (0.9-0.8)/(0.6+0.2)) < 1e-15, ...
    'suppression is the clipped TC-SMEW ramp');
check(Plant_Nutrient_Suppression(0.5,0,0,'N') == 1,'a vanishing denominator returns the reserve fraction');
fprintf('  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0; error('test_plant_reserve_fraction: %d check(s) failed',n_fail); end
end
