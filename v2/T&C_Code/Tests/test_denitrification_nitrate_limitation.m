function test_denitrification_nitrate_limitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Denitrification must retain the nitrate half-saturation in the active
%% rate law. Otherwise the per-unit-NO3 loss is independent of inventory
%% and wet maize events overdraw the event-scale N2O flux.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));
fprintf('\n=== test_denitrification_nitrate_limitation ===\n');

F = load(fullfile(here,'Fixtures','no3_immobilization_20190906.mat'),'fixture');
args = F.fixture.args;
args{end+1} = 0.8;
args{50} = Biogeochemistry_Parameter(1);
for k = 12:47
    args{k} = 0;
end
args{5}(:) = 0;
args{52}(:) = 0;
args{6} = 20;
args{7} = 20;
args{9} = 5;
args{10} = 1;
args{11} = 0.5;
args{55} = 1;

soil_mass = args{3}*args{4}*1000;
Kn = args{50}.Kn;
NO3_lo = 0.1*Kn;
NO3_hi = Kn;
DOC = 10*args{50}.Kc;

rate_lo = denitrification_rate(args,DOC,NO3_lo);
rate_hi = denitrification_rate(args,DOC,NO3_hi);
observed_ratio = (rate_lo/NO3_lo)/(rate_hi/NO3_hi);
expected_ratio = (NO3_lo/(NO3_lo+Kn))/(NO3_hi/(NO3_hi+Kn));
fprintf('soil_mass %.6g, Kn %.6g, rates %.6g %.6g, ratio %.6g\n', ...
    soil_mass,Kn,rate_lo,rate_hi,observed_ratio);

assert(rate_lo > 0 && rate_hi > 0, ...
    'The nitrate-limitation fixture must produce finite positive denitrification.');
assert(abs(observed_ratio-expected_ratio) < 1e-8, ...
    ['Denitrification does not follow NO3/(NO3+Kn): observed normalized ' ...
    'rate ratio %.6g, expected %.6g.'],observed_ratio,expected_ratio);
fprintf('PASS: nitrate-normalized rate ratio %.6g matches %.6g\n', ...
    observed_ratio,expected_ratio);

bio_source = fileread(fullfile(code,'BIOGEOCHEMISTRY_DYNAMIC3.m'));
assert(contains(bio_source,'B(31)= B(31)- NH4_Uptake') && ...
       contains(bio_source,'B(32)=B(32)-NO3_Uptake'), ...
    'Working mineral-N pools must reserve plant uptake once.');
assert(contains(bio_source,'(1-f_org_lea)*Ip_lea_nit - LEAK_NH4'), ...
    'NH4 reaction availability must use the post-uptake working pool.');
assert(contains(bio_source,'N2Oflx_nit - LEAK_NO3'), ...
    'NO3 reaction availability must use the post-uptake working pool.');
fix_source = fileread(fullfile(code,'Biogeo_Bio_fixation.m'));
assert(contains(fix_source,'REXMY(3)/Costfix'), ...
    'Biological N fixation must use only carbon allocated to fixation.');
fprintf('PASS: mineral-N reactions reserve uptake and fixation uses its C allocation\n');
end

function rate = denitrification_rate(args,DOC,NO3)
args{2}(1:34) = max(args{2}(1:34),1e-3);
args{2}(31:32) = 0;
args{2}(12) = DOC;
args{2}(13) = 0;
args{2}(32) = NO3;
[~,~,~,~,~,~,N2flx,~,~,~,~,~,~,~,N2Oflx] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(args{:});
rate = N2flx+N2Oflx;
end
