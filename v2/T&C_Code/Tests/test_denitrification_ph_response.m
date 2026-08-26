function test_denitrification_ph_response
% Gross denitrification and N2:N2O partition use published pH functions.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

pH = [3.5 4 5.5 7 8];
[gross_factor,n2_to_n2o] = Biogeo_Denitrification_PH_Response(pH);
expected_gross = 0.0016*exp(1.006*pH);
expected_ratio = [0.001 0.001 0.501 1 1];

assert(max(abs(gross_factor-expected_gross)) < 1e-12, ...
    'Gross denitrification pH response does not match Val Martin et al. (2023).');
assert(max(abs(n2_to_n2o-expected_ratio)) < 1e-12, ...
    'N2:N2O pH response does not match Val Martin et al. (2023).');
assert(all(diff(gross_factor) > 0) && all(diff(n2_to_n2o) >= 0), ...
    'Higher pH must not reduce gross denitrification or the N2:N2O ratio.');
assert(all(gross_factor > 0), ...
    'The gross denitrification pH factor must remain positive.');

assert_throws(@() Biogeo_Denitrification_PH_Response(0), ...
    'Biogeo_Denitrification_PH_Response:InvalidPH');
assert_throws(@() Biogeo_Denitrification_PH_Response(NaN), ...
    'Biogeo_Denitrification_PH_Response:InvalidPH');

fprintf('PASS: denitrification pH functions match Val Martin et al. (2023).\n');
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
