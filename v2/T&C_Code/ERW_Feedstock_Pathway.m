function [kinetic_mass_fraction, direct_carbonate_mass_fraction, rock_f_kinetic] = ...
    ERW_Feedstock_Pathway(mineral, rock_f_in)

if iscell(mineral) && isscalar(mineral); mineral = mineral{1}; end
if iscell(rock_f_in) && isscalar(rock_f_in); rock_f_in = rock_f_in{1}; end
if ((ischar(mineral) && isrow(mineral)) || ...
        (isstring(mineral) && isscalar(mineral))) && ...
        startsWith(strtrim(string(mineral)),"[")
    names = strtrim(string(jsondecode(char(mineral))));
else
    names = strtrim(string(cellstr(mineral)));
end
if ischar(rock_f_in) || isstring(rock_f_in)
    fractions = double(jsondecode(char(rock_f_in)));
else
    fractions = rock_f_in;
end
rock_f_kinetic = fractions;
if isempty(fractions) || (isnumeric(fractions) && all(isnan(fractions(:))))
    kinetic_mass_fraction = 0;
    direct_carbonate_mass_fraction = 0;
    return
end

fractions = fractions(:).';
if numel(names) ~= numel(fractions) || any(~isfinite(fractions)) || ...
        any(fractions < 0) || sum(fractions) > 1 + 1e-10
    error('ERW:invalidFeedstockMineralogy', ...
        'Mineral names and finite nonnegative mass fractions must match and sum to at most one.');
end

unsupported = ismember(lower(names),["aragonite","dolomite","magnesite"]);
if any(unsupported)
    error('ERW:unsupportedCarbonate', ...
        'Only calcite has a carbonate mass and stoichiometry pathway in this model.');
end

kinetic_mass_fraction = 1;
direct_carbonate_mass_fraction = 0;
rock_f_kinetic = fractions;

end
