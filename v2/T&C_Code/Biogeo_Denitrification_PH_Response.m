function [gross_factor,n2_to_n2o] = Biogeo_Denitrification_PH_Response(PH)
% Soil-pH factors used by the Val Martin et al. (2023) N2O scheme.

if ~isreal(PH) || any(~isfinite(PH(:))) || any(PH(:) <= 0 | PH(:) >= 14)
    error('Biogeo_Denitrification_PH_Response:InvalidPH', ...
        'Soil pH must be finite and satisfy 0 < pH < 14.');
end

gross_factor = 0.0016*exp(1.006*PH);
n2_to_n2o = ones(size(PH));
n2_to_n2o(PH <= 4) = 0.001;
middle = PH > 4 & PH < 7;
n2_to_n2o(middle) = 0.001 + (PH(middle)-4)/3;
end
