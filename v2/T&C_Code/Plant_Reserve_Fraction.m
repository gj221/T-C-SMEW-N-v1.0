function fraction = Plant_Reserve_Fraction(reserve,capacity,element_name)
% Retain the saturated TC-SMEW limit at the zero-capacity boundary.

if ~isfinite(reserve) || ~isfinite(capacity) || reserve < 0 || capacity < 0
    error('Plant_Reserve_Fraction:InvalidState', ...
        '%s reserve and storage capacity must be finite and nonnegative.',element_name);
end

tol = 100*eps(max([1,reserve,capacity]));
if capacity <= tol
    fraction = double(reserve > tol);
    return
end

fraction = min(max(reserve/capacity,0),1);
end
