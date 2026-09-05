function fraction = Plant_Reserve_Fraction(reserve,capacity,element_name)
% Retain the saturated TC-SMEW limit at the zero-capacity boundary.
% The reserve is a physical inventory (Plant_Reserve_Growth_Limit scales
% growth so it is never overdrawn), so a negative value is a broken
% state; VEG_DYN_RES clears the round-off of a balance that lands on
% zero before it reaches here.

if ~isfinite(reserve) || ~isfinite(capacity) || reserve < 0 || capacity < 0
    error('Plant_Reserve_Fraction:InvalidState', ...
        ['%s reserve and storage capacity must be finite and nonnegative ' ...
         '(reserve %.6g, capacity %.6g).'],element_name,reserve,capacity);
end

tol = 100*eps(max([1,reserve,capacity]));
if capacity <= tol
    fraction = double(reserve > tol);
    return
end

fraction = min(max(reserve/capacity,0),1);
end
