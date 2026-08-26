%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Date_out = Shift_Donor_Dates_To_Year(Date_in, target_year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(Date_in);
Date_vec = Date_in(:);
[~,MO,DA,HO,MI,SE] = datevec(Date_vec);
Date_out = datenum(target_year*ones(size(Date_vec)),MO,DA,HO,MI,SE);
Date_out = reshape(Date_out, sz);
end
