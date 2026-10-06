function cols = Read_Parameter_Columns(file_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A site's parameter-row columns from a one-row CSV, as the struct
%% Run_Field_Cell_Spinup takes ('MODParamColumns', 'ERWParamColumns').
%% A number stays a number, NaN included; anything else -- a name, a date,
%% a vector such as Zs -- is text, as the workbook rows hold it.  An empty
%% file name gives no columns; an empty value is refused.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cols = struct();
if isempty(file_name); return; end
file_name = char(file_name);
if ~isfile(file_name)
    error('Read_Parameter_Columns:MissingFile','%s not found.',file_name);
end
opts = detectImportOptions(file_name,'Delimiter',',','VariableNamingRule','preserve');
opts = setvartype(opts,'string');
T = readtable(file_name,opts);
if height(T) ~= 1
    error('Read_Parameter_Columns:BadFile','%s must hold one row of column values.',file_name);
end
for k = 1:width(T)
    v = strtrim(T{1,k});
    if ismissing(v) || strlength(v) == 0
        error('Read_Parameter_Columns:BadFile','%s: column %s has no value.', ...
            file_name,T.Properties.VariableNames{k});
    end
    x = str2double(v);
    if ~isnan(x) || strcmpi(v,"NaN"); v = x; end
    cols.(T.Properties.VariableNames{k}) = v;
end
end
