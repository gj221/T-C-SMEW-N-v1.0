%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bg_file = BG_State_File(Input_Directory, id_location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2 || isempty(id_location)
    error('BG_State_File:MissingSiteId', 'id_location is required for the BG state file.');
end
site_id = strtrim(char(string(id_location)));
if isempty(site_id)
    error('BG_State_File:MissingSiteId', 'id_location is required for the BG state file.');
end
site_id = regexprep(site_id, '\s+', '_');
%% Either path separator becomes '_'; char(92) keeps the source free of
%% backslash literals (test_bouchet_linux_portability).
site_id = regexprep(strrep(site_id, char(92), '/'), '/+', '_');
bg_file = fullfile(Input_Directory, ['Data_' site_id '_BG.mat']);
end
