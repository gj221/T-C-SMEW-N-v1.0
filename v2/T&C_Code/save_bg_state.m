function save_bg_state(bg_file, bg_state)
%%% Written beside the target and renamed into place, so a concurrent reader
%%% (another array task, a resubmission) sees the old file or the complete new
%%% one, never a partial state.
bg_dir = fileparts(bg_file);
if isempty(bg_dir)
    bg_dir = pwd;
end
if ~exist(bg_dir,'dir')
    mkdir(bg_dir);
end
tmp_file = [tempname(bg_dir) '.mat'];
try
    save(tmp_file, '-struct', 'bg_state', '-v7');
    movefile(tmp_file, bg_file, 'f');
catch ME
    if isfile(tmp_file)
        delete(tmp_file);
    end
    rethrow(ME);
end
end
