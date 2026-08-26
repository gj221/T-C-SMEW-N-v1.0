function save_bg_state(bg_file, bg_state)
bg_dir = fileparts(bg_file);
if ~exist(bg_dir,'dir')
    mkdir(bg_dir);
end
save(bg_file, '-struct', 'bg_state', '-v7');
end
