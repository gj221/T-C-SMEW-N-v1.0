function checkpoint_file = spinup_phase_checkpoint_filename(bg_file, phase_label)
%SPINUP_PHASE_CHECKPOINT_FILENAME  Per-phase spin-up checkpoint path.
%   Derives <bg>_checkpoint_<phase_label>.mat next to bg_file.
[bg_dir,bg_name,bg_ext] = fileparts(bg_file);
phase_label = regexprep(char(phase_label), '[^A-Za-z0-9]+', '_');
checkpoint_file = fullfile(bg_dir, [bg_name '_checkpoint_' phase_label bg_ext]);
end
