function bg_state = load_bg_state(bg_file)
if ~isfile(bg_file)
    error('load_bg_state:MissingBGState', ...
        'BG state not found: %s. Run TC_SMEW_N_run with RUN_MODE = ''spinup'' first.', bg_file);
end
bg_state = load(bg_file);
required = {'B_out','NavlI_out','NupI_L_out','O_out','RexmyI_out'};
for ii = 1:numel(required)
    if ~isfield(bg_state, required{ii})
        error('load_bg_state:InvalidBGState', 'BG state %s is missing %s.', bg_file, required{ii});
    end
end
if any(~isfinite(bg_state.B_out(:))) || any(~isfinite(bg_state.O_out(:)))
    error('load_bg_state:InvalidBGState', 'BG state %s contains non-finite B_out or O_out values.', bg_file);
end
end
