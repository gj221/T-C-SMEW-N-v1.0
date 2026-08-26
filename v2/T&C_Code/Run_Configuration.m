%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfg = Run_Configuration(run_chemistry, run_domain, buffer_planting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3 || isempty(run_chemistry) || isempty(run_domain) || isempty(buffer_planting)
    error('Run_Configuration:MissingInput', ...
        'Run chemistry, domain, and buffer planting must be provided.');
end

run_chemistry = lower(strrep(strtrim(char(run_chemistry)), '-', '_'));
if strcmp(run_chemistry, 'ph_only')
    run_chemistry = 'pH_only';
end

run_domain      = lower(strtrim(char(run_domain)));
buffer_planting = lower(strtrim(char(buffer_planting)));

run_domain      = strrep(run_domain,'-','_');
buffer_planting = strrep(buffer_planting,'-','_');

switch run_chemistry
    case 'off'
        cfg.OPT_pH = 0;
        cfg.OPT_ERW = 0;
        cfg.OPT_ROCK = 0;
    case 'pH_only'
        cfg.OPT_pH = 1;
        cfg.OPT_ERW = 1;
        cfg.OPT_ROCK = 0;
    case 'erw'
        cfg.OPT_pH = 1;
        cfg.OPT_ERW = 1;
        cfg.OPT_ROCK = 1;
    otherwise
        error('Run_Configuration:BadRunMode', ...
            'RUN_CHEMISTRY must be off, pH_only, or erw.');
end

switch run_domain
    case 'single'
        cfg.OPT_GRID = 0;
        cfg.domain_label = 'single';
    case {'field_buffer','grid'}
        cfg.OPT_GRID = 1;
        cfg.domain_label = 'field_buffer';
    otherwise
        error('Run_Configuration:BadDomain', ...
            'RUN_DOMAIN must be single or field_buffer.');
end

if cfg.OPT_GRID == 0
    cfg.OPT_BUFFER_PERENNIAL = 0;
    cfg.buffer_label = 'none';
else
    switch buffer_planting
        case {'crop','annual'}
            cfg.OPT_BUFFER_PERENNIAL = 0;
            cfg.buffer_label = 'crop';
        case {'native_perennial','perennial'}
            cfg.OPT_BUFFER_PERENNIAL = 1;
            cfg.buffer_label = 'native_perennial';
        otherwise
            error('Run_Configuration:BadBufferPlanting', ...
                'BUFFER_PLANTING must be crop or native_perennial.');
    end
end

cfg.run_chemistry = run_chemistry;
end
