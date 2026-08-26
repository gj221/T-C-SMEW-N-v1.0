%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Run_Descriptor  --  compact domain/amendment/vegetation label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function descriptor = Run_Descriptor(run_config, ERW_region_coverage, ERW_PARAM_input)

if nargin < 2 || isempty(ERW_region_coverage)
    ERW_region_coverage = struct('apply_main', true, 'apply_buffer', true);
end
if nargin < 3
    ERW_PARAM_input = [];
end

if isfield(run_config,'OPT_GRID') && run_config.OPT_GRID == 1
    domain_tag = 'domain_field_buffer';
else
    domain_tag = 'domain_single';
end

apply_main = true;
apply_buffer = true;
if isstruct(ERW_region_coverage)
    if isfield(ERW_region_coverage,'apply_main')
        apply_main = logical(ERW_region_coverage.apply_main);
    end
    if isfield(ERW_region_coverage,'apply_buffer')
        apply_buffer = logical(ERW_region_coverage.apply_buffer);
    end
end

rock_enabled = false;
if isfield(run_config,'OPT_ROCK')
    rock_enabled = run_config.OPT_ROCK == 1;
elseif isfield(run_config,'OPT_ERW')
    rock_enabled = run_config.OPT_ERW == 1;
end

if rock_enabled && is_calcite_run(ERW_PARAM_input)
    calcite_tag = ['calcite_' placement_tag(run_config, apply_main, apply_buffer)];
else
    calcite_tag = 'calcite_none';
end

native_tag = 'native_none';
if isfield(run_config,'OPT_GRID') && run_config.OPT_GRID == 1 && ...
        isfield(run_config,'OPT_BUFFER_PERENNIAL') && run_config.OPT_BUFFER_PERENNIAL == 1
    native_tag = 'native_buffer';
end

descriptor = string([domain_tag '__' calcite_tag '__' native_tag]);
end

function tf = is_calcite_run(ERW_PARAM_input)
tf = false;
if isempty(ERW_PARAM_input); return; end
M_rock_in = table_or_struct_value(ERW_PARAM_input,'M_rock_in',0);
if ~(isnumeric(M_rock_in) && M_rock_in > 0); return; end
mineral = lower(string(table_or_struct_value(ERW_PARAM_input,'mineral',"")));
description = lower(string(table_or_struct_value(ERW_PARAM_input,'ExperimentalDescription',"")));
tf = contains(mineral,'calcite') || contains(description,'calcite');
end

function tag = placement_tag(run_config, apply_main, apply_buffer)
if isfield(run_config,'OPT_GRID') && run_config.OPT_GRID == 1
    if apply_main && apply_buffer
        tag = 'main_buffer';
    elseif apply_main
        tag = 'main';
    elseif apply_buffer
        tag = 'buffer';
    else
        tag = 'none';
    end
else
    if apply_main
        tag = 'wholefield';
    else
        tag = 'none';
    end
end
end

function val = table_or_struct_value(x, name, default)
val = default;
if istable(x)
    if ~ismember(name,x.Properties.VariableNames); return; end
    raw = x.(name);
    if iscell(raw); raw = raw{1}; end
    if isstring(raw) && ~isempty(raw); val = raw(1); return; end
    if ischar(raw); val = raw; return; end
    if isnumeric(raw) && ~isempty(raw); val = raw(1); return; end
elseif isstruct(x) && isfield(x,name)
    raw = x.(name);
    if iscell(raw); raw = raw{1}; end
    if isstring(raw) && ~isempty(raw); val = raw(1); return; end
    if ischar(raw); val = raw; return; end
    if isnumeric(raw) && ~isempty(raw); val = raw(1); return; end
end
end
