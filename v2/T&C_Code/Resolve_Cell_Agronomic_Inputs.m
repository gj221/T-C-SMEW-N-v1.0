function resolved = Resolve_Cell_Agronomic_Inputs(topo,input_table)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolve measured soil properties or their interpolated cell baselines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(topo) || ~isscalar(topo) || ~isfield(topo,'N_cells') || ...
        ~isfield(topo,'cell_id')
    error('Resolve_Cell_Agronomic_Inputs:InvalidTopology', ...
        'The topology must contain N_cells and cell_id.');
end
if ~istable(input_table)
    error('Resolve_Cell_Agronomic_Inputs:InvalidInput', ...
        'Cell agronomic inputs must be supplied as a table.');
end
required = {'cell_id','baseline_bulk_pH','measured_bulk_pH', ...
    'use_measured_bulk_pH','baseline_sand_fraction', ...
    'baseline_clay_fraction','measured_sand_fraction', ...
    'measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC'};
missing = required(~ismember(required,input_table.Properties.VariableNames));
if ~isempty(missing)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'Cell agronomic input is missing: %s.',strjoin(missing,', '));
end

cell_id = input_table.cell_id;
valid_id = isnumeric(cell_id) && isreal(cell_id) && iscolumn(cell_id) && ...
    all(isfinite(cell_id)) && all(cell_id == fix(cell_id)) && all(cell_id >= 1);
if ~valid_id || numel(unique(cell_id)) ~= height(input_table)
    error('Resolve_Cell_Agronomic_Inputs:InvalidCellId', ...
        'cell_id must contain unique positive integers.');
end
topology_id = topo.cell_id(:);
if topo.N_cells ~= numel(topology_id) || height(input_table) ~= topo.N_cells || ...
        ~isequal(sort(cell_id),sort(topology_id))
    error('Resolve_Cell_Agronomic_Inputs:CellSetMismatch', ...
        'The agronomic table must contain one row for every topology cell.');
end
[~,order] = ismember(topology_id,cell_id);
T = input_table(order,:);

[cell_x,cell_y] = topology_centres(topo);
has_x = ismember('cell_x_projected_m',T.Properties.VariableNames);
has_y = ismember('cell_y_projected_m',T.Properties.VariableNames);
if xor(has_x,has_y)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'Projected cell coordinates must supply both x and y columns.');
end
if has_x
    x = T.cell_x_projected_m;
    y = T.cell_y_projected_m;
    if any(~isfinite(x)) || any(~isfinite(y)) || ...
            any(abs(x-cell_x) > 1e-6) || any(abs(y-cell_y) > 1e-6)
        error('Resolve_Cell_Agronomic_Inputs:CellCoordinateMismatch', ...
            'Projected cell coordinates do not match the selected grid topology.');
    end
end

use_pH = logical_flag(T.use_measured_bulk_pH,'use_measured_bulk_pH');
use_texture = logical_flag(T.use_measured_texture,'use_measured_texture');
use_OM = logical_flag(T.use_measured_OM,'use_measured_OM');
use_CEC = logical_flag(T.use_measured_CEC,'use_measured_CEC');

baseline_names = {'baseline_bulk_pH','baseline_sand_fraction', ...
    'baseline_clay_fraction','baseline_OM_fraction','baseline_CEC_cmolc_kg'};
for k = 1:numel(baseline_names)
    values = T.(baseline_names{k});
    if ~isnumeric(values) || ~isreal(values) || ~iscolumn(values) || ...
            any(~isfinite(values))
        error('Resolve_Cell_Agronomic_Inputs:MissingBaseline', ...
            '%s must contain one finite value for every simulated cell.',baseline_names{k});
    end
end
require_measurement(T.measured_bulk_pH,use_pH,'measured_bulk_pH');
require_measurement(T.measured_sand_fraction,use_texture,'measured_sand_fraction');
require_measurement(T.measured_clay_fraction,use_texture,'measured_clay_fraction');
require_measurement(T.measured_OM_fraction,use_OM,'measured_OM_fraction');
require_measurement(T.measured_CEC_cmolc_kg,use_CEC,'measured_CEC_cmolc_kg');

bulk_pH = choose(T.baseline_bulk_pH,T.measured_bulk_pH,use_pH);
sand = choose(T.baseline_sand_fraction,T.measured_sand_fraction,use_texture);
clay = choose(T.baseline_clay_fraction,T.measured_clay_fraction,use_texture);
OM = choose(T.baseline_OM_fraction,T.measured_OM_fraction,use_OM);
CEC = choose(T.baseline_CEC_cmolc_kg,T.measured_CEC_cmolc_kg,use_CEC);
silt = 1-sand-clay-OM;

if any(bulk_pH <= 0 | bulk_pH > 14)
    error('Resolve_Cell_Agronomic_Inputs:InvalidPH', ...
        'Resolved bulk pH must be greater than zero and no greater than 14.');
end
if any(sand < 0 | sand >= 1 | clay < 0 | clay >= 1 | OM < 0 | OM >= 1) || ...
        any(silt < -1e-12)
    error('Resolve_Cell_Agronomic_Inputs:TextureClosure', ...
        'Sand, clay, organic matter and derived silt fractions must be nonnegative and sum to one.');
end
if any(CEC <= 0)
    error('Resolve_Cell_Agronomic_Inputs:InvalidCEC', ...
        'Resolved CEC must be positive in cmol charge kg-1.');
end
silt = max(silt,0);

pH_source = source_label(use_pH);
texture_source = source_label(use_texture);
OM_source = source_label(use_OM);
CEC_source = source_label(use_CEC);
input_identity = strings(topo.N_cells,1);
for c = 1:topo.N_cells
    input_identity(c) = sprintf([ ...
        'cell=%d;x=%.17g;y=%.17g;pH=%.17g:%s;sand=%.17g;clay=%.17g:' ...
        '%s;OM=%.17g:%s;CEC=%.17g:%s'],topology_id(c),cell_x(c),cell_y(c), ...
        bulk_pH(c),pH_source(c),sand(c),clay(c),texture_source(c), ...
        OM(c),OM_source(c),CEC(c),CEC_source(c));
end

resolved = table(topology_id,cell_x,cell_y,bulk_pH,sand,clay,silt,OM,CEC, ...
    pH_source,texture_source,OM_source,CEC_source,input_identity, ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'bulk_pH','sand_fraction','clay_fraction','silt_fraction','OM_fraction', ...
    'CEC_cmolc_kg','bulk_pH_source','texture_source','OM_source','CEC_source', ...
    'input_identity'});
end

function [x,y] = topology_centres(topo)
if isfield(topo,'cell_xc_projected_m') && isfield(topo,'cell_yc_projected_m')
    x = topo.cell_xc_projected_m(:);
    y = topo.cell_yc_projected_m(:);
else
    x = nan(topo.N_cells,1);
    y = nan(topo.N_cells,1);
end
if numel(x) ~= topo.N_cells || numel(y) ~= topo.N_cells
    error('Resolve_Cell_Agronomic_Inputs:InvalidTopology', ...
        'Projected topology centres must contain one coordinate per cell.');
end
end

function flag = logical_flag(values,name)
valid = (islogical(values) || isnumeric(values)) && isreal(values) && ...
    iscolumn(values) && all(isfinite(double(values))) && ...
    all(double(values) == 0 | double(values) == 1);
if ~valid
    error('Resolve_Cell_Agronomic_Inputs:InvalidToggle', ...
        '%s must contain only logical or numeric zero/one values.',name);
end
flag = logical(values);
end

function require_measurement(values,use_measurement,name)
valid_type = isnumeric(values) && isreal(values) && iscolumn(values);
if ~valid_type || any(~isfinite(values(use_measurement)))
    error('Resolve_Cell_Agronomic_Inputs:MissingMeasurement', ...
        '%s must be finite wherever its measurement toggle is enabled.',name);
end
end

function value = choose(baseline,measured,use_measurement)
value = baseline;
value(use_measurement) = measured(use_measurement);
end

function label = source_label(use_measurement)
label = repmat("baseline",numel(use_measurement),1);
label(use_measurement) = "measured";
end
