function states = Load_Field_Spinup_States(state_file,topo,input_identities,start_dn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load, order and validate one automatic spin-up handoff per grid cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(ischar(state_file) || (isstring(state_file) && isscalar(state_file))) || ...
        ~isfile(state_file)
    error('Load_Field_Spinup_States:MissingFile', ...
        'Cell spin-up state file does not exist: %s.',char(string(state_file)));
end
if ~isstruct(topo) || ~isfield(topo,'cell_id') || ~isfield(topo,'N_cells')
    error('Load_Field_Spinup_States:InvalidTopology', ...
        'The topology must contain N_cells and cell_id.');
end
input_identities = string(input_identities(:));
if numel(input_identities) ~= topo.N_cells
    error('Load_Field_Spinup_States:IdentityCount', ...
        'One agronomic identity is required per topology cell.');
end
raw = load(state_file,'CellBGStates');
if ~isfield(raw,'CellBGStates') || ~isstruct(raw.CellBGStates) || ...
        numel(raw.CellBGStates) ~= topo.N_cells
    error('Load_Field_Spinup_States:InvalidFile', ...
        'The state file must contain CellBGStates with one struct per cell.');
end
states_in = raw.CellBGStates(:);
if ~all(isfield(states_in,{'cell_id','cell_agronomic_identity'}))
    error('Load_Field_Spinup_States:InvalidFile', ...
        'Every spin-up state requires cell_id and cell_agronomic_identity.');
end
state_ids = [states_in.cell_id].';
valid_ids = isnumeric(state_ids) && all(isfinite(state_ids)) && ...
    all(state_ids == fix(state_ids)) && numel(unique(state_ids)) == topo.N_cells;
if ~valid_ids || ~isequal(sort(state_ids),sort(topo.cell_id(:)))
    error('Load_Field_Spinup_States:CellSetMismatch', ...
        'Spin-up state cell ids do not match the selected grid topology.');
end
[~,order] = ismember(topo.cell_id(:),state_ids);
states = states_in(order).';
required = {'B_out','O_out','RexmyI_out','NavlI_out','NupI_L_out'};
for c = 1:topo.N_cells
    if ~strcmp(char(string(states(c).cell_agronomic_identity)),input_identities(c))
        error('Load_Field_Spinup_States:InputIdentityMismatch', ...
            'Cell %d spin-up state was produced from different agronomic inputs.', ...
            topo.cell_id(c));
    end
    missing = required(~isfield(states(c),required));
    if ~isempty(missing)
        error('Load_Field_Spinup_States:InvalidFile', ...
            'Cell %d spin-up state is missing: %s.',topo.cell_id(c),strjoin(missing,', '));
    end
    for k = 1:numel(required)
        values = states(c).(required{k});
        if ~isnumeric(values) || ~isreal(values) || any(~isfinite(values(:)))
            error('Load_Field_Spinup_States:InvalidFile', ...
                'Cell %d spin-up field %s must be finite and numeric.', ...
                topo.cell_id(c),required{k});
        end
    end
    %% A state seeds the transient that starts at the instant it handed off
    %% at, and no other: from an earlier one the run would repeat forcing the
    %% pools already integrated, from a later one it would skip some.
    if nargin >= 4 && ~isempty(start_dn)
        if ~isfield(states(c),'spinup_handoff_dn') || ...
                ~(isnumeric(states(c).spinup_handoff_dn) && ...
                isscalar(states(c).spinup_handoff_dn) && ...
                abs(states(c).spinup_handoff_dn - start_dn) < 1/48)
            handoff = 'an unrecorded instant';
            if isfield(states(c),'spinup_handoff_dn') && ...
                    isnumeric(states(c).spinup_handoff_dn) && ...
                    isscalar(states(c).spinup_handoff_dn)
                handoff = datestr(states(c).spinup_handoff_dn,'yyyy-mm-dd HH:MM');
            end
            error('Load_Field_Spinup_States:HandoffMismatch', ...
                ['Cell %d was spun up to hand off at %s; this run starts at ' ...
                 '%s.'],topo.cell_id(c),handoff,datestr(start_dn,'yyyy-mm-dd HH:MM'));
        end
    end
    Validate_BG_State_For_Transient(states(c),sprintf('%s cell %d', ...
        char(string(state_file)),topo.cell_id(c)));
end
end
