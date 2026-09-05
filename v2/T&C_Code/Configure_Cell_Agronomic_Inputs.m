function CellInputs = Configure_Cell_Agronomic_Inputs(MOD_base,ERW_base,resolved,CellBGStates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply resolved agronomic properties to each cell's model input rows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~istable(MOD_base) || ~istable(ERW_base) || ~istable(resolved)
    error('Configure_Cell_Agronomic_Inputs:InvalidInput', ...
        'MOD, ERW and resolved agronomic inputs must be tables.');
end
n = height(resolved);
if ~any(height(MOD_base) == [1 n]) || ~any(height(ERW_base) == [1 n])
    error('Configure_Cell_Agronomic_Inputs:RowCount', ...
        'MOD and ERW inputs must contain one row or one row per cell.');
end
required_MOD = {'PH','Psan','Pcla','Porg'};
required_ERW = {'CEC_tot_in'};
if ~all(ismember(required_MOD,MOD_base.Properties.VariableNames)) || ...
        ~all(ismember(required_ERW,ERW_base.Properties.VariableNames))
    error('Configure_Cell_Agronomic_Inputs:MissingColumns', ...
        'MOD requires PH, Psan, Pcla and Porg; ERW requires CEC_tot_in.');
end
if nargin < 4
    CellBGStates = [];
end
if ~isempty(CellBGStates) && (~isstruct(CellBGStates) || numel(CellBGStates) ~= n)
    error('Configure_Cell_Agronomic_Inputs:StateCount', ...
        'CellBGStates must contain one state per resolved cell.');
end

CellInputs = repmat(struct('cell_id',0,'MOD_PARAM',table(), ...
    'ERW_PARAM',table(),'BG_STATE',struct(),'agronomic_input',table(), ...
    'agronomic_identity',''),n,1);
for c = 1:n
    im = min(c,height(MOD_base));
    ie = min(c,height(ERW_base));
    M = MOD_base(im,:);
    E = ERW_base(ie,:);
    M.PH = resolved.bulk_pH(c);
    M.Psan = resolved.sand_fraction(c);
    M.Pcla = resolved.clay_fraction(c);
    M.Porg = resolved.OM_fraction(c);
    E.CEC_tot_in = resolved.CEC_cmolc_kg(c);
    %% Measured exchangeable bases become spin-up observation targets for
    %% this cell (BG_Observation_Targets reads these MOD_PARAM columns and
    %% maps them onto the exchangeable pools); absent -> the spin-up's own
    %% exchanger stands, exactly as before.
    if ismember('exch_Ca_cmolc_kg',resolved.Properties.VariableNames) && ...
            isfinite(resolved.exch_Ca_cmolc_kg(c))
        M.InitExchCa_cmolc_kg = resolved.exch_Ca_cmolc_kg(c);
        M.InitExchMg_cmolc_kg = resolved.exch_Mg_cmolc_kg(c);
        M.InitExchK_cmolc_kg  = resolved.exch_K_cmolc_kg(c);
    end
    CellInputs(c).cell_id = resolved.cell_id(c);
    CellInputs(c).MOD_PARAM = M;
    CellInputs(c).ERW_PARAM = E;
    CellInputs(c).agronomic_input = resolved(c,:);
    CellInputs(c).agronomic_identity = char(resolved.input_identity(c));
    if ~isempty(CellBGStates)
        if CellBGStates(c).cell_id ~= resolved.cell_id(c) || ...
                ~strcmp(char(string(CellBGStates(c).cell_agronomic_identity)), ...
                CellInputs(c).agronomic_identity)
            error('Configure_Cell_Agronomic_Inputs:StateIdentity', ...
                'Spin-up state %d does not match resolved cell %d.',c,resolved.cell_id(c));
        end
        CellInputs(c).BG_STATE = CellBGStates(c);
    end
end
end
