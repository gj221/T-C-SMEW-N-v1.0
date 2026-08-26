function test_field_spinup_states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Per-cell spin-up states must be converged and match resolved inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_field_spinup_states ===\n');
topo = struct('N_cells',2,'cell_id',[1;2]);
identities = ["cell-one";"cell-two"];
CellBGStates = [valid_state(2,"cell-two") valid_state(1,"cell-one")]; %#ok<NASGU>
state_file = [tempname '.mat'];
cleanup = onCleanup(@() delete_if_present(state_file));
save(state_file,'CellBGStates');

S = Load_Field_Spinup_States(state_file,topo,identities);
assert(isequal([S.cell_id],[1 2]));
assert(S(1).B_out(1) == 1 && S(2).B_out(1) == 2);

CellBGStates(1).cell_agronomic_identity = 'wrong'; %#ok<NASGU>
save(state_file,'CellBGStates');
assert_error(@() Load_Field_Spinup_States(state_file,topo,identities), ...
    'Load_Field_Spinup_States:InputIdentityMismatch');

CellBGStates = [valid_state(2,"cell-two") valid_state(1,"cell-one")]; %#ok<NASGU>
CellBGStates(2).spinup_converged = false;
save(state_file,'CellBGStates');
assert_error(@() Load_Field_Spinup_States(state_file,topo,identities), ...
    'TC_SMEW_N_run:UnconvergedBGState');

fprintf('  state identity, ordering and convergence checks passed\n');
clear cleanup
delete_if_present(state_file)
end

function S = valid_state(cell_id,identity)
B = ones(1,79);
B(1) = cell_id;
budget = struct('closed',true,'residual',0,'tolerance',1e-9);
budgets = struct('N',budget,'P',budget,'K',budget, ...
    'Ca',budget,'Mg',budget,'Si',budget);
S = struct('cell_id',cell_id,'cell_agronomic_identity',char(identity), ...
    'B_out',B,'O_out',zeros(1,8),'RexmyI_out',zeros(1,3), ...
    'NavlI_out',zeros(1,3),'NupI_L_out',zeros(1,3,1), ...
    'spinup_converged',true,'spinup_complete',true, ...
    'spinup_final_pass',true,'spinup_erw_handoff_evaluated',true, ...
    'spinup_erw_handoff_ready',true,'spinup_nutrient_budget_vers',12, ...
    'spinup_nutrient_budgets',budgets, ...
    'spinup_biogeochemical_soil_mass_kg_m2',388.5, ...
    'spinup_water_loss_mean_l_m2_d',1.2);
end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end

function delete_if_present(file_name)
if isfile(file_name)
    delete(file_name)
end
end
