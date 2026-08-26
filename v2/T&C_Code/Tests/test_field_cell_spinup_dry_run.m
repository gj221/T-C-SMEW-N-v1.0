function test_field_cell_spinup_dry_run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The shapefile topology and cell soil rows produce spin-up configurations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
addpath(code_directory,fullfile(v2_directory,'spinup'),fullfile(v2_directory,'run'));
shape_file = fullfile(v2_directory,'Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
setups = Build_Shapefile_Grid_Setups(shape_file,[],[1 2]);
topo = setups(1);
n = topo.N_cells;
T = table(topo.cell_id,topo.cell_xc_projected_m,topo.cell_yc_projected_m, ...
    repmat(6.1,n,1),nan(n,1),false(n,1), ...
    repmat(0.10,n,1),repmat(0.35,n,1),nan(n,1),nan(n,1),false(n,1), ...
    repmat(0.025,n,1),nan(n,1),false(n,1), ...
    repmat(8.6,n,1),nan(n,1),false(n,1), ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction', ...
    'measured_sand_fraction','measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC'});
input_file = [tempname '.csv'];
cleanup = onCleanup(@() delete_if_present(input_file));
writetable(T,input_file);

[state_file,CellInputs] = Run_Field_Cell_Spinup(shape_file,input_file, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1,'DryRun',true);
assert(numel(CellInputs) == topo.N_cells);
assert(all([CellInputs.cell_id].' == topo.cell_id));
assert(all(arrayfun(@(x) x.MOD_PARAM.PH,CellInputs) == 6.1));
assert(contains(state_file,fullfile('polygon_01','CellBGStates.mat')));
assert(~isfile(state_file));
expect_error(@() Run_Field_Cell_Spinup(shape_file,input_file, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',0,'DryRun',true), ...
    'Run_Field_Cell_Spinup:BadPolygonIndex');

template_directory = tempname;
template_cleanup = onCleanup(@() remove_directory(template_directory));
files = Write_Field_Agronomic_Templates(shape_file,template_directory,[],[1 2]);
assert(numel(files) == 2);
assert(height(readtable(files{1})) == setups(1).N_cells);
assert(height(readtable(files{2})) == setups(2).N_cells);
state_files = {fullfile(template_directory,'state_1.mat'), ...
    fullfile(template_directory,'state_2.mat')};
dummy = true; %#ok<NASGU>
save(state_files{1},'dummy');
save(state_files{2},'dummy');
plan = Run_All_12_Parallel('Runs',1,'NDays',1,'DryRun',true, ...
    'FieldShapefile',shape_file,'FieldFeatureIndices',[1 2], ...
    'CellAgronomicFiles',files,'CellSpinupStateFiles',state_files, ...
    'PreferenceDirectoryRoot',template_directory);
assert(numel(plan) == 2 && isequal([plan.polygon_index],[1 2]));
assert(contains(plan(1).batch,files{1}) && contains(plan(2).batch,files{2}));
assert(contains(plan(1).result_file,Field_Input_Tag(files{1},state_files{1})));
assert(contains(plan(2).result_file,Field_Input_Tag(files{2},state_files{2})));
[~,runtime_leaf_1] = fileparts(plan(1).runtime_directory);
[~,runtime_leaf_2] = fileparts(plan(2).runtime_directory);
assert(strcmp(runtime_leaf_1,'run_01_polygon01'));
assert(strcmp(runtime_leaf_2,'run_01_polygon02'));

fprintf('test_field_cell_spinup_dry_run: all checks passed\n');
clear template_cleanup
remove_directory(template_directory)
clear cleanup
delete_if_present(input_file)
end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end

function remove_directory(path_name)
if isfolder(path_name)
    rmdir(path_name,'s')
end
end

function delete_if_present(file_name)
if isfile(file_name)
    delete(file_name)
end
end
