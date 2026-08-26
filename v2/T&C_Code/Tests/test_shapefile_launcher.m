function test_shapefile_launcher()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shapefile options survive the parallel-launcher worker boundary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
v2_directory = fileparts(code_directory);
run_directory = fullfile(v2_directory,'run');
addpath(code_directory,run_directory);

shape_file = fullfile(v2_directory,'Spatial inputs','Bruce Farm Virginia', ...
    'applied_yield_shp.shp');
plan = Run_All_12_Parallel('Runs',1,'DryRun',true, ...
    'FieldShapefile',shape_file);

assert(numel(plan) == 1);
assert(contains(plan.result_file,'field_applied_yield_shp__geometry'));
assert(contains(plan.batch,'''FieldShapefile'''));
assert(contains(plan.batch,strrep(shape_file,'''','''''')));

explicit = Run_All_12_Parallel('Runs',1,'DryRun',true, ...
    'FieldShapefile',shape_file,'GridCellSizeM',100);
assert(contains(explicit.result_file,'field_applied_yield_shp__cell100m'));
assert(contains(explicit.batch,'''GridCellSizeM'',100'));

old_directory = cd(v2_directory);
directory_cleanup = onCleanup(@() cd(old_directory));
relative_shape = fullfile('Spatial inputs','Bruce Farm Virginia', ...
    'applied_yield_shp.shp');
relative_plan = Run_All_12_Parallel('Runs',1,'DryRun',true, ...
    'FieldShapefile',relative_shape);
assert(contains(relative_plan.batch,strrep(shape_file,'''','''''')));
clear directory_cleanup
cd(old_directory);

expect_error(@() Run_All_12_Parallel('Runs',4,'DryRun',true, ...
    'FieldShapefile',shape_file,'GridCellSizeM',100), ...
    'Run_All_12_Parallel:BadFieldLayout');
expect_error(@() Run_All_12_Parallel('Runs',1,'DryRun',true, ...
    'FieldShapefile',shape_file,'GridCellSizeM',100,'GridNx',2), ...
    'Run_All_12_Parallel:ConflictingGridOptions');
expect_error(@() Run_All_12_Parallel('Runs',1,'DryRun',true, ...
    'GridCellSizeM',100), ...
    'Run_All_12_Parallel:IncompleteFieldGrid');
expect_error(@() Run_One_Subfield(4,'FieldShapefile',shape_file), ...
    'Run_One_Subfield:BadFieldLayout');
expect_error(@() Run_One_Subfield(1,'GridCellSizeM',100), ...
    'Run_One_Subfield:IncompleteFieldGrid');
expect_error(@() Run_One_Subfield(1,'FieldShapefile',shape_file, ...
    'GridCellSizeM',100,'GridNx',2), ...
    'Run_One_Subfield:ConflictingGridOptions');

fprintf('test_shapefile_launcher: all checks passed\n');

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
