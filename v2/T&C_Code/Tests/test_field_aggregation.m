function test_field_aggregation()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In-field totals use polygon overlap, while full cells remain simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
addpath(code_directory);

topo.cell_area_m2 = [100;100;100];
topo.cell_overlap_area_m2 = [75;25;0];
cell_values = [1 2 100; 3 4 200];

[field_mean,field_total] = Field_Aggregate_Areal(cell_values,topo);
[accounting_area,weights,fraction] = Field_Aggregation_Areas(topo);

assert(max(abs(accounting_area-[75;25;0])) < 1e-12);
assert(max(abs(weights-[0.75 0.25 0])) < 1e-12);
assert(max(abs(fraction-[0.75;0.25;0])) < 1e-12);
assert(max(abs(field_mean-[1.25;3.25])) < 1e-12);
assert(max(abs(field_total-[125;325])) < 1e-12);

field_pH = pH_Area_Weighted([6 7 9],weights);
expected_pH = -log10(0.75*10^-6+0.25*10^-7);
assert(abs(field_pH-expected_pH) < 1e-12);

legacy.cell_area_m2 = [3;1];
[legacy_mean,legacy_total] = Field_Aggregate_Areal([2 6],legacy);
assert(abs(legacy_mean-3) < 1e-12);
assert(abs(legacy_total-12) < 1e-12);

bad = topo;
bad.cell_overlap_area_m2 = [0;0;0];
expect_error(@() Field_Aggregation_Areas(bad), ...
    'Field_Aggregation_Areas:EmptyField');
bad.cell_overlap_area_m2 = [101;0;0];
expect_error(@() Field_Aggregation_Areas(bad), ...
    'Field_Aggregation_Areas:BadOverlapArea');
expect_error(@() Field_Aggregate_Areal(ones(2,2),topo), ...
    'Field_Aggregate_Areal:SizeMismatch');

driver_source = fileread(fullfile(code_directory,'MAIN_FRAME_BG_ERW_2D.m'));
assert(contains(driver_source, ...
    '[Domain.N2Oflx,Domain.N2O_total_gN] = Field_Aggregate_Areal'));
assert(contains(driver_source,'Domain.cell_overlap_area_m2'));
assert(contains(driver_source,'Domain.exit_water_field_total_m3'));

plot_source = fileread(fullfile(fileparts(code_directory), ...
    'plotting','Plot_1D_Then_All_12.m'));
assert(contains(plot_source,'Field_Aggregation_Areas(S.topo)'));

fprintf('test_field_aggregation: all checks passed\n');

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
