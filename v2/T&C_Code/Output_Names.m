%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Output_Names  --  scenario directory and result-file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function names = Output_Names(base_scenario, id_location, run_config, ERW_region_coverage, ERW_PARAM_input)

base_scenario = clean_name(base_scenario);
id_location = clean_name(id_location);
descriptor = Run_Descriptor(run_config, ERW_region_coverage, ERW_PARAM_input);

if contains(base_scenario, "_domain_")
    scenario = base_scenario;
else
    scenario = base_scenario + "_" + descriptor;
end

names = struct();
names.descriptor = descriptor;
names.scenario = scenario;
names.result_file = "TC_ERW_results_" + id_location + "__" + descriptor + ".mat";
end

function txt = clean_name(txt)
txt = string(txt);
txt = strip(txt);
txt = regexprep(txt, '^[\\/]+', '');
txt = regexprep(txt, '[\\/]+$', '');
txt = regexprep(txt, '\s+', '_');
end
