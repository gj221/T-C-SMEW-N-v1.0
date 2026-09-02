function test_kriged_cell_parameterization()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_kriged_cell_parameterization
%%
%% Unit checks on the kriging -> per-cell parameterization bridge
%% (Krige_Attribute_At_Points, Apply_Kriged_Cell_Values). Uses synthetic
%% samples throughout so it does not depend on the Bruce Farm workbook.
%%
%%   1  Krige_Attribute_At_Points matches Krige_Attribute_Grid (the user's
%%      own function, untouched) at coincident points
%%   2  ordinary kriging is exact at a sample location
%%   3  the convex-hull flag correctly separates interpolation from
%%      extrapolation
%%   4  Apply_Kriged_Cell_Values promotes only cells that are both inside
%%      the sample hull and below the relative-std threshold
%%   5  a rejected cell keeps its template baseline untouched
%%   6  Apply_Kriged_Cell_Values refuses to silently overwrite an existing
%%      output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_kriged_cell_parameterization ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
v2   = fileparts(code);
addpath(code, fullfile(v2,'interpolation'), fullfile(v2,'spinup'));
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% --- synthetic sample cloud + fitted model --------------------------
rng(7);
xs = 40*rand(25,1); ys = 40*rand(25,1);
vs = 6 + 0.05*xs + 0.03*ys + 0.3*randn(25,1);
model = Fit_Ordinary_Kriging(xs,ys,vs,'exponential',struct('range',30));

%% ---------------------------------------------------------------------
%% 1. Matches Krige_Attribute_Grid at coincident points
%% ---------------------------------------------------------------------
gx = 0:5:40; gy = 0:5:40;
[pred_grid,var_grid] = Krige_Attribute_Grid(model,gx,gy,[]);
[GX,GY] = meshgrid(gx,gy);
[pred_pts,var_pts] = Krige_Attribute_At_Points(model,GX(:),GY(:));
check('Krige_Attribute_At_Points matches Krige_Attribute_Grid (prediction)', ...
    max(abs(pred_pts - pred_grid(:))) < 1e-9, ...
    sprintf('max|d| = %.3e',max(abs(pred_pts-pred_grid(:)))));
check('Krige_Attribute_At_Points matches Krige_Attribute_Grid (variance)', ...
    max(abs(var_pts - var_grid(:))) < 1e-9);

%% ---------------------------------------------------------------------
%% 2. Exactness at a sample location
%% ---------------------------------------------------------------------
[pred_at_sample,var_at_sample] = Krige_Attribute_At_Points(model,xs(3),ys(3));
check('ordinary kriging is exact at a sample location', ...
    abs(pred_at_sample - vs(3)) < 1e-8, ...
    sprintf('got %.6f expected %.6f',pred_at_sample,vs(3)));
check('kriging variance at a sample location is ~0', var_at_sample < 1e-6);

%% ---------------------------------------------------------------------
%% 3. Convex-hull extrapolation flag
%% ---------------------------------------------------------------------
[~,~,near_m,in_hull] = Krige_Attribute_At_Points(model, ...
    [mean(xs); 1000], [mean(ys); 1000]);
check('a point among the samples is inside the convex hull', in_hull(1));
check('a point far outside the sample cloud is flagged as extrapolation', ...
    ~in_hull(2));
check('nearest-sample distance is reported and sane', ...
    near_m(1) < 30 && near_m(2) > 900);

%% ---------------------------------------------------------------------
%% 4-6. Apply_Kriged_Cell_Values
%% ---------------------------------------------------------------------
tmp = tempname; mkdir(tmp);
n = 6;
Ttemplate = table((1:n)',nan(n,1),nan(n,1),false(n,1), ...
    'VariableNames',{'cell_id','baseline_bulk_pH','measured_bulk_pH', ...
    'use_measured_bulk_pH'});
Ttemplate.baseline_bulk_pH(:) = 6.0;    % pre-filled baseline, all cells
template_file = fullfile(tmp,'Cell_Agronomic_test.csv');
writetable(Ttemplate,template_file);

K = table((1:n)', ...
    [6.5;6.6;6.4;6.5;6.5;6.5], ...            % predicted
    [0.05;0.05;3.0;0.05;0.05;0.05], ...       % std: cell 3 far too noisy
    [true;true;true;true;false;true], ...     % in_sample_hull: cell 5 outside
    'VariableNames',{'cell_id','bulk_pH_predicted','bulk_pH_std', ...
    'bulk_pH_in_sample_hull'});

[out_file,log] = Apply_Kriged_Cell_Values(template_file,K, ...
    'MaxRelativeStd',0.1,'RequireInHull',true);
R = readtable(out_file,'VariableNamingRule','preserve');

check('cells inside the hull with low std are promoted to measured', ...
    all(R.use_measured_bulk_pH([1 2 4 6])) && ...
    isequal(R.measured_bulk_pH([1 2 4 6]),K.bulk_pH_predicted([1 2 4 6])));
check('the noisy in-hull cell is rejected on relative std', ...
    ~R.use_measured_bulk_pH(3));
check('the out-of-hull cell is rejected regardless of its std', ...
    ~R.use_measured_bulk_pH(5));
check('rejected cells keep their template baseline untouched', ...
    R.baseline_bulk_pH(3) == 6.0 && R.baseline_bulk_pH(5) == 6.0 && ...
    isnan(R.measured_bulk_pH(3)) && isnan(R.measured_bulk_pH(5)));
check('the log reports the right accept/reject counts', ...
    log.n_accepted(1) == 4 && log.n_rejected_std(1) == 1 && ...
    log.n_rejected_hull(1) == 1);

assert_error(@() Apply_Kriged_Cell_Values(template_file,K, ...
    'OutputFile',out_file), 'Apply_Kriged_Cell_Values:OutputExists');
check('refuses to silently overwrite an existing output file (see above)', true);

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_kriged_cell_parameterization: %d check(s) failed',n_fail);
end
end

%% =======================================================================
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
