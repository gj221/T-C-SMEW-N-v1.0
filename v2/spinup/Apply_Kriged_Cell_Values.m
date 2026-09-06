function [output_file,log] = Apply_Kriged_Cell_Values( ...
    template_file, kriged, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply_Kriged_Cell_Values
%%
%% Deliberate, threshold-gated merge of Krige_Cell_Agronomic_Values'
%% review table into a Cell_Agronomic_*.csv template (Write_Field_
%% Agronomic_Templates / Resolve_Cell_Agronomic_Inputs). A cell is
%% promoted from its baseline value to the kriged measured value only
%% when the estimate is NOT an extrapolation (inside the sample convex
%% hull) AND its standard deviation is acceptably small relative to the
%% prediction. That <stem>_std is Krige_Cell_Agronomic_Values' TOTAL
%% uncertainty -- interpolation error plus the variogram-range
%% uncertainty integrated out by Monte Carlo -- not the classical
%% kriging variance at one assumed range. Every other cell is left
%% exactly as the template had it (baseline, or a hand-entered
%% measurement) -- this never silently overwrites a value the user
%% already set.
%%
%% INPUTS
%%   template_file  a Cell_Agronomic_*.csv from Write_Field_Agronomic_Templates
%%   kriged         a table from Krige_Cell_Agronomic_Values (.table), or
%%                  the path to its written review CSV
%%
%% Name/value
%%   'MaxRelativeStd'  struct or scalar: reject a cell/attribute when
%%                     std/|prediction| exceeds this (default 0.35, i.e.
%%                     the kriging standard error is under 35% of the
%%                     predicted value). Pass a struct with a field named
%%                     after each template stem (e.g. .bulk_pH) to set it
%%                     per attribute.
%%   'RequireInHull'   reject an extrapolated cell regardless of std
%%                     (default true)
%%   'OutputFile'      default: <template stem>_parameterized.csv
%%   'Overwrite'       overwrite an existing OutputFile (default false)
%%
%% OUTPUTS
%%   output_file  path written
%%   log          table: template, n_cells, n_accepted, n_rejected_hull,
%%                n_rejected_std, n_missing_baseline_and_measured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('MaxRelativeStd',0.35,'RequireInHull',true, ...
    'OutputFile','','Overwrite',false, ...
    'OverwriteMeasured',false);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

if ~isfile(template_file)
    error('Apply_Kriged_Cell_Values:MissingTemplate', ...
        'Template not found: %s',char(string(template_file)));
end
T = readtable(template_file,'VariableNamingRule','preserve');
if ~ismember('cell_id',T.Properties.VariableNames)
    error('Apply_Kriged_Cell_Values:BadTemplate', ...
        '%s has no cell_id column.',template_file);
end

if istable(kriged)
    K = kriged;
elseif isfile(char(string(kriged)))
    K = readtable(char(string(kriged)),'VariableNamingRule','preserve');
else
    error('Apply_Kriged_Cell_Values:BadKriged', ...
        'kriged must be a table or an existing review CSV path.');
end
if ~ismember('cell_id',K.Properties.VariableNames)
    error('Apply_Kriged_Cell_Values:BadKriged', ...
        'The kriged table has no cell_id column.');
end

%% Attribute stems present as <stem>_predicted in K, matched to
%% measured_<stem> (always the full stem, Write_Field_Agronomic_Templates)
%% and use_measured_<toggle> in T.  baseline_<stem> is optional: the
%% base-saturation fractions are measured-only attributes (there is no
%% pedotransfer baseline for the exchanger; a cell without a measurement
%% has its exchanger derived from pH in the spin-up), so their absence of
%% a baseline column must not exclude them.
%% The toggle column uses a shorter alias for some stems (one
%% use_measured_texture covers both sand_fraction and clay_fraction; OM
%% and CEC drop their unit/kind suffix) -- resolve it rather than assume
%% it matches the stem verbatim.
pred_cols = K.Properties.VariableNames(endsWith( ...
    K.Properties.VariableNames,'_predicted'));
stems = erase(pred_cols,'_predicted');
toggle_of = containers.Map();
keep = false(size(stems));
for a = 1:numel(stems)
    s = stems{a};
    if ~ismember(sprintf('measured_%s',s),T.Properties.VariableNames)
        continue
    end
    tog = resolve_toggle_column(s,T.Properties.VariableNames);
    if isempty(tog); continue; end
    toggle_of(s) = tog;
    keep(a) = true;
end
stems = stems(keep);
if isempty(stems)
    error('Apply_Kriged_Cell_Values:NoMatchingAttributes', ...
        ['No kriged attribute matches a measured_/use_measured_ ' ...
         'pair in the template.']);
end

[~,it,ik] = intersect(T.cell_id,K.cell_id,'stable');
%% cell_id is 1..N within ONE polygon setup and is regenerated for every
%% (shapefile, cell size, feature) combination, so the same integer can
%% name different ground in two tables.  Refuse to merge unless the rows
%% also agree on the polygon and, when both sides carry them, on the
%% projected cell centres (the tolerance Resolve_Cell_Agronomic_Inputs
%% itself applies).
if all(ismember({'polygon_index'},T.Properties.VariableNames)) && ...
        all(ismember({'polygon_index'},K.Properties.VariableNames))
    if any(T.polygon_index(it) ~= K.polygon_index(ik))
        error('Apply_Kriged_Cell_Values:PolygonMismatch', ...
            ['The template and the kriged table disagree on polygon_index ' ...
             'for matched cells; they were built from different polygons.']);
    end
end
xy_cols = {'cell_x_projected_m','cell_y_projected_m'};
if all(ismember(xy_cols,T.Properties.VariableNames)) && ...
        all(ismember(xy_cols,K.Properties.VariableNames))
    dx = abs(T.cell_x_projected_m(it) - K.cell_x_projected_m(ik));
    dy = abs(T.cell_y_projected_m(it) - K.cell_y_projected_m(ik));
    if any(dx > 1e-6 | dy > 1e-6)
        error('Apply_Kriged_Cell_Values:CellCoordinateMismatch', ...
            ['Matched cell_ids sit at different projected centres (max ' ...
             '%.3g m apart): the kriged table was built on a different ' ...
             'grid than the template.'],max([dx;dy]));
    end
elseif ~all(ismember(xy_cols,K.Properties.VariableNames))
    warning('Apply_Kriged_Cell_Values:NoCoordinateCheck', ...
        ['The kriged table carries no projected cell centres, so the merge ' ...
         'is matched on cell_id alone; regenerate it with the current ' ...
         'Krige_Cell_Agronomic_Values to enable the geometry check.']);
end
%% Per stem: the gate.  A value is only WRITTEN in the group pass below,
%% once every attribute sharing the stem's toggle has passed, so the
%% three base-saturation fractions (or sand and clay) always come from one
%% source; and a cell whose toggle is already on carries a value someone
%% entered deliberately, which the merge never overwrites unless told to.
accept_of_stem = struct(); pred_of_stem = struct(); counts = struct();
for a = 1:numel(stems)
    s = stems{a};
    vn = matlab.lang.makeValidName(s);
    max_rel = threshold_for(opt.MaxRelativeStd,s);
    pred = K.(sprintf('%s_predicted',s))(ik);
    sd   = K.(sprintf('%s_std',s))(ik);
    in_hull_col = sprintf('%s_in_sample_hull',s);
    if ismember(in_hull_col,K.Properties.VariableNames)
        in_hull = logical(K.(in_hull_col)(ik));
    else
        in_hull = true(numel(ik),1);
    end
    rel_ok = isfinite(pred) & isfinite(sd) & ...
        (sd <= max_rel*max(abs(pred),eps));
    hull_ok = ~opt.RequireInHull | in_hull;
    range_ok = physical_range_ok(s,pred);
    pre_use = logical(T.(toggle_of(s))(it));   % toggle state before this call
    keep_existing = pre_use & ~opt.OverwriteMeasured;
    accept = rel_ok & hull_ok & range_ok & ~keep_existing;
    accept_of_stem.(vn) = accept;
    pred_of_stem.(vn) = pred;
    baseline_col = sprintf('baseline_%s',s);
    has_baseline = ismember(baseline_col,T.Properties.VariableNames);
    if has_baseline
        still_bare = ~accept & ~pre_use & isnan(T.(baseline_col)(it));
    else
        still_bare = false(numel(it),1);       % measured-only: the model's own value stands
    end
    counts.(vn) = struct('n_cells',numel(ik),'n_accepted',nnz(accept), ...
        'n_rejected_hull',nnz(~hull_ok),'n_rejected_std',nnz(hull_ok & ~rel_ok), ...
        'n_rejected_range',nnz(hull_ok & rel_ok & ~range_ok), ...
        'n_kept_existing',nnz(keep_existing),'n_written',0, ...
        'n_missing_baseline_and_measured',nnz(still_bare),'has_baseline',has_baseline);
end

%% A toggle shared by several stems (use_measured_texture covers both
%% sand_fraction and clay_fraction; use_measured_base_saturation the
%% three fractions) is only set true, and its values only written, for
%% cells every co-toggled attribute accepted -- a good sand estimate does
%% not license turning on an unreviewed clay estimate.  The co-toggle
%% group is taken from the TEMPLATE, not from whichever stems happen to
%% be in the kriged table: if the template has both texture stems but the
%% sample file carried only sand, the toggle must stay off, because
%% switching it on would also switch on an unreviewed (NaN) measured clay.
unique_toggles = unique(cellfun(@(s) toggle_of(s),stems,'UniformOutput',false));
for u = 1:numel(unique_toggles)
    use_col = unique_toggles{u};
    group = template_toggle_group(use_col,T.Properties.VariableNames);
    missing_from_K = setdiff(group,stems);
    if ~isempty(missing_from_K)
        warning('Apply_Kriged_Cell_Values:IncompleteToggleGroup', ...
            ['%s also covers %s, which the kriged table does not carry; ' ...
             'the toggle is left unchanged for every cell.'], ...
            use_col,strjoin(missing_from_K,', '));
        continue
    end
    joint_accept = true(numel(ik),1);
    group_sum = zeros(numel(ik),1);
    for c = 1:numel(group)
        gvn = matlab.lang.makeValidName(group{c});
        joint_accept = joint_accept & accept_of_stem.(gvn);
        group_sum = group_sum + pred_of_stem.(gvn);
    end
    if strcmp(use_col,'use_measured_base_saturation')
        %% the three fractions must leave room for the acid saturation
        joint_accept = joint_accept & group_sum < 1;
    end
    for c = 1:numel(group)
        gvn = matlab.lang.makeValidName(group{c});
        measured_col = sprintf('measured_%s',group{c});
        Tv = T.(measured_col);
        Tv(it(joint_accept)) = pred_of_stem.(gvn)(joint_accept);
        T.(measured_col) = Tv;
        counts.(gvn).n_written = nnz(joint_accept);
    end
    Uv = T.(use_col);
    Uv(it(joint_accept)) = true;
    T.(use_col) = Uv;
end

log_rows = cell(numel(stems),9);
for a = 1:numel(stems)
    s = stems{a}; cs = counts.(matlab.lang.makeValidName(s));
    log_rows(a,:) = {s,cs.n_cells,cs.n_accepted,cs.n_written,cs.n_rejected_hull, ...
        cs.n_rejected_std,cs.n_rejected_range,cs.n_kept_existing, ...
        cs.n_missing_baseline_and_measured};
    if cs.has_baseline
        tail = sprintf('%d cell(s) still have neither a baseline nor a measured value', ...
            cs.n_missing_baseline_and_measured);
    else
        tail = 'no baseline column: unaccepted cells keep the model''s own value';
    end
    fprintf(['Apply_Kriged_Cell_Values: %-14s accepted %4d / %4d cells, wrote %4d ' ...
        '(rejected: %d outside sample hull, %d over relative-std %.0f%%, ' ...
        '%d outside physical range; %d kept as already measured); %s\n'], ...
        s,cs.n_accepted,cs.n_cells,cs.n_written,cs.n_rejected_hull,cs.n_rejected_std, ...
        100*threshold_for(opt.MaxRelativeStd,s),cs.n_rejected_range,cs.n_kept_existing,tail);
end
log = cell2table(log_rows,'VariableNames',{'attribute','n_cells', ...
    'n_accepted','n_written','n_rejected_hull','n_rejected_std','n_rejected_range', ...
    'n_kept_existing','n_missing_baseline_and_measured'});

if numel(it) < height(T)
    warning('Apply_Kriged_Cell_Values:UnmatchedCells', ...
        ['%d of %d template cells have no matching row in the kriged ' ...
         'table and were left untouched.'],height(T)-numel(it),height(T));
end

output_file = char(opt.OutputFile);
if isempty(output_file)
    [d,n,e] = fileparts(template_file);
    output_file = fullfile(d,[n '_parameterized' e]);
end
if isfile(output_file) && ~opt.Overwrite
    error('Apply_Kriged_Cell_Values:OutputExists', ...
        ['%s already exists. Pass Overwrite=true to replace it, or an ' ...
         'explicit OutputFile.'],output_file);
end
writetable(T,output_file);
fprintf('Apply_Kriged_Cell_Values: wrote %s\n',output_file);
end

%% =======================================================================
function v = threshold_for(spec,stem)
if isstruct(spec)
    if isfield(spec,stem); v = spec.(stem); else; v = 0.35; end
else
    v = double(spec);
end
end

%% =======================================================================
function group = template_toggle_group(use_col,template_columns)
%% Every stem the TEMPLATE routes to this toggle: the canonical aliases
%% plus any exact-name stem, kept only if the template carries its
%% measured_ column (a baseline_ column is optional, see above).
alias = containers.Map( ...
    {'bulk_pH','sand_fraction','clay_fraction','OM_fraction','CEC_cmolc_kg', ...
     'bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction'}, ...
    {'bulk_pH','texture',      'texture',      'OM',         'CEC', ...
     'base_saturation','base_saturation','base_saturation'});
short = erase(use_col,'use_measured_');
group = {};
for k = keys(alias)
    stem = k{1};
    if strcmp(alias(stem),short) || strcmp(stem,short)
        if ismember(sprintf('measured_%s',stem),template_columns)
            group{end+1} = stem; %#ok<AGROW>
        end
    end
end
if isempty(group) && ismember(sprintf('measured_%s',short),template_columns)
    group = {short};
end
end

%% =======================================================================
function col = resolve_toggle_column(stem,template_columns)
%% Write_Field_Agronomic_Templates names its toggle columns shorter than
%% the baseline_/measured_ stem for three of the five canonical
%% attributes: OM_fraction -> use_measured_OM, CEC_cmolc_kg ->
%% use_measured_CEC, and {sand_fraction,clay_fraction} share one
%% use_measured_texture. bulk_pH's toggle matches its stem exactly.
alias = containers.Map( ...
    {'bulk_pH','sand_fraction','clay_fraction','OM_fraction','CEC_cmolc_kg', ...
     'bs_Ca_fraction','bs_Mg_fraction','bs_K_fraction'}, ...
    {'bulk_pH','texture',      'texture',      'OM',         'CEC', ...
     'base_saturation','base_saturation','base_saturation'});
exact = sprintf('use_measured_%s',stem);
if ismember(exact,template_columns)
    col = exact;
elseif isKey(alias,stem) && ...
        ismember(sprintf('use_measured_%s',alias(stem)),template_columns)
    col = sprintf('use_measured_%s',alias(stem));
else
    col = '';
end
end

%% =======================================================================
function ok = physical_range_ok(stem,pred)
%% The gate a kriged value must clear before it can become an input:
%% fractions in (0,1), pH in (0,14), CEC positive.  Ordinary kriging
%% weights can be negative, so a small-valued attribute (K saturation of
%% 1-3 %) can be predicted below zero with a small std and would otherwise
%% pass the relative-std gate and then abort the whole polygon downstream.
if endsWith(stem,'_fraction')
    ok = isfinite(pred) & pred > 0 & pred < 1;
elseif strcmp(stem,'bulk_pH')
    ok = isfinite(pred) & pred > 0 & pred < 14;
elseif strcmp(stem,'CEC_cmolc_kg')
    ok = isfinite(pred) & pred > 0;
else
    ok = isfinite(pred);
end
end
