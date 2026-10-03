function out = Run_RF_Validation_Experiment(shape_file, agronomic_file, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_RF_Validation_Experiment
%%
%% Score the random forest of the flux-tracking field spin-up against the
%% truth it stands in for.  The truth is a spin-up of EVERY cell of the
%% polygon, assembled: the flux-tracking spin-up with SubsampleFraction 1
%% (every cell tracked and replayed; Run_Field_Cell_Spinup), or Method
%% per_cell for the 1-D reference -- a cluster array makes either
%% affordable.  Three questions, three tables:
%%
%%   1  RF regression.  For each subsample fraction and seed: draw the
%%      topographically stratified subsample, train the forest on the
%%      subsample's truth totals, predict the held-out cells, and report
%%      the held-out relative error per conserved total (Lian et al.
%%      2025 target: < 10 %) beside the forest's own k-fold estimate.
%%      Fractions and seeds are the design; the seed spread is the
%%      sampling uncertainty of the method at that fraction.
%%   2  Base-saturation attainment.  For every cell that carried measured
%%      base saturation into its spin-up, the exchangeable Ca, Mg, K of
%%      the handoff exchanger seed (spinup_erw_handoff_diagnostic.f_CEC
%%      times the cell's CEC), which the transient's first ERW step
%%      installs, against the targets (cmolc kg-1).  The spin-up's own
%%      T&C exchange pools B(61/69/53) are not scored: that step
%%      re-initialises them from the seed.
%%   3  Reconstructed run (optional).  If a directory of
%%      cellstate_<id>.mat files from a flux-tracking run at a production
%%      fraction is given, every cell's conserved totals and exchangeable
%%      pools against its own truth, with the handoff-gate result and how
%%      the cell was obtained (replay, rf, fallback).
%%
%% Name/value
%%   'FieldFeatureIndices','FieldPolygonIndex','GridCellSizeM','MinCoverage'
%%                        as Run_Field_Cell_Spinup
%%   'FieldDEM'           required (the forest's covariates)
%%   'CaseNumber'         workbook case (default 5)
%%   'TruthStateFile'     assembled CellBGStates.mat of the truth (every
%%                        cell replayed, or per_cell)
%%   'Fractions'          subsample fractions (default [0.10 0.15 0.20 0.30])
%%   'Seeds'              subsample / forest seeds (default 1:5)
%%   'KFold'              forest cross-validation folds (default 5)
%%   'RFStatesDirectory'  cellstate files of a flux-tracking run at a
%%                        production fraction (optional)
%%   'OutputDirectory'    where the CSVs and figure go (default: the
%%                        truth file's directory)
%%   'MODParamColumns'    provenance only (recorded)
%%
%% OUTPUT struct: .regression (table), .base_saturation (table),
%%   .rf_path (table or []), .files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'), fullfile(v2,'T&C_Code','Inputs'), ...
    fullfile(v2,'T&C_Code','Parameters'), fullfile(v2,'interpolation'));
opt = struct('FieldFeatureIndices',[],'FieldPolygonIndex',1,'GridCellSizeM',[], ...
    'MinCoverage',[],'FieldDEM','','CaseNumber',5,'TruthStateFile','', ...
    'Fractions',[0.10 0.15 0.20 0.30],'Seeds',1:5,'KFold',5, ...
    'RFStatesDirectory','','OutputDirectory','','MODParamColumns',[]);
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_RF_Validation_Experiment:BadOption','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if isempty(opt.FieldDEM)
    error('Run_RF_Validation_Experiment:NeedsDEM','FieldDEM is required.');
end
if ~isfile(opt.TruthStateFile)
    error('Run_RF_Validation_Experiment:NoTruth', ...
        'TruthStateFile not found: %s', opt.TruthStateFile);
end
if isempty(opt.OutputDirectory); opt.OutputDirectory = fileparts(opt.TruthStateFile); end
if ~exist(opt.OutputDirectory,'dir'); mkdir(opt.OutputDirectory); end

%% --- grid, inputs, truth -----------------------------------------------
setups = Build_Shapefile_Grid_Setups(shape_file,opt.GridCellSizeM, ...
    opt.FieldFeatureIndices,opt.MinCoverage);
topo = setups(opt.FieldPolygonIndex);
dem = Load_Field_DEM(char(opt.FieldDEM),topo,'OutResolutionM', ...
    default_dem_resolution(opt.GridCellSizeM,topo));
topo = Attach_Topography_To_Topology(topo,dem);
resolved = Resolve_Cell_Agronomic_Inputs(topo, ...
    readtable(agronomic_file,'VariableNamingRule','preserve'));
workbook = fullfile(v2,'T&C_Code','Inputs','ERW_MOD_PARAM_parameters.xlsx');
[MOD,ERW] = Read_Parameter_Workbook(workbook);
MOD = MOD(MOD.Number == opt.CaseNumber,:); ERW = ERW(ERW.Number == opt.CaseNumber,:);
CellInputs = Configure_Cell_Agronomic_Inputs(MOD,ERW,resolved);
truth = Load_Field_Spinup_States(opt.TruthStateFile,topo, ...
    string({CellInputs.agronomic_identity}));
N = topo.N_cells;
[X, feature_names] = Cell_Topographic_Features(topo,resolved);
[~, target_names] = Conserved_Pool_Totals(truth(1).B_out(end,:));
Y = zeros(N,numel(target_names));
for c = 1:N; Y(c,:) = Conserved_Pool_Totals(truth(c).B_out(end,:),target_names); end
[~,shp] = fileparts(shape_file);
stem = sprintf('%s_polygon_%02d',shp,opt.FieldPolygonIndex);
files = struct();

%% --- 1. RF regression across fractions and seeds ------------------------
rows = {};
for f = opt.Fractions(:).'
    for sd = opt.Seeds(:).'
        sub = Select_Topographic_Subsample(topo,resolved,f,'Seed',sd);
        sub = sub(:); held = setdiff((1:N).',sub);
        row = struct('fraction',f,'seed',sd,'n_subsample',numel(sub), ...
            'n_heldout',numel(held),'trained',false);
        for t = 1:numel(target_names)
            row.(['heldout_rel_err_' target_names{t}]) = NaN;
            row.(['heldout_max_cell_rel_err_' target_names{t}]) = NaN;
            row.(['heldout_p90_cell_rel_err_' target_names{t}]) = NaN;
            row.(['heldout_skill_vs_mean_' target_names{t}]) = NaN;
            row.(['heldout_skill_vs_nn_' target_names{t}]) = NaN;
            row.(['cv_rel_err_' target_names{t}]) = NaN;
        end
        if numel(sub) >= 2*opt.KFold && ~isempty(held)
            rf = Train_Pool_Random_Forest(X(sub,:),Y(sub,:),feature_names, ...
                target_names,'Seed',sd,'KFold',opt.KFold);
            pred = rf.predict(X(held,:));
            %% Baselines the forest has to beat: the training mean, and the
            %% nearest subsample cell in standardised covariate space (the
            %% template the RF path would otherwise copy).
            mu = mean(X(sub,:),1); sg = std(X(sub,:),0,1); sg(sg == 0) = 1;
            [~,nn] = min(pdist2((X(held,:)-mu)./sg,(X(sub,:)-mu)./sg),[],2);
            pred_mean = repmat(mean(Y(sub,:),1),numel(held),1);
            pred_nn = Y(sub(nn),:);
            for t = 1:numel(target_names)
                obs = Y(held,t); scale = max(mean(abs(obs)),eps);
                row.(['heldout_rel_err_' target_names{t}]) = mean(abs(pred(:,t)-obs))/scale;
                cell_err = abs(pred(:,t)-obs)./max(abs(obs),eps);
                row.(['heldout_max_cell_rel_err_' target_names{t}]) = max(cell_err);
                row.(['heldout_p90_cell_rel_err_' target_names{t}]) = prctile(cell_err,90);
                mae = mean(abs(pred(:,t)-obs));
                row.(['heldout_skill_vs_mean_' target_names{t}]) = ...
                    1 - mae/max(mean(abs(pred_mean(:,t)-obs)),eps);
                row.(['heldout_skill_vs_nn_' target_names{t}]) = ...
                    1 - mae/max(mean(abs(pred_nn(:,t)-obs)),eps);
                row.(['cv_rel_err_' target_names{t}]) = rf.cv_rel_error(t);
            end
            row.trained = true;
        else
            fprintf(['Run_RF_Validation_Experiment: fraction %.2f seed %d gives %d ' ...
                'subsample cells, fewer than 2 x KFold = %d; not trained\n'], ...
                f,sd,numel(sub),2*opt.KFold);
        end
        rows{end+1} = row; %#ok<AGROW>
    end
end
regression = struct2table([rows{:}]);
files.regression = fullfile(opt.OutputDirectory,sprintf('RF_Validation_%s.csv',stem));
writetable(regression,files.regression);

%% --- 2. base-saturation attainment on the truth --------------------------
bs = table();
if ismember('base_saturation_source',resolved.Properties.VariableNames)
    meas = find(resolved.base_saturation_source == "measured");
    if ~isempty(meas)
        state = nan(numel(meas),3); target = zeros(numel(meas),3);
        for i = 1:numel(meas)
            c = meas(i);
            target(i,:) = [resolved.exch_Ca_cmolc_kg(c) resolved.exch_Mg_cmolc_kg(c) ...
                resolved.exch_K_cmolc_kg(c)];
            %% f_CEC = [Ca Mg K Na Al H] of the exchanger the seed built
            if isfield(truth(c),'spinup_erw_handoff_diagnostic') && ...
                    isstruct(truth(c).spinup_erw_handoff_diagnostic) && ...
                    isfield(truth(c).spinup_erw_handoff_diagnostic,'f_CEC') && ...
                    numel(truth(c).spinup_erw_handoff_diagnostic.f_CEC) == 6
                f = truth(c).spinup_erw_handoff_diagnostic.f_CEC;
                state(i,:) = f(1:3)*resolved.CEC_cmolc_kg(c);
            end
        end
        if any(isnan(state(:)))
            warning('Run_RF_Validation_Experiment:NoHandoffSeed', ...
                '%d measured cell(s) carry no handoff exchanger seed; scored NaN.', ...
                nnz(any(isnan(state),2)));
        end
        bs = table(resolved.cell_id(meas),target(:,1),state(:,1),target(:,2),state(:,2), ...
            target(:,3),state(:,3),(state-target)./max(abs(target),eps), ...
            'VariableNames',{'cell_id','Ca_target_cmolc_kg','Ca_state_cmolc_kg', ...
            'Mg_target_cmolc_kg','Mg_state_cmolc_kg','K_target_cmolc_kg','K_state_cmolc_kg', ...
            'rel_err_Ca_Mg_K'});
        files.base_saturation = fullfile(opt.OutputDirectory, ...
            sprintf('Base_Saturation_Attainment_%s.csv',stem));
        writetable(bs,files.base_saturation);
        fprintf(['Run_RF_Validation_Experiment: base saturation on %d measured cells, ' ...
            'max |rel. error| Ca %.3g, Mg %.3g, K %.3g\n'],numel(meas), ...
            max(abs(bs.rel_err_Ca_Mg_K(:,1))),max(abs(bs.rel_err_Ca_Mg_K(:,2))), ...
            max(abs(bs.rel_err_Ca_Mg_K(:,3))));
    end
end

%% --- 3. RF path against the truth (optional) ---------------------------
rf_path = [];
if ~isempty(opt.RFStatesDirectory)
    rrows = {};
    for c = 1:N
        f = fullfile(opt.RFStatesDirectory,sprintf('cellstate_%04d.mat',resolved.cell_id(c)));
        if ~isfile(f); continue; end
        S = load(f,'st','method');
        Yr = Conserved_Pool_Totals(S.st.B_out(end,:),target_names);
        r = struct('cell_id',resolved.cell_id(c),'method',string(S.method), ...
            'handoff_ready',isfield(S.st,'spinup_erw_handoff_ready') && ...
            logical(S.st.spinup_erw_handoff_ready));
        for t = 1:numel(target_names)
            r.(['rel_err_' target_names{t}]) = (Yr(t)-Y(c,t))/max(abs(Y(c,t)),eps);
        end
        Bt = truth(c).B_out(end,:); Br = S.st.B_out(end,:);
        r.rel_err_exch_Ca = (Br(61)-Bt(61))/max(abs(Bt(61)),eps);
        r.rel_err_exch_Mg = (Br(69)-Bt(69))/max(abs(Bt(69)),eps);
        r.rel_err_exch_K  = (Br(53)-Bt(53))/max(abs(Bt(53)),eps);
        rrows{end+1} = r; %#ok<AGROW>
    end
    if ~isempty(rrows)
        rf_path = struct2table([rrows{:}]);
        files.rf_path = fullfile(opt.OutputDirectory,sprintf('RF_Path_vs_Truth_%s.csv',stem));
        writetable(rf_path,files.rf_path);
        fprintf('Run_RF_Validation_Experiment: RF path scored on %d cells (gate pass %.0f%%)\n', ...
            height(rf_path),100*mean(rf_path.handoff_ready));
    else
        fprintf('Run_RF_Validation_Experiment: no cellstate files in %s\n',opt.RFStatesDirectory);
    end
end

%% --- figure: held-out error vs fraction, one line per target -------------
fig = figure('Color','w','Visible','off','Position',[80 80 900 420]);
try; fig.Theme = 'light'; catch; end
ax = axes(fig); hold(ax,'on');
fr = unique(regression.fraction);
for t = 1:numel(target_names)
    col = ['heldout_rel_err_' target_names{t}];
    m = arrayfun(@(f) mean(regression.(col)(regression.fraction == f),'omitnan'),fr);
    s = arrayfun(@(f) std(regression.(col)(regression.fraction == f),'omitnan'),fr);
    errorbar(ax,100*fr,100*m,100*s,'o-','LineWidth',1.2,'DisplayName', ...
        strrep(target_names{t},'_',' '));
end
yline(ax,10,'--','Lian et al. 10 %','LabelHorizontalAlignment','left');
xlabel(ax,'subsample fraction [%]'); ylabel(ax,'held-out relative error [%]');
title(ax,sprintf('%s: RF reconstruction against %d spun-up cells (mean +/- std over %d seeds)', ...
    stem,N,numel(opt.Seeds)),'Interpreter','none');
legend(ax,'Location','northeastoutside'); grid(ax,'on');
files.figure = fullfile(opt.OutputDirectory,sprintf('RF_Validation_%s.png',stem));
exportgraphics(fig,files.figure,'Resolution',150); close(fig);

out = struct('regression',regression,'base_saturation',bs,'rf_path',rf_path, ...
    'files',files,'target_names',{target_names},'feature_names',{feature_names}, ...
    'n_cells',N);
fprintf('Run_RF_Validation_Experiment: wrote %s\n',files.regression);
end

function r = default_dem_resolution(requested,topo)
r = requested;
if isempty(r) && isfield(topo,'cell_size_m') && ~isempty(topo.cell_size_m)
    r = topo.cell_size_m;
end
end
