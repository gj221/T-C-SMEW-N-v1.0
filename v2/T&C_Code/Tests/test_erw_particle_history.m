function test_erw_particle_history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The particle-class histories of an amended column (rows x classes x
%% applications) are kept in daily blocks during the run and assembled at
%% its end, so that ERW_Daily_Step_2D does not copy every full history
%% array each day. Everything ERW_function reads must stay exactly what
%% the full arrays would give: the previous row, and the first row of each
%% application whose last class is positive (its initial size
%% distribution, which ERW_function finds with
%% psd_0(find(psd_0(:,end,k)>0,1),:,k)). The assembled arrays must equal
%% the full arrays written in place.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

nn = 60; ncl = 3; napp = 4;
names = {'psd','d','lamb','SSA','delta_d'};
full = struct();
for k = 1:numel(names)
    full.(names{k}) = zeros(nn,ncl,napp);
end
%% Row 1 is the initial state; application 1 is present from the start.
for k = 1:numel(names)
    full.(names{k})(1,:,1) = k + [0.1 0.2 0.3];
end
full.EW = zeros(nn,2);            % two-dimensional fields stay in place
full.psd(40,:,4) = [0.5 0.6 0.7];
full.lamb(40,:,4) = [1.5 1.6 1.7];
[H,weath] = ERW_Particle_History('start',full);
assert(~any(isfield(weath,names)), ...
    'The live weathering state must not retain the full particle histories.');
assert(isequal(sort(H.names),sort(names)), ...
    'Only the three-dimensional histories are kept in blocks.');

%% Before the first application ERW_function returns scalar zeros; the
%% original array assignment expands them across all scheduled cohorts.
scalar_out = struct();
for k = 1:numel(names); scalar_out.(names{k}) = 0; end
before_day = ERW_Particle_History('open',H,2,1);
before_day = ERW_Particle_History('write',H,before_day,scalar_out,2);
before_first = ERW_Particle_History('track',H,before_day,scalar_out,2,weath);
assert(isequal(before_first.first,H.first) && ...
    isequal(before_first.first_row,H.first_row), ...
    'Scalar zeros before an application must preserve the first-row cache.');
assert(isequal(ERW_Particle_History('row',H,before_day,weath,'psd',2), ...
    zeros(1,ncl,napp)), ...
    'A scalar solver output must expand as in the full history.');

rng(7);
%% Days of five rows from row 2, then overwrites of rows 20-24 (the final
%% flush can revisit a window) and of rows 12-14, whose row 13 is the first
%% positive row of application 3 and loses its sign, so the history must
%% be searched again.
windows = {2:6, 7:11, 12:16, 17:21, 22:26, 20:24, 12:14, 27:31, 32:36};
for w = 1:numel(windows)
    rows = windows{w};
    day = ERW_Particle_History('open',H,rows(1),numel(rows));
    for r = rows
        %% What ERW_function would receive at this step.
        for k = {'psd','d','lamb','SSA'}
            expected = full.(k{1})(r-1,:,:);
            got = ERW_Particle_History('row',H,day,weath,k{1},r-1);
            assert(isequal(got,expected), ...
                'Previous row of %s at row %d differs from the full history.',k{1},r);
        end
        for k = {'psd','lamb'}
            first = ERW_Particle_History('first',H,day,k{1});
            for a = 1:napp
                expected = full.(k{1})(find(full.(k{1})(:,end,a)>0,1),:,a);
                got = first(find(first(:,end,a)>0,1),:,a);
                assert(isequal(got,expected), ...
                    'The initial %s of application %d at row %d differs.',k{1},a,r);
            end
        end
        %% This step's output: applications start at rows 5, 13 and 29,
        %% and application 2 goes non-positive in the last class at row 33.
        out = struct();
        for k = 1:numel(names)
            v = zeros(1,ncl,napp);
            v(1,:,1) = rand(1,ncl);
            if r >= 5;  v(1,:,2) = rand(1,ncl); end
            if r >= 13; v(1,:,3) = rand(1,ncl); end
            if r >= 29; v(1,:,4) = rand(1,ncl); end
            if w == numel(windows) && r == 33; v(1,end,2) = 0; end
            out.(names{k}) = v;
        end
        if w == 6 && r == 21
            out.psd(1,end,2) = -1;     % an overwritten row may lose its sign
        end
        if w == 7 && r == 13
            out.psd(1,end,3) = 0; out.lamb(1,end,3) = -2;
        end
        for k = 1:numel(names)
            full.(names{k})(r,:,:) = out.(names{k});
        end
        day = ERW_Particle_History('write',H,day,out,r);
        H = ERW_Particle_History('track',H,day,out,r,weath);
    end
    H = ERW_Particle_History('close',H,day);
    if w == 4
        saved = struct('weath',weath,'hist',H);
        checkpoint = [tempname '.mat'];
        cleanup_file = onCleanup(@() delete(checkpoint));
        save(checkpoint,'saved','-v7.3');
        loaded = load(checkpoint,'saved');
        assert(isequaln(saved,loaded.saved), ...
            'Saving and reloading the history must preserve the live state.');
        H = loaded.saved.hist; weath = loaded.saved.weath;
    end
end

E = struct('weath',weath,'hist',H);
E = ERW_Particle_History('assemble',E);
assert(~isfield(E,'hist'),'Assembly must leave no block store in the column state.');
for k = 1:numel(names)
    assert(isequal(E.weath.(names{k}),full.(names{k})), ...
        'Assembled %s differs from the history written in place.',names{k});
end
assert(isequal(E.weath.EW,full.EW),'Fields kept in place must be untouched.');
E2 = ERW_Particle_History('assemble',struct('weath',full));
assert(isequal(E2.weath,full),'A column without blocks is left as it is.');
empty_day = ERW_Particle_History('open',H,nn+1,0);
assert(isequal(ERW_Particle_History('close',H,empty_day),H), ...
    'An empty final window must not add a history block.');
try
    ERW_Particle_History('write',H,empty_day,out,nn+1);
    error('test_erw_particle_history:MissingError','An out-of-window write succeeded.');
catch ME
    assert(strcmp(ME.identifier,'ERW_Particle_History:RowOutsideDay'), ...
        'An out-of-window write must report its row bounds.');
end
for r = [0 nn+1]
    try
        ERW_Particle_History('row',H,empty_day,weath,'psd',r);
        error('test_erw_particle_history:MissingError','An invalid history row was read.');
    catch ME
        assert(strcmp(ME.identifier,'ERW_Particle_History:RowOutsideHistory'), ...
            'An invalid history row must report the history bounds.');
    end
end
single_application = struct('psd',full.psd(:,:,1),'lamb',full.lamb(:,:,1));
[no_blocks,unchanged] = ERW_Particle_History('start',single_application);
assert(isempty(no_blocks.names) && isequal(unchanged,single_application), ...
    'A single application must retain the existing two-dimensional storage.');

%% The daily step and the frame use it; ERW_function is unchanged.
step = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
frame = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
assert(contains(step,"ERW_Particle_History('start',") && ...
    contains(step,"ERW_Particle_History('first',") && ...
    contains(step,"ERW_Particle_History('write',") && ...
    contains(step,"ERW_Particle_History('close',") && ...
    contains(frame,"ERW{c} = ERW_Particle_History('assemble',ERW{c});"), ...
    'The daily step must keep the histories in blocks and the frame assemble them.');
fprintf('PASS: particle-class histories kept in blocks read and assemble as in place.\n');
end
