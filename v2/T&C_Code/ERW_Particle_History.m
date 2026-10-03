function varargout = ERW_Particle_History(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Daily blocks for the particle-class histories of an amended column.
%%
%% psd, d, lamb, SSA and delta_d are rows x size classes x applications:
%% one row per ERW step of the whole run, one slice per scheduled rock
%% application. Over 35 years (1840896 ten-minute rows, 35 applications)
%% each is about 3.6 GB. ERW_Daily_Step_2D used to write its rows into these
%% arrays in place, and since the frame still holds the column state while
%% the step runs, MATLAB copied every whole array on each daily call; and
%% ERW_function searched each whole history for every application at every
%% step. Neither changes a number, but together they made an amended
%% 35-year column several times slower than an unamended one. Here the
%% rows a day writes collect in a block, the frame writes the blocks into
%% the full arrays once, after the run, and the rows ERW_function reads
%% come from the blocks:
%%   [H,weath] = ERW_Particle_History('start',weath,weath_names)
%%   day = ERW_Particle_History('open',H,first_row,n_rows)
%%   v   = ERW_Particle_History('row',H,day,weath,name,r)    row r, 1 x classes x apps
%%   F   = ERW_Particle_History('first',H,day,name)          psd or lamb
%%   day = ERW_Particle_History('write',H,day,out,r)
%%   H   = ERW_Particle_History('track',H,day,out,r,weath)
%%   H   = ERW_Particle_History('close',H,day)
%%   E   = ERW_Particle_History('assemble',E)
%% 'first' returns, for each application, the first row whose last size
%% class is positive (zeros where none is): ERW_function's
%% psd_0(find(psd_0(:,end,k)>0,1),:,k) gives the same row from it as from
%% the whole history, which is that application's initial distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch action
    case 'start'
        [varargout{1},varargout{2}] = start(varargin{:});
    case 'open'
        [H,first_row,n_rows] = varargin{:};
        day = struct('first',first_row,'n',n_rows,'written',false(n_rows,1),'data',struct());
        for k = 1:numel(H.names)
            name = H.names{k};
            day.data.(name) = zeros(n_rows,H.size.(name)(2),H.size.(name)(3), ...
                'like',H.initial.(name));
        end
        varargout{1} = day;
    case 'row'
        varargout{1} = row(varargin{:});
    case 'first'
        [H,~,name] = varargin{:};
        varargout{1} = H.first.(name);
    case 'write'
        [H,day,out,r] = varargin{:};
        idx = r - day.first + 1;
        if idx < 1 || idx > day.n
            error('ERW_Particle_History:RowOutsideDay', ...
                'Row %d is outside the day''s rows %d-%d.',r,day.first,day.first+day.n-1);
        end
        for k = 1:numel(H.names)
            name = H.names{k};
            if isfield(out,name)
                day.data.(name)(idx,:,:) = out.(name);
            end
        end
        day.written(idx) = true;
        varargout{1} = day;
    case 'track'
        varargout{1} = track(varargin{:});
    case 'close'
        [H,day] = varargin{:};
        if any(day.written)
            H.blocks{end+1} = day;
        end
        varargout{1} = H;
    case 'assemble'
        E = varargin{1};
        if isstruct(E) && isfield(E,'hist')
            for k = 1:numel(E.hist.names)
                name = E.hist.names{k};
                E.weath.(name) = zeros(E.hist.size.(name),'like',E.hist.initial.(name));
                E.weath.(name)(E.hist.initial_rows.(name),:,:) = E.hist.initial.(name);
            end
            for b = 1:numel(E.hist.blocks)
                block = E.hist.blocks{b};
                rows = block.first + find(block.written) - 1;
                for k = 1:numel(E.hist.names)
                    name = E.hist.names{k};
                    E.weath.(name)(rows,:,:) = block.data.(name)(block.written,:,:);
                end
            end
            E = rmfield(E,'hist');
        end
        varargout{1} = E;
    otherwise
        error('ERW_Particle_History:UnknownAction','Unknown action %s.',action);
end
end

function [H,weath] = start(weath,weath_names)
if nargin < 2; weath_names = fieldnames(weath); end
candidates = {'psd','d','lamb','SSA','delta_d'};
keep = cellfun(@(n) isfield(weath,n) && any(strcmp(weath_names,n)) && ...
    size(weath.(n),3) > 1,candidates);
H.names = candidates(keep);
H.blocks = {};
H.size = struct();
H.initial = struct();
H.initial_rows = struct();
H.first = struct();
H.first_row = struct();
for k = 1:numel(H.names)
    name = H.names{k};
    H.size.(name) = size(weath.(name));
    rows = find(any(reshape(weath.(name),size(weath.(name),1),[]),2));
    H.initial_rows.(name) = rows;
    H.initial.(name) = weath.(name)(rows,:,:);
end
for name = intersect({'psd','lamb'},H.names)
    n = name{1};
    napp = size(weath.(n),3);
    H.first.(n) = zeros(1,size(weath.(n),2),napp);
    H.first_row.(n) = inf(1,napp);
    for a = 1:napp
        r = find(weath.(n)(:,end,a) > 0,1);
        if ~isempty(r)
            H.first.(n)(1,:,a) = weath.(n)(r,:,a);
            H.first_row.(n)(a) = r;
        end
    end
end
weath = rmfield(weath,H.names);
end

function v = row(H,day,~,name,r)
%% The latest value written to row r: today's block, then earlier blocks
%% from the most recent, then the array as initialised.
if r < 1 || r > H.size.(name)(1)
    error('ERW_Particle_History:RowOutsideHistory', ...
        'Row %d is outside the history''s rows 1-%d.',r,H.size.(name)(1));
end
idx = r - day.first + 1;
if idx >= 1 && idx <= day.n && day.written(idx)
    v = day.data.(name)(idx,:,:);
    return
end
for b = numel(H.blocks):-1:1
    block = H.blocks{b};
    idx = r - block.first + 1;
    if idx >= 1 && idx <= block.n && block.written(idx)
        v = block.data.(name)(idx,:,:);
        return
    end
end
idx = find(H.initial_rows.(name) == r,1);
if isempty(idx)
    v = zeros(1,H.size.(name)(2),H.size.(name)(3),'like',H.initial.(name));
else
    v = H.initial.(name)(idx,:,:);
end
end

function H = track(H,day,out,r,weath)
%% Keep each application's first row with a positive last class.
for name = intersect({'psd','lamb'},H.names)
    n = name{1};
    if ~isfield(out,n); continue; end
    value = day.data.(n)(r-day.first+1,:,:);
    for a = 1:numel(H.first_row.(n))
        positive = value(1,end,a) > 0;
        if positive && r < H.first_row.(n)(a)
            H.first.(n)(1,:,a) = value(1,:,a);
            H.first_row.(n)(a) = r;
        elseif r == H.first_row.(n)(a)
            if positive
                H.first.(n)(1,:,a) = value(1,:,a);
            else
                %% The first positive row was overwritten with a value that
                %% is not: search the whole history again.
                [H.first.(n)(1,:,a),H.first_row.(n)(a)] = rescan(H,day,weath,n,a);
            end
        end
    end
end
end

function [value,first_row] = rescan(H,day,weath,name,a)
column = zeros(H.size.(name)(1),1,'like',H.initial.(name));
column(H.initial_rows.(name)) = H.initial.(name)(:,end,a);
blocks = [H.blocks {day}];
for b = 1:numel(blocks)
    block = blocks{b};
    rows = block.first + find(block.written) - 1;
    column(rows) = block.data.(name)(block.written,end,a);
end
first_row = find(column > 0,1);
value = zeros(1,H.size.(name)(2),'like',H.initial.(name));
if isempty(first_row)
    first_row = inf;
else
    v = row(H,day,weath,name,first_row);
    value = v(1,:,a);
end
end
