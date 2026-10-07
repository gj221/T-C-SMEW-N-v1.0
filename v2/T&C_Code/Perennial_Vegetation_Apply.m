%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perennial_Vegetation_Apply  (script)
%%
%% Run by a perennial stand's vegetation file after it sets its initial
%% stand.  When the run follows another -- a spin-up year, the handoff, the
%% transient -- the stand starts from that run's end state instead: Veg_out
%% (Perennial_Vegetation_Capture, or a handed-off state) sets row 1 of every
%% state Perennial_Vegetation_States lists.  The 2-D driver loads each cell's
%% parameters in a workspace of their own, where nothing is preallocated: a
%% state the parameter files have not created, or created as its first row
%% alone (x(1,:) = 0 makes a scalar there), becomes that row.  Without
%% Veg_out nothing changes: a cold start begins from the stand the file sets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('Veg_out','var') && isstruct(Veg_out) && isscalar(Veg_out) && ...
        ~isempty(fieldnames(Veg_out))
    veg_names = Perennial_Vegetation_States();
    veg_missing = veg_names(~isfield(Veg_out,veg_names));
    if ~isempty(veg_missing)
        error('Perennial_Vegetation_Apply:IncompleteState', ...
            'The vegetation state handed in lacks %s.',strjoin(veg_missing,', '));
    end
    for veg_k = 1:numel(veg_names)
        veg_v = Veg_out.(veg_names{veg_k});
        if ~(isnumeric(veg_v) && isreal(veg_v) && all(isfinite(veg_v(:))) && size(veg_v,1) == 1)
            error('Perennial_Vegetation_Apply:InvalidState', ...
                'The vegetation state %s handed in must be one finite row.',veg_names{veg_k});
        end
        if exist(veg_names{veg_k},'var') && size(eval(veg_names{veg_k}),1) > 1
            veg_x = eval(veg_names{veg_k});
            veg_x(1,:,:) = veg_v;
        else
            veg_x = veg_v;
        end
        eval([veg_names{veg_k} ' = veg_x;']);
    end
    clear veg_names veg_missing veg_k veg_v veg_x
end
