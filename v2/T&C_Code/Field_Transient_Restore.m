%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Transient_Restore  (script)
%%
%% Run by MAIN_FRAME_BG_ERW_2D, after it has allocated and seeded this
%% segment, when the caller hands it the end state of the segment before
%% (FIELD_TRANSIENT_CARRY, Field_Transient_Carry 'capture'). The carried
%% days become this segment's first hours, the last daily rows its first
%% daily rows, the columns' carried state and their ERW chemistry (on this
%% segment's axes) replace what the parameter files seeded, and the time
%% loop resumes at the first hour after the carried days, with the day
%% counter on the last of them. The segment's first hours must be the
%% carried hours (Field_Transient_Carry 'check').
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Field_Transient_Carry('check',FIELD_TRANSIENT_CARRY,Date,N_col);
ftc_names = Field_Transient_Carry('names');
ftc_L = FIELD_TRANSIENT_CARRY.overlap_hours;
ftc_D = FIELD_TRANSIENT_CARRY.overlap_days + 1;     % daily rows carried
for ftc_k = 1:numel(ftc_names.hourly)
    eval(sprintf('%s(1:%d,:) = FIELD_TRANSIENT_CARRY.hourly.%s;', ...
        ftc_names.hourly{ftc_k},ftc_L,ftc_names.hourly{ftc_k}));
end
for ftc_k = 1:numel(ftc_names.daily)
    eval(sprintf('%s(1:%d,:) = FIELD_TRANSIENT_CARRY.daily.%s;', ...
        ftc_names.daily{ftc_k},ftc_D,ftc_names.daily{ftc_k}));
end
for ftc_k = 1:numel(ftc_names.column)
    ftc_value = FIELD_TRANSIENT_CARRY.column.(ftc_names.column{ftc_k});
    if ~isequal(size(eval(ftc_names.column{ftc_k})),size(ftc_value))
        error('Field_Transient_Restore:ColumnState', ...
            '%s is [%s] here and [%s] in the carried state.',ftc_names.column{ftc_k}, ...
            num2str(size(eval(ftc_names.column{ftc_k}))),num2str(size(ftc_value)));
    end
    eval([ftc_names.column{ftc_k} ' = ftc_value;']);
end
for ftc_field = fieldnames(FIELD_TRANSIENT_CARRY.lateral_diag).'
    lateral_diag.(ftc_field{1})(1:ftc_D,:) = FIELD_TRANSIENT_CARRY.lateral_diag.(ftc_field{1});
end
jDay = FIELD_TRANSIENT_CARRY.jDay;
L_day = FIELD_TRANSIENT_CARRY.L_day;
for ftc_c = 1:N_col
    ERW{ftc_c} = Field_Transient_Carry('erw',FIELD_TRANSIENT_CARRY.ERW{ftc_c},NN,NNd,Ca);
end
j = ftc_D;
field_transient_first_step = ftc_L + 1;
fprintf('  field transient: carried state from %s, resuming at %s\n', ...
    datestr(Date(1),'yyyy-mm-dd'),datestr(Date(field_transient_first_step),'yyyy-mm-dd HH:MM'));
clear ftc_names ftc_L ftc_D ftc_k ftc_value ftc_field ftc_c
