function test_erw_seed_water_loss()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERW_Seed_Water_Loss: the mean water loss [L m-2 d-1] the rainwater Na
%% seed divides the rain Na flux by (ERW_MOD_PARAM_excel).
%%   1  a value handed in (seed_water_loss_l_m2_d) wins
%%   2  else mean(Tw + Lw) of the forcing, Tw and Lw each day's total [m/d]:
%%      from a struct, and -- as the plot-scale driver passes it -- from a
%%      timetable or table, to the same value
%%   3  missing, non-finite or non-positive forcing gives NaN, which the
%%      Na seed refuses rather than use
%%   4  the plot-scale driver does pass a timetable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_directory = fileparts(mfilename('fullpath'));
code = fileparts(test_directory);
addpath(code);
Tw = [1.2e-3; 2.0e-3; 0.5e-3];   % [m/d]
Lw = [0.3e-3; 0.0;    1.1e-3];
expected = mean(Tw+Lw)*1000;     % [L m-2 d-1]

%% --- 1
assert(ERW_Seed_Water_Loss(struct('seed_water_loss_l_m2_d',2.5),struct('Tw',Tw,'Lw',Lw)) == 2.5, ...
    'A water loss handed in must win.');

%% --- 2
S = struct('Tw',Tw,'Lw',Lw);
t = datetime(2020,1,1) + days(0:2)';
TT = timetable(t,Tw,Lw,'VariableNames',{'Tw','Lw'});
TT.s = [0.4;0.5;0.6];             % other forcing the driver carries
T = timetable2table(TT);
assert(abs(ERW_Seed_Water_Loss(struct(),S)-expected) < 1e-12, ...
    'From a struct: mean(Tw + Lw) in L m-2 d-1.');
assert(abs(ERW_Seed_Water_Loss(struct(),TT)-expected) < 1e-12 && ...
    abs(ERW_Seed_Water_Loss([],T)-expected) < 1e-12, ...
    'A timetable or table must give the struct''s value.');

%% --- 3
assert(isnan(ERW_Seed_Water_Loss(struct(),removevars(TT,'Lw'))) && ...
    isnan(ERW_Seed_Water_Loss(struct(),struct('Tw',Tw))) && ...
    isnan(ERW_Seed_Water_Loss(struct(),[])), 'Without Tw and Lw there is no water loss.');
TTn = TT; TTn.Tw(2) = NaN;
TTz = TT; TTz.Tw(:) = 0; TTz.Lw(:) = 0;
assert(isnan(ERW_Seed_Water_Loss(struct(),TTn)) && isnan(ERW_Seed_Water_Loss(struct(),TTz)), ...
    'Non-finite or zero forcing gives no water loss.');

%% --- 4
frame = fileread(fullfile(code,'MAIN_FRAME_BG_ERW_10mins.m'));
assert(contains(frame,'TT_m = synchronize(') && ...
    contains(frame,'ERW_MOD_PARAM_excel(t, TT_m, ERW_TC_const'), ...
    'The plot-scale driver hands ERW_MOD_PARAM_excel a timetable.');

fprintf('test_erw_seed_water_loss: all checks passed\n');
end
