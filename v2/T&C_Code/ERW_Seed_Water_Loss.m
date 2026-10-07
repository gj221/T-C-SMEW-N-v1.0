function Q = ERW_Seed_Water_Loss(ERW_TC_const_input,TT_m_input)
Q = NaN;
if isstruct(ERW_TC_const_input) && ...
        isfield(ERW_TC_const_input,'seed_water_loss_l_m2_d')
    q_input = ERW_TC_const_input.seed_water_loss_l_m2_d;
    if isscalar(q_input) && isfinite(q_input) && q_input > 0
        Q = q_input;
        return
    end
end

%%% v1.0: the plot-scale driver hands its forcing as a timetable (TT_m, Tw
%%% and Lw each day's total [m/d] at every step), which isfield does not
%%% see; a table or timetable with Tw and Lw is read as a struct is.
has_series = (isstruct(TT_m_input) && isfield(TT_m_input,'Tw') && ...
        isfield(TT_m_input,'Lw')) || ...
    ((istable(TT_m_input) || istimetable(TT_m_input)) && ...
        all(ismember({'Tw','Lw'},TT_m_input.Properties.VariableNames)));
if has_series
    Tw = TT_m_input.Tw(:);
    Lw = TT_m_input.Lw(:);
    if ~isempty(Tw) && isequal(size(Tw),size(Lw)) && ...
            all(isfinite(Tw)) && all(isfinite(Lw))
        q_forcing = mean(Tw+Lw)*1000;
        if isfinite(q_forcing) && q_forcing > 0
            Q = q_forcing;
        end
    end
end
end
