function [t, names] = Conserved_Pool_Totals(B, names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conserved_Pool_Totals
%%
%% The conserved element totals of a 79-pool biogeochemical state, in
%% g element m-2 (the units BG_Observation_Targets produces).  These are
%% the quantities the random forest of the field spin-up predicts
%% (Lian et al. 2025, Fig. 3e) and the quantities its validation scores,
%% so both read them from here.
%%
%%   [t, names] = Conserved_Pool_Totals(B)         default eight totals
%%   t          = Conserved_Pool_Totals(B, names)  a chosen subset / order
%%
%% Pool indices follow BIOGEOCHEMISTRY_DYNAMIC3: 9-10 particulate organic
%% C, 11 mineral-associated organic C, 12-13 dissolved organic C, 18-19
%% microbial C, 23-34 the nitrogen pools (31 NH4, 32 NO3), 35-47 the
%% phosphorus pools.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2 || isempty(names)
    names = {'POC_gCm2','MAOC_gCm2','DOC_gCm2','MBC_gCm2', ...
        'total_N_gNm2','NH4_gNm2','NO3_gNm2','total_P_gPm2'};
end
B = B(:).';
if numel(B) < 47
    error('Conserved_Pool_Totals:ShortState', ...
        'A 79-pool state is required (got %d pools).', numel(B));
end
t = zeros(1, numel(names));
for k = 1:numel(names)
    switch names{k}
        case 'POC_gCm2',      t(k) = sum(B(9:10));
        case 'MAOC_gCm2',     t(k) = B(11);
        case 'DOC_gCm2',      t(k) = sum(B(12:13));
        case 'MBC_gCm2',      t(k) = sum(B(18:19));
        case 'total_N_gNm2',  t(k) = sum(B(23:34));
        case 'NH4_gNm2',      t(k) = B(31);
        case 'NO3_gNm2',      t(k) = B(32);
        case 'total_P_gPm2',  t(k) = sum(B(35:47));
        otherwise
            error('Conserved_Pool_Totals:UnknownTotal', ...
                'Unknown conserved total "%s".', names{k});
    end
end
end
