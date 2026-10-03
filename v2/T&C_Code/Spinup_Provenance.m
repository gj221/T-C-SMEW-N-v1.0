function key = Spinup_Provenance(cell_id,agronomic_identity,case_number,MOD_row,ERW_row,seed_label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spinup_Provenance
%%
%% The identity of one cell's spin-up: the cell, its soil inputs, the
%% workbook case, the MOD and ERW rows it runs (which carry the forcing,
%% meteo and transient windows) and what it starts from ('cold', or the
%% RF seed's digest).  Run_BG_Spinup stamps it on every phase checkpoint and
%% on the final state, and resumes or reuses only what carries the same key.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txt = sprintf('cell=%d;agronomic=%s;case=%d;mod=%s;erw=%s;seed=%s', ...
    cell_id, char(string(agronomic_identity)), case_number, ...
    ERW_Runtime_Parameter_Identity(MOD_row), ...
    ERW_Runtime_Parameter_Identity(ERW_row), char(string(seed_label)));
key = sha256_hex(txt);
end

function hex = sha256_hex(txt)
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(unicode2native(txt,'UTF-8'),'int8'));
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end
