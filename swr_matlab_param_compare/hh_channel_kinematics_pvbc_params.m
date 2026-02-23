function hh = hh_channel_kinematics_pvbc_params(json_path)
% Load Channel_Kinematics PVBC parameters from SONATA dynamics JSON.
% Uses soma conductances and mod-file kinetics implemented in hh_compute_gating()
% kind='ck_pvbc' corresponding to na3n.mod + kdrbca1.mod + kdb.mod.

if nargin < 1 || isempty(json_path)
    json_path = fullfile('..','PVBC_dynamics_params_sonata.json');
end

raw = jsondecode(fileread(json_path));

hh.kind = 'ck_pvbc';
hh.name = 'Channel Kinematics PVBC (soma)';

hh.celsius = raw.conditions(1).celsius;
rev = raw.conditions(1).erev(1);
hh.ena = rev.ena;
hh.ek = rev.ek;

hh.gbar_na3  = find_param(raw.genome, 'soma', 'gbar_na3');
hh.gbar_kdrb = find_param(raw.genome, 'soma', 'gkdrbar_kdrb');
hh.gbar_kdb  = find_param(raw.genome, 'soma', 'gkdbar_kdb');

% default shifts from mod files
hh.sh_na3 = 0;
hh.sh_kdrb = 0;
hh.sh_kdb = 0;

end

function val = find_param(genome, section_name, param_name)
idx = find(arrayfun(@(g) strcmp(g.section, section_name) && strcmp(g.name, param_name), genome), 1, 'first');
if isempty(idx)
    error('Could not find %s in section %s.', param_name, section_name);
end
val = genome(idx).value;
end
