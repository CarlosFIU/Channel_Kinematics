function hh = hh_channel_kinematics_pc_params(json_path)
% Load Channel_Kinematics PC parameters from SONATA dynamics JSON.
% Uses soma conductances and mod-file kinetics implemented in hh_compute_gating()
% kind='ck_pc' corresponding to naxn.mod + kdrca1.mod.

if nargin < 1 || isempty(json_path)
    json_path = fullfile('..','PC_dynamics_params_sonata.json');
end

raw = jsondecode(fileread(json_path));

hh.kind = 'ck_pc';
hh.name = 'Channel Kinematics PC (soma)';

hh.celsius = raw.conditions(1).celsius;
rev = raw.conditions(1).erev(1);
hh.ena = rev.ena;
hh.ek = rev.ek;

hh.gbar_nax = find_param(raw.genome, 'soma', 'gbar_nax');
hh.gbar_kdr = find_param(raw.genome, 'soma', 'gkdrbar_kdr');

% default shifts from mod files
hh.sh_nax = 0;
hh.sh_kdr = 0;

end

function val = find_param(genome, section_name, param_name)
idx = find(arrayfun(@(g) strcmp(g.section, section_name) && strcmp(g.name, param_name), genome), 1, 'first');
if isempty(idx)
    error('Could not find %s in section %s.', param_name, section_name);
end
val = genome(idx).value;
end
