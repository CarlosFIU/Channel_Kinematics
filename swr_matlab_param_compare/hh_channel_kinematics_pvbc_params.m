function hh = hh_channel_kinematics_pvbc_params(json_path, section_name)
% Load Channel_Kinematics PVBC parameters from SONATA dynamics JSON.
% Uses conductances for the requested compartment section and
% mod-file kinetics implemented in hh_compute_gating().
% kind='ck_pvbc' corresponding to na3n.mod + kdrbca1.mod + kdb.mod.

if nargin < 1 || isempty(json_path)
    this_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(this_dir);
    json_path = fullfile(repo_root, 'PVBC_dynamics_params_sonata.json');
elseif ~isfile(json_path)
    this_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(this_dir);
    alt_path = fullfile(repo_root, json_path);
    if isfile(alt_path)
        json_path = alt_path;
    end
end

if ~isfile(json_path)
    error('Could not open PVBC dynamics JSON: %s', json_path);
end
if nargin < 2 || isempty(section_name)
    section_name = 'soma';
end

raw = jsondecode(fileread(json_path));

hh.kind = 'ck_pvbc';
hh.name = sprintf('Channel Kinematics PVBC (%s)', section_name);
hh.section = section_name;

hh.celsius = raw.conditions(1).celsius;
rev = raw.conditions(1).erev(1);
hh.ena = rev.ena;
hh.ek = rev.ek;

hh.gbar_na3  = find_param_or_nan(raw.genome, section_name, 'gbar_na3');
hh.gbar_kdrb = find_param(raw.genome, section_name, 'gkdrbar_kdrb');

% gkdb is present only in soma in this model.
hh.gbar_kdb  = find_param_or_nan(raw.genome, section_name, 'gkdbar_kdb');

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

function val = find_param_or_nan(genome, section_name, param_name)
idx = find(arrayfun(@(g) strcmp(g.section, section_name) && strcmp(g.name, param_name), genome), 1, 'first');
if isempty(idx)
    val = nan;
    return;
end
val = genome(idx).value;
end
