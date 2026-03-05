function swr_plot_pc_pvbc_targeted_kinematics_compartment_series(V)
% Generate complete targeted-kinematics figure sets for requested compartments.
% Figure 18: PC soma vs PVBC soma
% Figure 19: PC apic only
% Figure 20: PC dend (basal) only
% Figure 21: PVBC dend only

if nargin < 1 || isempty(V)
    V = linspace(-100, 50, 600);
end

runs = { ...
    struct('pc', 'soma', 'pv', 'soma', 'fig_num', 18, 'scope', 'both'), ...
    struct('pc', 'apic', 'pv', '',     'fig_num', 19, 'scope', 'pc_only'), ...
    struct('pc', 'dend', 'pv', '',     'fig_num', 20, 'scope', 'pc_only'), ...
    struct('pc', '',     'pv', 'dend', 'fig_num', 21, 'scope', 'pv_only') ...
};

table_out_dir = fullfile(pwd, 'exports_model_tables');
master_tbl = table();

for k = 1:numel(runs)
    r = runs{k};
    run_tag = build_run_tag(r.pc, r.pv, r.scope);
    run_label = build_run_label(r.pc, r.pv, r.scope);
    fprintf('\nRunning targeted kinematics for %s -> Figure %d\n', run_label, r.fig_num);
    swr_plot_pc_pvbc_targeted_kinematics(V, r.pc, r.pv, r.fig_num, r.scope);

    run_csv = fullfile(table_out_dir, sprintf('bluebrain_channel_similarity_residuals_%s.csv', run_tag));
    if ~isfile(run_csv)
        warning('Missing run residual table: %s', run_csv);
        continue;
    end

    t = readtable(run_csv);
    if isempty(t)
        continue;
    end
    t.figure_number = repmat(r.fig_num, height(t), 1);
    t.run_tag = repmat(string(run_tag), height(t), 1);
    t.run_label = repmat(string(run_label), height(t), 1);
    t.same_compartment = strcmpi(string(t.bluebrain_compartment), string(t.other_compartment));
    master_tbl = [master_tbl; t]; %#ok<AGROW>
end

export_master_similarity_table(master_tbl, table_out_dir);
end

function export_master_similarity_table(master_tbl, out_dir)
if isempty(master_tbl)
    fprintf('No per-run residual tables were found; master ranking was not created.\n');
    return;
end

master_tbl = master_tbl(master_tbl.same_compartment, :);
if isempty(master_tbl)
    fprintf('No same-compartment residual rows available for master ranking.\n');
    return;
end

if any(strcmp(master_tbl.Properties.VariableNames, 'rank'))
    master_tbl = removevars(master_tbl, 'rank');
end

% Keep one row per channel/model/compartment pairing across all runs.
keys = strcat( ...
    string(master_tbl.channel), '|', ...
    string(master_tbl.bluebrain_compartment), '|', ...
    string(master_tbl.other_model), '|', ...
    string(master_tbl.other_compartment));
[~, ia] = unique(keys, 'stable');
master_tbl = master_tbl(sort(ia), :);

master_tbl = sortrows(master_tbl, {'residual_rmse','residual_mae'}, {'ascend','ascend'});
master_tbl.master_rank = (1:height(master_tbl)).';
master_tbl = movevars(master_tbl, 'master_rank', 'Before', 'channel');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_csv = fullfile(out_dir, 'bluebrain_channel_similarity_residuals_master_ranked.csv');
writetable(master_tbl, out_csv);

fprintf('\nMaster BlueBrain similarity ranking (same-compartment, lower RMSE is more similar):\n');
disp(master_tbl(:, {'master_rank','channel','bluebrain_compartment','other_model','residual_rmse','similarity','run_tag'}));
fprintf('Saved master residual table: %s\n', out_csv);
end

function out = sanitize_token(txt)
out = lower(strtrim(char(txt)));
out = regexprep(out, '[^a-z0-9]+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'na';
end
end

function run_tag = build_run_tag(pc_section, pv_section, scope)
switch lower(scope)
    case 'both'
        run_tag = sprintf('pc_%s__pv_%s', sanitize_token(pc_section), sanitize_token(pv_section));
    case 'pc_only'
        run_tag = sprintf('pc_%s', sanitize_token(pc_section));
    case 'pv_only'
        run_tag = sprintf('pv_%s', sanitize_token(pv_section));
    otherwise
        error('Unsupported scope: %s', scope);
end
end

function run_label = build_run_label(pc_section, pv_section, scope)
switch lower(scope)
    case 'both'
        run_label = sprintf('PC %s | PVBC %s', pc_section, pv_section);
    case 'pc_only'
        run_label = sprintf('PC %s', pc_section);
    case 'pv_only'
        run_label = sprintf('PVBC %s', pv_section);
    otherwise
        error('Unsupported scope: %s', scope);
end
end
