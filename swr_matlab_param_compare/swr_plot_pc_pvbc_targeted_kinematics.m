function swr_plot_pc_pvbc_targeted_kinematics(V, pc_section, pv_section, complete_fig_num, plot_scope)
% Targeted PC/PVBC channel kinematics across Canakci, BlueBrain, Bezaire.
% - Voltage-dependent gate/open probability comparison
% - Ca_i sweeps only for calcium-dependent channels

if nargin < 1 || isempty(V)
    V = linspace(-100, 50, 600);
end
if nargin < 2 || isempty(pc_section)
    pc_section = 'soma';
end
if nargin < 3 || isempty(pv_section)
    pv_section = 'soma';
end
if nargin < 4 || isempty(complete_fig_num)
    complete_fig_num = 18;
end
if nargin < 5 || isempty(plot_scope)
    plot_scope = 'both';
end

pc_section = normalize_section_name(pc_section);
pv_section = normalize_section_name(pv_section);
plot_scope = lower(strtrim(char(plot_scope)));
switch plot_scope
    case 'both'
        run_tag = sprintf('pc_%s__pv_%s', sanitize_token(pc_section), sanitize_token(pv_section));
        run_label = sprintf('PC %s | PVBC %s', pc_section, pv_section);
        pc_sections = {pc_section};
        pv_sections = {pv_section};
    case 'pc_only'
        run_tag = sprintf('pc_%s', sanitize_token(pc_section));
        run_label = sprintf('PC %s', pc_section);
        pc_sections = {pc_section};
        pv_sections = {};
    case 'pv_only'
        run_tag = sprintf('pv_%s', sanitize_token(pv_section));
        run_label = sprintf('PVBC %s', pv_section);
        pc_sections = {};
        pv_sections = {pv_section};
    otherwise
        error('Unsupported plot_scope: %s (expected both|pc_only|pv_only)', plot_scope);
end

% Canakci kinetics (single mechanism set in this script). We still label
% each compartment explicitly in legends for side-by-side section figures.
g_can_pc = hh_compute_gating(V, hh_canakci_pyr_params());
g_can_pv = hh_compute_gating(V, hh_canakci_pv_params(34));
p_can_pc_na = g_can_pc.m_inf.^3 .* g_can_pc.h_inf .* g_can_pc.i_inf;
p_can_pv_na = g_can_pv.m_inf.^3 .* g_can_pv.h_inf;

% BlueBrain kinetics by compartment from SONATA parameters.
bb_pc_by_sec = load_bluebrain_sections(V, 'pc', pc_sections);
bb_pv_by_sec = load_bluebrain_sections(V, 'pvbc', pv_sections);

% Bezaire kinetics by compartment from mbezaire/ca1 equations in this file.
bz_pc_by_sec = struct();
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    bz_pc_by_sec.(sec) = bezaire_pc_rates_by_section(V, 34, sec);
end

bz_pv_by_sec = struct();
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    bz_pv_by_sec.(sec) = bezaire_pvbc_rates_by_section(V, 34, 5e-6, sec);
end

col_can = [0.00 0.45 0.74];
col_bb  = [0.85 0.33 0.10];
col_bz  = [0.20 0.60 0.20];

out_dir = fullfile(pwd, 'exports_channel_kinematics');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
fig_handles = gobjects(0,1);

% Compute BlueBrain channels beyond na/kdr so missing channels are included.
bb_extra = bluebrain_extra_rates(V, 34, 50e-6);
bb_pc_json = load_bluebrain_json('pc');
bb_pv_json = load_bluebrain_json('pvbc');

% -------------------------------------------------------------------------
% 1) Shared/similar channels first
% -------------------------------------------------------------------------
shared_defs = {};

% PC Na
entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    entries = add_channel_entry(entries, 'Canakci', sec, ...
        {'m','h','i'}, ...
        {g_can_pc.m_inf, g_can_pc.h_inf, g_can_pc.i_inf}, ...
        {g_can_pc.tau_m_ms, g_can_pc.tau_h_ms, g_can_pc.tau_i_ms}, ...
        p_can_pc_na);

    if isfield(bb_pc_by_sec, sec) && section_has_gbar(bb_pc_by_sec.(sec), 'gbar_Na')
        g = bb_pc_by_sec.(sec);
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h'}, ...
            {g.m_inf, g.h_inf}, ...
            {g.tau_m_ms, g.tau_h_ms}, ...
            g.p_open_Na);
    end

    bz = bz_pc_by_sec.(sec);
    if strcmpi(sec, 'axon')
        entries = add_channel_entry(entries, 'Bezaire', sec, ...
            {'m','h'}, ...
            {bz.m_nav, bz.h_nav}, ...
            {bz.tau_m_nav, bz.tau_h_nav}, ...
            bz.p_nav);
    else
        entries = add_channel_entry(entries, 'Bezaire', sec, ...
            {'m','h','s'}, ...
            {bz.m_nav, bz.h_nav, bz.s_nav}, ...
            {bz.tau_m_nav, bz.tau_h_nav, bz.tau_s_nav}, ...
            bz.p_nav);
    end
end
shared_defs = append_channel_def(shared_defs, 'PC Na', entries);

% PC KDR
entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    entries = add_channel_entry(entries, 'Canakci', sec, ...
        {'n'}, {g_can_pc.n_inf}, {g_can_pc.tau_n_ms}, g_can_pc.n_inf);

    if isfield(bb_pc_by_sec, sec) && section_has_gbar(bb_pc_by_sec.(sec), 'gbar_Kdr')
        g = bb_pc_by_sec.(sec);
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n'}, {g.n_inf}, {g.tau_n_ms}, g.p_open_Kdr);
    end

    bz = bz_pc_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'n'}, {bz.n_kdr}, {bz.tau_n_kdr}, bz.p_kdr);
end
shared_defs = append_channel_def(shared_defs, 'PC KDR', entries);

% PVBC Na
entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    entries = add_channel_entry(entries, 'Canakci', sec, ...
        {'m','h'}, ...
        {g_can_pv.m_inf, g_can_pv.h_inf}, ...
        {g_can_pv.tau_m_ms, g_can_pv.tau_h_ms}, ...
        p_can_pv_na);

    if isfield(bb_pv_by_sec, sec) && section_has_gbar(bb_pv_by_sec.(sec), 'gbar_Na')
        g = bb_pv_by_sec.(sec);
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h','s'}, ...
            {g.m_inf, g.h_inf, g.s_inf}, ...
            {g.tau_m_ms, g.tau_h_ms, g.tau_s_ms}, ...
            g.p_open_Na);
    end

    bz = bz_pv_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'m','h'}, ...
        {bz.m_nav, bz.h_nav}, ...
        {bz.tau_m_nav, bz.tau_h_nav}, ...
        bz.p_nav);
end
shared_defs = append_channel_def(shared_defs, 'PVBC Na', entries);

% PVBC KDR
entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    entries = add_channel_entry(entries, 'Canakci', sec, ...
        {'n'}, {g_can_pv.n_inf}, {g_can_pv.tau_n_ms}, g_can_pv.n_inf);

    if isfield(bb_pv_by_sec, sec) && section_has_gbar(bb_pv_by_sec.(sec), 'gbar_Kdrb')
        g = bb_pv_by_sec.(sec);
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n'}, {g.n_inf}, {g.tau_n_ms}, g.p_open_Kdrb);
    end

    bz = bz_pv_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'n'}, {bz.n_kdrfast}, {bz.tau_n_kdrfast}, bz.p_kdrfast);
end
shared_defs = append_channel_def(shared_defs, 'PVBC KDR', entries);

% PC A-type K (family-matched)
entries = struct([]);
bb_pc_kap_secs = pc_sections;
for s = 1:numel(bb_pc_kap_secs)
    sec = bb_pc_kap_secs{s};
    if isfield(bb_pc_by_sec, sec) && json_has_param(bb_pc_json, sec, 'gkabar_kap')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n','l'}, ...
            {bb_extra.pc.kap.n, bb_extra.pc.kap.l}, ...
            {bb_extra.pc.kap.tau_n, bb_extra.pc.kap.tau_l}, ...
            bb_extra.pc.kap.p_open);
    end
end
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    bz = bz_pc_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'n','l'}, {bz.n_kvap, bz.l_kvap}, ...
        {bz.tau_n_kvap, bz.tau_l_kvap}, ...
        bz.p_kvap);
end
shared_defs = append_channel_def(shared_defs, 'PC A-type K', entries);

% PVBC A-type K (family-matched; kap in soma/axon, kad in dend/apic)
entries = struct([]);
bb_pv_a_secs = pv_sections;
for s = 1:numel(bb_pv_a_secs)
    sec = bb_pv_a_secs{s};
    [a_rates, a_tag] = bluebrain_pv_a_rates_for_section(bb_extra, bb_pv_json, sec);
    if ~isempty(a_tag) && isfield(bb_pv_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n','l'}, ...
            {a_rates.n, a_rates.l}, ...
            {a_rates.tau_n, a_rates.tau_l}, ...
            a_rates.p_open);
    end
end
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    bz = bz_pv_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'n','l'}, ...
        {bz.n_kva, bz.l_kva}, ...
        {bz.tau_n_kva, bz.tau_l_kva}, ...
        bz.p_kva);
end
shared_defs = append_channel_def(shared_defs, 'PVBC A-type K', entries);

% Ca-dependent families are reported as cai-vs-probability only (below).

for k = 1:numel(shared_defs)
    fig_handles(end+1) = plot_channel_definition(V, shared_defs{k}, out_dir, 'shared', run_tag, run_label, col_can, col_bb, col_bz); %#ok<AGROW>
end

% -------------------------------------------------------------------------
% 2) Individual model-specific channels after shared ones
% -------------------------------------------------------------------------
ind_defs = {};

% BlueBrain PC-only families
entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    if isfield(bb_pc_by_sec, sec) && json_has_param(bb_pc_json, sec, 'gbar_kmb')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m'}, {bb_extra.pc.kmb.m}, {bb_extra.pc.kmb.tau_m}, bb_extra.pc.kmb.p_open);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PC M-type K', entries);

entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    if isfield(bb_pc_by_sec, sec) && json_has_param(bb_pc_json, sec, 'gcatbar_cat')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h'}, ...
            {bb_extra.pc.cat.m, bb_extra.pc.cat.h}, ...
            {bb_extra.pc.cat.tau_m, bb_extra.pc.cat.tau_h}, ...
            bb_extra.pc.cat.p_open);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PC T-type Ca', entries);

% BlueBrain PVBC-only families
entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    if isfield(bb_pv_by_sec, sec) && section_has_gbar(bb_pv_by_sec.(sec), 'gbar_Kdb')
        g = bb_pv_by_sec.(sec);
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n'}, {g.n_kdb_inf}, {g.tau_n_kdb_ms}, g.p_open_Kdb);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PVBC KDB', entries);

entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    if isfield(bb_pv_by_sec, sec) && json_has_param(bb_pv_json, sec, 'gbar_kmb')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m'}, {bb_extra.pv.kmb.m}, {bb_extra.pv.kmb.tau_m}, bb_extra.pv.kmb.p_open);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PVBC M-type K', entries);

entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    if isfield(bb_pv_by_sec, sec) && json_has_param(bb_pv_json, sec, 'ghdbar_hd')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'l'}, {bb_extra.pv.hd.l}, {bb_extra.pv.hd.tau_l}, bb_extra.pv.hd.p_open);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PVBC HCN', entries);

entries = struct([]);
for s = 1:numel(pv_sections)
    sec = pv_sections{s};
    if isfield(bb_pv_by_sec, sec) && json_has_param(bb_pv_json, sec, 'gcatbar_cat')
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h'}, ...
            {bb_extra.pv.cat.m, bb_extra.pv.cat.h}, ...
            {bb_extra.pv.cat.tau_m, bb_extra.pv.cat.tau_h}, ...
            bb_extra.pv.cat.p_open);
    end
end
ind_defs = append_channel_def(ind_defs, 'BlueBrain PVBC T-type Ca', entries);

% Bezaire PC channels that do not have BlueBrain counterparts in these sections.
entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    if strcmpi(sec, 'soma')
        continue;
    end
    bz = bz_pc_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'n','l'}, ...
        {bz.n_kvad, bz.l_kvad}, ...
        {bz.tau_n_kvad, bz.tau_l_kvad}, ...
        bz.p_kvad);
end
ind_defs = append_channel_def(ind_defs, 'Bezaire PC A-type K distal', entries);

entries = struct([]);
for s = 1:numel(pc_sections)
    sec = pc_sections{s};
    if strcmpi(sec, 'soma')
        continue;
    end
    bz = bz_pc_by_sec.(sec);
    entries = add_channel_entry(entries, 'Bezaire', sec, ...
        {'l'}, ...
        {bz.l_hcn}, ...
        {bz.tau_l_hcn}, ...
        bz.p_hcn);
end
ind_defs = append_channel_def(ind_defs, 'Bezaire PC HCN', entries);

for k = 1:numel(ind_defs)
    fig_handles(end+1) = plot_channel_definition(V, ind_defs{k}, out_dir, 'individual', run_tag, run_label, col_can, col_bb, col_bz); %#ok<AGROW>
end

all_defs = [shared_defs ind_defs];
sim_tbl = build_bluebrain_similarity_table(all_defs);
table_out_dir = fullfile(pwd, 'exports_model_tables');
export_bluebrain_similarity_table(sim_tbl, table_out_dir, run_tag, run_label);

% --------------------------
% Calcium-only sweeps for calcium-dependent channels (cai vs P_open only)
% --------------------------
cai = logspace(log10(1e-6), log10(2e-3), 350); % mM
Vfix = -40;

f_cai = figure('Name', sprintf('Calcium-dependent channels (cai vs probability) [%s]', run_label), 'Color','w');
hold on;
cai_panel = init_cai_panel();
if ~isempty(pv_sections)
    cai_panel.model_names = {'BlueBrain','Bezaire'};
else
    cai_panel.model_names = {'BlueBrain'};
end

if ~isempty(pc_sections)
    [m_bb_cal_fix, ~] = bb_cal_rates(Vfix, 34, 5, 0.2, 0.1, 2, 4, 0.1);
    [m_bb_can_fix, ~, h_bb_can_fix, ~] = bb_can_rates(Vfix, 34, 5, 0.2, 3, 0.03, 2, -14, 0.1);
    p_bb_kca = (cai ./ 0.00035).^4 ./ (1 + (cai ./ 0.00035).^4);
    p_bb_kca = p_bb_kca.^3;
    [p_bb_cagk_cai, ~] = bb_cagk_rates(Vfix .* ones(size(cai)), 34, cai, 0.84, 1, 0.48e-3, 0.13e-6, 0.28, 0.48);
    p_bb_cal_cai = (m_bb_cal_fix.^2) .* (cai ./ (cai + 50e-6));
    p_bb_can_cai = (m_bb_can_fix.^2) .* h_bb_can_fix .* (cai ./ (cai + 50e-6));

    h = plot(cai*1e3, p_bb_kca, 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PC %s kca', pc_section));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_cagk_cai, 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PC %s cagk @ V=%g', pc_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_cal_cai, 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PC %s cal @ V=%g', pc_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_can_cai, 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PC %s can @ V=%g', pc_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
end

if ~isempty(pv_sections)
    % Bezaire PV calcium-dependent channels.
    [p_bz_kcas_cai, ~, p_bz_cavl_h2_cai] = bezaire_ca_sweeps(cai, Vfix, 34);
    [p_bz_kvcab_cai, ~] = kvcab_rates(Vfix .* ones(size(cai)), cai, 34);
    [m_bz_cavl_fix, ~] = cavl_m_rates(Vfix);
    p_bz_cavl_cai = (m_bz_cavl_fix.^2) .* p_bz_cavl_h2_cai;

    % BlueBrain PV uses same mechanism forms as PC in this local set.
    [m_bb_cal_fix, ~] = bb_cal_rates(Vfix, 34, 5, 0.2, 0.1, 2, 4, 0.1);
    [m_bb_can_fix, ~, h_bb_can_fix, ~] = bb_can_rates(Vfix, 34, 5, 0.2, 3, 0.03, 2, -14, 0.1);
    p_bb_kca = (cai ./ 0.00035).^4 ./ (1 + (cai ./ 0.00035).^4);
    p_bb_kca = p_bb_kca.^3;
    [p_bb_cagk_cai, ~] = bb_cagk_rates(Vfix .* ones(size(cai)), 34, cai, 0.84, 1, 0.48e-3, 0.13e-6, 0.28, 0.48);
    p_bb_cal_cai = (m_bb_cal_fix.^2) .* (cai ./ (cai + 50e-6));
    p_bb_can_cai = (m_bb_can_fix.^2) .* h_bb_can_fix .* (cai ./ (cai + 50e-6));

    h = plot(cai*1e3, p_bb_kca, '--', 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PVBC %s kca', pv_section));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_cagk_cai, '--', 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PVBC %s cagk @ V=%g', pv_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_cal_cai, '--', 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PVBC %s cal @ V=%g', pv_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bb_can_cai, '--', 'LineWidth', 1.7, 'DisplayName', sprintf('BlueBrain PVBC %s can @ V=%g', pv_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bz_kcas_cai, ':', 'LineWidth', 1.9, 'DisplayName', sprintf('Bezaire PVBC %s KCaS', pv_section));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bz_kvcab_cai, ':', 'LineWidth', 1.9, 'DisplayName', sprintf('Bezaire PVBC %s KvCaB @ V=%g', pv_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
    h = plot(cai*1e3, p_bz_cavl_cai, ':', 'LineWidth', 1.9, 'DisplayName', sprintf('Bezaire PVBC %s CavL @ V=%g', pv_section, Vfix));
    cai_panel = add_cai_curve(cai_panel, h);
end

set(gca,'XScale','log');
xlabel('cai (uM)');
ylabel('P_{open}');
grid on;
place_legend_safely(gca);
set(gcf,'Position',[120 180 1200 420]);
exportgraphics(f_cai, fullfile(out_dir, sprintf('paper_compare_PC_PVBC_cai_sweeps_%s.png', run_tag)), ...
    'Resolution',180);
fig_handles(end+1) = f_cai; %#ok<AGROW>

f_complete = build_complete_figure(V, all_defs, cai_panel, out_dir, complete_fig_num, run_tag, run_label, col_can, col_bb, col_bz);
if ishghandle(f_complete)
    fig_handles(end+1) = f_complete; %#ok<AGROW>
end

export_target_figures_to_powerpoint(fig_handles, out_dir, ...
    sprintf('paper_compare_PC_PVBC_targeted_kinematics_%s.pptx', run_tag));
end

function by_sec = load_bluebrain_sections(V, cell_type, sections)
by_sec = struct();
for s = 1:numel(sections)
    sec = sections{s};
    try
        switch lower(cell_type)
            case 'pc'
                hh = hh_channel_kinematics_pc_params([], sec);
            case 'pvbc'
                hh = hh_channel_kinematics_pvbc_params([], sec);
            otherwise
                error('Unsupported cell_type: %s', cell_type);
        end
        by_sec.(sec) = hh_compute_gating(V, hh);
    catch ME
        warning('BlueBrain load failed for %s/%s: %s', cell_type, sec, ME.message);
    end
end
end

function raw = load_bluebrain_json(cell_type)
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
switch lower(cell_type)
    case {'pc','pyramidal'}
        json_path = fullfile(repo_root, 'PC_dynamics_params_sonata.json');
    case {'pv','pvbc'}
        json_path = fullfile(repo_root, 'PVBC_dynamics_params_sonata.json');
    otherwise
        error('Unsupported cell_type for JSON load: %s', cell_type);
end

if ~isfile(json_path)
    error('Could not open SONATA dynamics JSON: %s', json_path);
end
raw = jsondecode(fileread(json_path));
end

function tf = json_has_param(raw, section_name, param_name)
if ~isfield(raw, 'genome') || isempty(raw.genome)
    tf = false;
    return;
end

idx = find(arrayfun(@(g) strcmp(g.section, section_name) && strcmp(g.name, param_name), raw.genome), 1, 'first');
if isempty(idx)
    tf = false;
    return;
end

val = raw.genome(idx).value;
tf = ~(isempty(val) || isnan(val) || val <= 0);
end

function [rates, tag] = bluebrain_pv_a_rates_for_section(bb_extra, bb_pv_json, sec)
rates = struct();
tag = '';
if json_has_param(bb_pv_json, sec, 'gkabar_kap')
    rates = bb_extra.pv.kap;
    tag = 'kap';
    return;
end
if json_has_param(bb_pv_json, sec, 'gkabar_kad')
    rates = bb_extra.pv.kad;
    tag = 'kad';
    return;
end
end

function tf = section_has_gbar(g, field_name)
if ~isfield(g, field_name)
    tf = false;
    return;
end
v = g.(field_name);
tf = ~(isempty(v) || all(isnan(v)) || all(v <= 0));
end

function out = bezaire_pc_rates_by_section(V, celsius, section_name)
out = bezaire_pc_rates(V, celsius);
if strcmpi(section_name, 'axon')
    % Axonal sodium in Bezaire PC uses ch_Navaxonp kinetics (m,h).
    [m,h,tm,th] = navaxon_rates(V, celsius, 15);
    out.m_nav = m;
    out.h_nav = h;
    out.s_nav = nan(size(V));
    out.tau_m_nav = tm;
    out.tau_h_nav = th;
    out.tau_s_nav = nan(size(V));
    out.p_nav = m.^3 .* h;
end
end

function out = bezaire_pvbc_rates_by_section(V, celsius, cai, ~)
out = bezaire_pvbc_rates(V, celsius, cai);
end

function out = bezaire_pc_rates(V, celsius)
[m,h,s,tm,th,ts] = navp_rates(V, celsius, 15, 1);
[n_kdr, t_kdr] = kdrp_rates(V, celsius);
[n_kvap,l_kvap,t_n_kvap,t_l_kvap] = kva_prox_rates(V, celsius);
[n_kvad,l_kvad,t_n_kvad,t_l_kvad] = kva_dist_rates(V, celsius);
[l_hcn, t_l_hcn] = hcn_rates(V, celsius, -82);

out.m_nav = m;
out.h_nav = h;
out.s_nav = s;
out.tau_m_nav = tm;
out.tau_h_nav = th;
out.tau_s_nav = ts;

out.n_kdr = n_kdr;
out.tau_n_kdr = t_kdr;

out.n_kvap = n_kvap;
out.l_kvap = l_kvap;
out.tau_n_kvap = t_n_kvap;
out.tau_l_kvap = t_l_kvap;

out.n_kvad = n_kvad;
out.l_kvad = l_kvad;
out.tau_n_kvad = t_n_kvad;
out.tau_l_kvad = t_l_kvad;

out.l_hcn = l_hcn;
out.tau_l_hcn = t_l_hcn;

out.p_nav = m.^3 .* h .* s;
out.p_kdr = n_kdr;
out.p_kvap = n_kvap .* l_kvap;
out.p_kvad = n_kvad .* l_kvad;
out.p_hcn = l_hcn;
end

function out = bezaire_pvbc_rates(V, celsius, cai)
[m,h,tm,th] = navaxon_rates(V, celsius, 15);
[n_kdrfast,t_kdrfast] = kdrfast_rates(V, celsius);
[n_kva,l_kva,t_n_kva,t_l_kva] = kva_fast_rates(V, celsius);
[c,d,t_c,t_d] = cavn_rates(V, celsius);
[m_l,t_m_l] = cavl_m_rates(V);
h2_l = 0.001 ./ (0.001 + cai);
[q,t_q] = kcas_rates(cai, celsius);
[o,t_o] = kvcab_rates(V, cai, celsius);

out.m_nav = m;
out.h_nav = h;
out.tau_m_nav = tm;
out.tau_h_nav = th;

out.n_kdrfast = n_kdrfast;
out.tau_n_kdrfast = t_kdrfast;

out.n_kva = n_kva;
out.l_kva = l_kva;
out.tau_n_kva = t_n_kva;
out.tau_l_kva = t_l_kva;

out.c_cavn = c;
out.d_cavn = d;
out.tau_c_cavn = t_c;
out.tau_d_cavn = t_d;

out.m_cavl = m_l;
out.tau_m_cavl = t_m_l;

out.q_kcas = q;
out.tau_q_kcas = t_q;

out.o_kvcab = o;
out.tau_o_kvcab = t_o;

out.p_nav = m.^3 .* h;
out.p_kdrfast = n_kdrfast.^4;
out.p_kva = n_kva .* l_kva;
out.p_cavn = c.^2 .* d;
out.p_cavl = m_l.^2 .* h2_l;
out.p_kcas = q.^2;
out.p_kvcab = o;
end

function [pkcas, pkvcab, ph2] = bezaire_ca_sweeps(cai, Vfix, celsius)
[q,~] = kcas_rates(cai, celsius);
pkcas = q.^2;
[pkvcab,~] = kvcab_rates(Vfix .* ones(size(cai)), cai, celsius);
ph2 = 0.001 ./ (0.001 + cai);
end

function [m,h,s,tau_m,tau_h,tau_s] = navp_rates(V, celsius, sh, ar2)
% Bezaire ch_Navp.mod (m,h,s gates)
qt = 2.^((celsius-24)./10);

a = trap0(V, -30+sh, 0.4, 7.2);
b = trap0(-V, 30-sh, 0.124, 7.2);
m = a./(a+b);
tau_m = max(1./(a+b)./qt, 0.02);

a = trap0(V, -45+sh, 0.03, 1.5);
b = trap0(-V, 45-sh, 0.01, 1.5);
h = 1./(1+exp((V+50-sh)./4));
tau_h = max(1./(a+b)./qt, 0.5);

alpv = 1./(1+exp((V+58-sh)./2));
alps = exp(1e-3.*12.*(V+60-sh).*9.648e4./(8.315.*(273.16+celsius)));
bets = exp(1e-3.*12.*0.2.*(V+60-sh).*9.648e4./(8.315.*(273.16+celsius)));
s = alpv + ar2.*(1-alpv);
tau_s = max(bets./(0.0003.*(1+alps)), 10);
end

function [m,h,tau_m,tau_h] = navaxon_rates(V, celsius, sh)
% Bezaire ch_Navaxonp.mod (m,h gates)
qt = 2.^((celsius-24)./10);

a = trap0(V, -30+sh, 0.4, 7.2);
b = trap0(-V, 30-sh, 0.124, 7.2);
m = a./(a+b);
tau_m = max(1./(a+b)./qt, 0.02);

a = trap0(V, -45+sh, 0.03, 1.5);
b = trap0(-V, 45-sh, 0.01, 1.5);
h = 1./(1+exp((V+50-sh)./4));
tau_h = max(1./(a+b)./qt, 0.5);
end

function [n, tau_n] = kdrp_rates(V, celsius)
% Bezaire ch_Kdrp.mod (n gate)
a = exp(1e-3*(-3).*(V-13).*9.648e4./(8.315.*(273.16+celsius)));
b = exp(1e-3*(-3)*0.7.*(V-13).*9.648e4./(8.315.*(273.16+celsius)));
n = 1./(1+a);
tau_n = max(b./(0.02.*(1+a)), 2);
end

function [n,l,tau_n,tau_l] = kva_prox_rates(V, celsius)
% Bezaire ch_KvAproxp.mod
zeta = -1.5 + (-1)./(1+exp((V+40)./5));
an = exp(1e-3.*zeta.*(V-11).*9.648e4./(8.315.*(273.16+celsius)));
bn = exp(1e-3.*zeta.*0.55.*(V-11).*9.648e4./(8.315.*(273.16+celsius)));
al = exp(1e-3.*3.*(V+56).*9.648e4./(8.315.*(273.16+celsius)));
n = 1./(1+an);
l = 1./(1+al);
qt = 5.^((celsius-24)./10);
tau_n = max(bn./(qt.*0.05.*(1+an)), 0.1);
tau_l = max(0.26.*(V+50), 2);
end

function [n,l,tau_n,tau_l] = kva_dist_rates(V, celsius)
% Bezaire ch_KvAdistp.mod
zeta = -1.8 + (-1)./(1+exp((V+40)./5));
an = exp(1e-3.*zeta.*(V+1).*9.648e4./(8.315.*(273.16+celsius)));
bn = exp(1e-3.*zeta.*0.39.*(V+1).*9.648e4./(8.315.*(273.16+celsius)));
al = exp(1e-3.*3.*(V+56).*9.648e4./(8.315.*(273.16+celsius)));
n = 1./(1+an);
l = 1./(1+al);
qt = 5.^((celsius-24)./10);
tau_n = max(bn./(qt.*0.1.*(1+an)), 0.2);
tau_l = max(0.26.*(V+50), 2);
end

function [l, tau_l] = hcn_rates(V, celsius, vhalfl)
% Bezaire ch_HCNp.mod
alpl = exp(0.0378*4.*(V-vhalfl));
alpt = exp(0.0378*2.2.*(V+75));
bett = exp(0.0378*2.2*0.4.*(V+75));
qt = 4.5.^((celsius-33)./10);
l = 1./(1 + alpl);
tau_l = bett./(qt.*0.011.*(1+alpt));
end

function [n, tau_n] = kdrfast_rates(V, celsius)
% Bezaire ch_Kdrfast.mod
alpha = -0.07.*vtrap(V+18, -6);
beta = 0.264./exp((V+43)./40);
n = alpha./(alpha+beta);
q10 = 3.^((celsius-34)./10);
tau_n = (1./(alpha+beta)) ./ q10;
end

function [n,l,tau_n,tau_l] = kva_fast_rates(V, celsius)
% Bezaire ch_KvA.mod
an = exp(1e-3*(-3).*(V+33.6).*9.648e4./(8.315.*(273.16+celsius)));
bn = exp(1e-3*(-3)*0.6.*(V+33.6).*9.648e4./(8.315.*(273.16+celsius)));
al = exp(1e-3*(4).*(V+83).*9.648e4./(8.315.*(273.16+celsius)));
bl = exp(1e-3*(4)*1.*(V+83).*9.648e4./(8.315.*(273.16+celsius)));
n = 1./(1+an);
l = 1./(1+al);
q10 = 3.^((celsius-30)./10);
tau_n = bn./(q10.*0.02.*(1+an));
tau_l = bl./(q10.*0.08.*(1+al));
end

function [c,d,tau_c,tau_d] = cavn_rates(V, celsius)
% Bezaire ch_CavN.mod
ac = -0.19.*vtrap(V-19.88, -10);
bc = 0.046.*exp(-V./20.73);
c = ac./(ac+bc);
ad = 1.6e-4.*exp(-V./48.4);
bd = 1./(exp((-V+39)./10)+1);
d = ad./(ad+bd);
q10 = 3.^((celsius-34)./10);
tau_c = (1./(ac+bc)) ./ q10;
tau_d = (1./(ad+bd)) ./ q10;
end

function [m, tau_m] = cavl_m_rates(V)
% Bezaire ch_CavL.mod
a = 15.69.*(-V+81.5)./(exp((-V+81.5)./10)-1);
b = 0.29.*exp(-V./10.86);
m = a./(a+b);
tau_m = 1./(a+b);
end

function [q, tau_q] = kcas_rates(cai, celsius)
% Bezaire ch_KCaS.mod (calcium-only gate)
alpha = 15 .* cai.^2;
beta = 2.5e-4;
q = alpha./(alpha+beta);
q10 = 3.^((celsius-34)./10);
tau_q = 1./(alpha+beta)./q10;
end

function [o, tau_o] = kvcab_rates(V, cai, celsius)
% Bezaire ch_KvCaB.mod
exp1a = 0.48e-3 .* exp(-2*0.84*96485.*V./(8.313424.*(273.15+celsius)));
exp1b = 0.13e-6 .* exp(-2*1.0*96485.*V./(8.313424.*(273.15+celsius)));
alp = cai.*0.28 ./ (cai + exp1a);
bet = 0.48 ./ (1 + cai./exp1b);
o = alp ./ (alp + bet);
tau_o = 1./(alp + bet);
end

function d = make_channel_def(name, entries)
d = struct();
d.name = name;
d.entries = entries;
end

function entries = add_channel_entry(entries, model_name, compartment, gate_names, gate_open, gate_tau, p_open)
e = struct();
e.model = model_name;
e.compartment = compartment;
e.gate_names = gate_names;
e.gate_open = gate_open;
e.gate_tau = gate_tau;
e.p_open = p_open;

if isempty(entries)
    entries = e;
else
    entries(end+1) = e; %#ok<AGROW>
end
end

function defs = append_channel_def(defs, name, entries)
if isempty(entries)
    return;
end
defs{end+1} = make_channel_def(name, entries); %#ok<AGROW>
end

function f = plot_channel_definition(V, d, out_dir, prefix, run_tag, run_label, col_can, col_bb, col_bz)
gate_styles = {'-','--',':','-.'};
popen_styles = {'-','--',':','-.'};
f = figure('Name', sprintf('Channel Comparison - %s [%s]', d.name, run_label), 'Color', 'w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% opening
nexttile; hold on;
for m = 1:numel(d.entries)
    e = d.entries(m);
    c = model_color(e.model, col_can, col_bb, col_bz);
    for g = 1:numel(e.gate_open)
        ls = gate_styles{mod(g-1, numel(gate_styles))+1};
        plot(V, e.gate_open{g}, ls, 'Color', c, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s %s', e.model, e.gate_names{g}));
    end
end
xlabel('V (mV)'); ylabel('x_{\infty}');
grid on; place_legend_safely(gca);

% closing
nexttile; hold on;
for m = 1:numel(d.entries)
    e = d.entries(m);
    c = model_color(e.model, col_can, col_bb, col_bz);
    for g = 1:numel(e.gate_open)
        ls = gate_styles{mod(g-1, numel(gate_styles))+1};
        plot(V, 1 - e.gate_open{g}, ls, 'Color', c, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s 1-%s', e.model, e.gate_names{g}));
    end
end
xlabel('V (mV)'); ylabel('1-x_{\infty}');
grid on; place_legend_safely(gca);

% tau
nexttile; hold on;
for m = 1:numel(d.entries)
    e = d.entries(m);
    c = model_color(e.model, col_can, col_bb, col_bz);
    for g = 1:numel(e.gate_tau)
        if all(isnan(e.gate_tau{g}))
            continue;
        end
        ls = gate_styles{mod(g-1, numel(gate_styles))+1};
        plot(V, e.gate_tau{g}, ls, 'Color', c, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s tau_%s', e.model, e.gate_names{g}));
    end
end
xlabel('V (mV)'); ylabel('\tau (ms)');
grid on; place_legend_safely(gca);

% opening probability
nexttile; hold on;
for m = 1:numel(d.entries)
    e = d.entries(m);
    c = model_color(e.model, col_can, col_bb, col_bz);
    ls = popen_styles{mod(m-1, numel(popen_styles))+1};
    plot(V, e.p_open, ls, 'Color', c, 'LineWidth', 1.8, ...
        'DisplayName', e.model);
end
xlabel('V (mV)'); ylabel('P_{open}');
grid on; place_legend_safely(gca);

set(f, 'Position', [130 80 1300 760]);

safe = regexprep([prefix '_' run_tag '_' d.name], '[^A-Za-z0-9_\- ]', '_');
safe = strrep(safe, ' ', '_');
exportgraphics(f, fullfile(out_dir, ['paper_compare_channel_' safe '.png']), ...
    'Resolution', 180);
end

function place_legend_safely(ax)
lgd = legend(ax, 'Location', 'best', 'Interpreter', 'none');
lgd.Box = 'off';
lgd.Color = 'none';

line_count = numel(findobj(ax, 'Type', 'line'));
entry_count = numel(lgd.String);
if line_count >= 7 || entry_count >= 7
    try
        lgd.Location = 'bestoutside';
    catch
        lgd.Location = 'best';
    end
end
end

function export_target_figures_to_powerpoint(fig_handles, out_dir, pptx_name)
if isempty(fig_handles)
    return;
end

if exist('mlreportgen.ppt.Presentation', 'class') ~= 8
    warning(['MATLAB Report Generator is not available. ', ...
        'Cannot export %s using mlreportgen.ppt workflow.'], pptx_name);
    return;
end

import mlreportgen.ppt.*

ppt_path = fullfile(out_dir, pptx_name);
if exist(ppt_path, 'file')
    try
        delete(ppt_path);
    catch ME
        warning('Could not overwrite PPT file %s: %s', ppt_path, ME.message);
        return;
    end
end

ppt = Presentation(ppt_path);
try
    open(ppt);
catch ME
    warning('Could not open PPT file %s for writing: %s', ppt_path, ME.message);
    return;
end

slide_count = 0;
for k = 1:numel(fig_handles)
    f = fig_handles(k);
    if ~ishandle(f) || ~strcmp(get(f, 'Type'), 'figure')
        continue;
    end

    fig_name = get(f, 'Name');
    if isempty(strtrim(fig_name))
        fig_name = sprintf('Figure %d', k);
    end
    safe_name = regexprep(fig_name, '[^A-Za-z0-9_\- ]', '_');
    if isempty(strtrim(safe_name))
        safe_name = sprintf('figure_%02d', k);
    end
    asset_path = export_figure_asset_for_ppt(f, out_dir, k, safe_name);

    slide = add(ppt, 'Title and Content');
    replace(slide, 'Title', fig_name);
    replace(slide, 'Content', Picture(asset_path));
    slide_count = slide_count + 1;

    % Add one caption slide immediately after the complete-figure slide.
    if contains(lower(fig_name), 'complete targeted channel set')
        [cap_title, cap_body] = supplementary_caption_for_complete_figure(fig_name);
        if ~isempty(cap_title)
            slide = add(ppt, 'Title and Content');
            replace(slide, 'Title', cap_title);
            replace(slide, 'Content', cap_body);
            slide_count = slide_count + 1;
        end
    end
end

close(ppt);
fprintf('PowerPoint saved with %d slides: %s\n', slide_count, ppt_path);
end

function asset_path = export_figure_asset_for_ppt(fig_handle, out_dir, idx, safe_name)
% Prefer vector export for sharper PowerPoint rendering on Windows.
asset_path = '';

if ispc
    emf_path = fullfile(out_dir, sprintf('ppt_%02d_%s.emf', idx, safe_name));
    try
        print(fig_handle, emf_path, '-dmeta', '-painters');
        if exist(emf_path, 'file')
            asset_path = emf_path;
            return;
        end
    catch
    end
end

png_path = fullfile(out_dir, sprintf('ppt_%02d_%s.png', idx, safe_name));
try
    exportgraphics_highres(fig_handle, png_path, 600);
catch
    print(fig_handle, png_path, '-dpng', '-r600');
end
asset_path = png_path;
end

function [cap_title, cap_body] = supplementary_caption_for_complete_figure(fig_name)
cap_title = '';
cap_body = '';

tok = regexp(fig_name, 'Figure\s+(\d+)\s*-\s*Complete targeted channel set\s*\[(.*?)\]\s*$', 'tokens', 'once');
if isempty(tok)
    return;
end

fig_num = str2double(tok{1});
run_label = strtrim(tok{2});
if isnan(fig_num)
    return;
end

% Requested numbering: Figure 18 -> Supplementary Figure 1, etc.
sup_num = fig_num - 17;
if sup_num < 1
    sup_num = fig_num;
end

cap_title = sprintf('Supplementary Figure %d Caption', sup_num);
cap_body = sprintf([ ...
    'Supplementary Figure %d. Complete channel-kinematics comparison for %s.\n' ...
    'Composite comparison of model-derived CA1 channel kinetics for the selected compartment set. ' ...
    'A. Lettered panels (A., B., C., ...) each correspond to one channel family and report gate steady-state variables (x_inf), complementary closing variables (1-x_inf), time constants (tau), and channel opening probability (P_open) across membrane voltage. ' ...
    'Panel headers indicate the channel identity and which models are included (Canakci, BlueBrain, Bezaire, where available).\n' ...
    'B. Calcium-dependent channels are summarized separately as calcium concentration versus opening probability (cai vs P_open), aligned to the same compartment set. ' ...
    'Residual-based BlueBrain-versus-other similarity rankings for matched compartments are provided in the associated supplementary residual tables.'], ...
    sup_num, run_label);
end

function sim_tbl = build_bluebrain_similarity_table(defs)
rows = struct( ...
    'channel', {}, ...
    'bluebrain_compartment', {}, ...
    'other_model', {}, ...
    'other_compartment', {}, ...
    'residual_rmse', {}, ...
    'residual_mae', {}, ...
    'similarity', {});

for i = 1:numel(defs)
    d = defs{i};
    if ~isfield(d, 'entries') || isempty(d.entries)
        continue;
    end
    entries = d.entries;
    model_names = arrayfun(@(e) e.model, entries, 'UniformOutput', false);
    bb_idx = find(strcmpi(model_names, 'BlueBrain'));
    if isempty(bb_idx)
        continue;
    end
    other_idx = find(~strcmpi(model_names, 'BlueBrain'));
    if isempty(other_idx)
        continue;
    end

    for j = 1:numel(other_idx)
        e_other = entries(other_idx(j));
        bb_match = pick_bluebrain_entry(entries, bb_idx, e_other.compartment);
        [rmse, mae] = residual_metrics(bb_match.p_open, e_other.p_open);
        if isnan(rmse)
            continue;
        end

        r = struct();
        r.channel = d.name;
        r.bluebrain_compartment = bb_match.compartment;
        r.other_model = e_other.model;
        r.other_compartment = e_other.compartment;
        r.residual_rmse = rmse;
        r.residual_mae = mae;
        r.similarity = residual_similarity_label(rmse);
        rows(end+1,1) = r; %#ok<AGROW>
    end
end

if isempty(rows)
    sim_tbl = table();
    return;
end

sim_tbl = struct2table(rows);
sim_tbl = sortrows(sim_tbl, {'residual_rmse','residual_mae'}, {'ascend','ascend'});
sim_tbl.rank = (1:height(sim_tbl)).';
sim_tbl = movevars(sim_tbl, 'rank', 'Before', 'channel');
end

function bb_entry = pick_bluebrain_entry(entries, bb_idx, target_compartment)
if nargin < 3
    target_compartment = '';
end

bb_entry = entries(bb_idx(1));
if isempty(target_compartment)
    return;
end

bb_compartments = arrayfun(@(e) e.compartment, entries(bb_idx), 'UniformOutput', false);
match_local = find(strcmpi(bb_compartments, target_compartment), 1, 'first');
if ~isempty(match_local)
    bb_entry = entries(bb_idx(match_local));
end
end

function [rmse, mae] = residual_metrics(a, b)
n = min(numel(a), numel(b));
if n == 0
    rmse = nan;
    mae = nan;
    return;
end

a = a(1:n);
b = b(1:n);
mask = ~(isnan(a) | isnan(b));
if ~any(mask)
    rmse = nan;
    mae = nan;
    return;
end

d = a(mask) - b(mask);
rmse = sqrt(mean(d.^2));
mae = mean(abs(d));
end

function label = residual_similarity_label(rmse)
if rmse < 0.05
    label = 'similar';
elseif rmse < 0.15
    label = 'moderately different';
else
    label = 'different';
end
end

function export_bluebrain_similarity_table(sim_tbl, out_dir, run_tag, run_label)
if isempty(sim_tbl)
    fprintf('No BlueBrain-to-other channel overlaps were found for residual ranking [%s].\n', run_label);
    return;
end

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_csv = fullfile(out_dir, sprintf('bluebrain_channel_similarity_residuals_%s.csv', run_tag));
writetable(sim_tbl, out_csv);

fprintf('\nBlueBrain vs other models residual ranking (P_open RMSE; lower is more similar) [%s]:\n', run_label);
disp(sim_tbl(:, {'rank','channel','other_model','other_compartment','residual_rmse','similarity'}));
fprintf('Saved residual table: %s\n', out_csv);
end

function f_complete = build_complete_figure(V, defs, cai_panel, out_dir, complete_fig_num, run_tag, run_label, col_can, col_bb, col_bz)
if nargin < 5 || isempty(complete_fig_num)
    complete_fig_num = 18;
end
if nargin < 6 || isempty(run_tag)
    run_tag = 'pc_soma__pv_soma';
end
if nargin < 7 || isempty(run_label)
    run_label = 'PC soma | PVBC soma';
end
if nargin < 8 || isempty(col_can)
    col_can = [0.00 0.45 0.74];
end
if nargin < 9 || isempty(col_bb)
    col_bb  = [0.85 0.33 0.10];
end
if nargin < 10 || isempty(col_bz)
    col_bz  = [0.20 0.60 0.20];
end

has_cai = isstruct(cai_panel) && isfield(cai_panel, 'curves') && ~isempty(cai_panel.curves);
n_panels = numel(defs) + double(has_cai);
if n_panels == 0
    f_complete = gobjects(0,1);
    return;
end

final_export_dpi = 900;
f_complete = figure('Visible', 'on'); clf(f_complete);
set(f_complete, 'Name', sprintf('Figure %d - Complete targeted channel set [%s]', complete_fig_num, run_label), ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [20 20 5600 3400]);

n_cols = ceil(sqrt(n_panels));
n_rows = ceil(n_panels / n_cols);
tl = tiledlayout(f_complete, n_rows, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
tl.Padding = 'compact';
tl.TileSpacing = 'compact';
set(tl, 'Units', 'normalized');
tl.Position = [0.01 0.01 0.98 0.98];
drawnow;

panel_idx = 0;
for k = 1:numel(defs)
    d = defs{k};
    panel_idx = panel_idx + 1;
    host_ax = nexttile(tl);
    tile_pos = host_ax.Position;
    configure_placeholder_tile(host_ax);

    model_names = unique(arrayfun(@(e) e.model, d.entries, 'UniformOutput', false), 'stable');
    panel_title = sprintf('%s. %s | Models: %s', ...
        figure18_panel_letter(panel_idx), d.name, strjoin(model_names, ', '));
    draw_channel_panel_in_tile(f_complete, tile_pos, V, d, panel_title, col_can, col_bb, col_bz);
end

if has_cai
    panel_idx = panel_idx + 1;
    host_ax = nexttile(tl);
    tile_pos = host_ax.Position;
    configure_placeholder_tile(host_ax);

    model_names = {'BlueBrain'};
    if isfield(cai_panel, 'model_names') && ~isempty(cai_panel.model_names)
        model_names = cai_panel.model_names;
    end
    panel_title = sprintf('%s. Calcium-dependent channels (cai vs probability) | Models: %s', ...
        figure18_panel_letter(panel_idx), strjoin(model_names, ', '));
    draw_cai_panel_in_tile(f_complete, tile_pos, cai_panel, panel_title);
end

out_png = fullfile(out_dir, sprintf('paper_compare_complete_figure_%02d_%s.png', complete_fig_num, run_tag));
exportgraphics_highres(f_complete, out_png, final_export_dpi);
end

function configure_placeholder_tile(ax)
% Keep tile occupancy so nexttile advances correctly while hiding the host axes.
set(ax, 'Visible', 'off');
ax.XTick = [];
ax.YTick = [];
ax.Box = 'off';
end

function draw_channel_panel_in_tile(fig_handle, tile_pos, V, d, panel_title, col_can, col_bb, col_bz)
title_h = 0.11 * tile_pos(4);
mx = 0.05 * tile_pos(3);
my = 0.08 * tile_pos(4);
gx = 0.09 * tile_pos(3);
gy = 0.18 * tile_pos(4);

plot_box = [tile_pos(1) + mx, tile_pos(2) + my, tile_pos(3) - 2*mx, tile_pos(4) - title_h - my];
if plot_box(3) <= 0 || plot_box(4) <= 0
    return;
end

w = (plot_box(3) - gx) / 2;
h = (plot_box(4) - gy) / 2;
if w <= 0 || h <= 0
    return;
end

p1 = [plot_box(1),          plot_box(2) + h + gy, w, h];
p2 = [plot_box(1) + w + gx, plot_box(2) + h + gy, w, h];
p3 = [plot_box(1),          plot_box(2),          w, h];
p4 = [plot_box(1) + w + gx, plot_box(2),          w, h];

ax1 = axes('Parent', fig_handle, 'Position', p1); hold(ax1, 'on');
ax2 = axes('Parent', fig_handle, 'Position', p2); hold(ax2, 'on');
ax3 = axes('Parent', fig_handle, 'Position', p3); hold(ax3, 'on');
ax4 = axes('Parent', fig_handle, 'Position', p4); hold(ax4, 'on');

gate_styles = {'-','--',':','-.'};
popen_styles = {'-','--',':','-.'};
for m = 1:numel(d.entries)
    e = d.entries(m);
    c = model_color(e.model, col_can, col_bb, col_bz);
    model_lbl = compact_model_label(e.model);

    for g = 1:numel(e.gate_open)
        ls = gate_styles{mod(g-1, numel(gate_styles))+1};
        plot(ax1, V, e.gate_open{g}, ls, 'Color', c, 'LineWidth', 0.9, ...
            'DisplayName', sprintf('%s %s', model_lbl, e.gate_names{g}));
        plot(ax2, V, 1 - e.gate_open{g}, ls, 'Color', c, 'LineWidth', 0.9, ...
            'DisplayName', sprintf('%s 1-%s', model_lbl, e.gate_names{g}));
    end

    for g = 1:numel(e.gate_tau)
        if all(isnan(e.gate_tau{g}))
            continue;
        end
        ls = gate_styles{mod(g-1, numel(gate_styles))+1};
        plot(ax3, V, e.gate_tau{g}, ls, 'Color', c, 'LineWidth', 0.9, ...
            'DisplayName', sprintf('%s tau_%s', model_lbl, e.gate_names{g}));
    end

    ls = popen_styles{mod(m-1, numel(popen_styles))+1};
    plot(ax4, V, e.p_open, ls, 'Color', c, 'LineWidth', 1.0, ...
        'DisplayName', model_lbl);
end

ylabel(ax1, 'x_{\infty}', 'FontSize', 7); title(ax1, 'Open', 'FontSize', 8, 'FontWeight', 'bold');
title(ax2, 'Closed', 'FontSize', 8, 'FontWeight', 'bold');
xlabel(ax3, 'V (mV)', 'FontSize', 7); ylabel(ax3, '\tau (ms)', 'FontSize', 7); title(ax3, '\tau', 'FontSize', 8, 'FontWeight', 'bold');
xlabel(ax4, 'V (mV)', 'FontSize', 7); title(ax4, 'P_{open}', 'FontSize', 8, 'FontWeight', 'bold');
set([ax1 ax2 ax3 ax4], 'FontSize', 7, 'Box', 'on', 'LineWidth', 0.75, ...
    'TickDir', 'out', 'TickLength', [0.012 0.012], 'Layer', 'top');

% Avoid duplicated labels that crowd the 2x2 mini-grid.
ax1.XTickLabel = [];
ax2.XTickLabel = [];
ax2.YTickLabel = [];
ax4.YTickLabel = [];

grid(ax1, 'on'); grid(ax2, 'on'); grid(ax3, 'on'); grid(ax4, 'on');
style_complete_axis(ax1);
style_complete_axis(ax2);
style_complete_axis(ax3);
style_complete_axis(ax4);

pad_axis_limits(ax1);
pad_axis_limits(ax2);
pad_axis_limits(ax3);
pad_axis_limits(ax4);
set_compact_tick_density(ax1, 4, 4);
set_compact_tick_density(ax2, 4, 4);
set_compact_tick_density(ax3, 4, 4);
set_compact_tick_density(ax4, 4, 4);
place_compact_legend(ax1);
place_compact_legend(ax2);
place_compact_legend(ax3);
place_compact_legend(ax4);

annotation(fig_handle, 'textbox', ...
    [tile_pos(1) + 0.002, tile_pos(2) + tile_pos(4) - title_h + 0.002, tile_pos(3) - 0.004, title_h - 0.003], ...
    'String', panel_title, ...
    'Interpreter', 'none', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 9, ...
    'FontWeight', 'bold');
end

function draw_cai_panel_in_tile(fig_handle, tile_pos, cai_panel, panel_title)
title_h = 0.11 * tile_pos(4);
mx = 0.06 * tile_pos(3);
my = 0.10 * tile_pos(4);
ax_pos = [tile_pos(1) + mx, tile_pos(2) + my, tile_pos(3) - 2*mx, tile_pos(4) - title_h - my];
if ax_pos(3) <= 0 || ax_pos(4) <= 0
    return;
end

ax = axes('Parent', fig_handle, 'Position', ax_pos); hold(ax, 'on');
for i = 1:numel(cai_panel.curves)
    c = cai_panel.curves(i);
    plot(ax, c.x, c.y, ...
        'LineStyle', c.line_style, ...
        'Color', c.color, ...
        'LineWidth', max(0.9, c.line_width), ...
        'DisplayName', compact_legend_label(c.label));
end
set(ax, 'XScale', 'log', 'FontSize', 7, 'Box', 'on', ...
    'LineWidth', 0.75, 'TickDir', 'out', 'TickLength', [0.012 0.012], 'Layer', 'top');
xlabel(ax, 'cai (uM)', 'FontSize', 7);
ylabel(ax, 'P_{open}', 'FontSize', 7);
grid(ax, 'on');
style_complete_axis(ax);
pad_axis_limits(ax);
set_compact_tick_density(ax, 4, 4);
place_compact_legend(ax);

annotation(fig_handle, 'textbox', ...
    [tile_pos(1) + 0.002, tile_pos(2) + tile_pos(4) - title_h + 0.002, tile_pos(3) - 0.004, title_h - 0.003], ...
    'String', panel_title, ...
    'Interpreter', 'none', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 9, ...
    'FontWeight', 'bold');
end

function panel = init_cai_panel()
panel = struct();
panel.model_names = {};
panel.curves = struct('x', {}, 'y', {}, 'label', {}, 'line_style', {}, 'color', {}, 'line_width', {});
end

function panel = add_cai_curve(panel, line_handle)
if isempty(line_handle) || ~ishghandle(line_handle)
    return;
end
c = struct();
c.x = get(line_handle, 'XData');
c.y = get(line_handle, 'YData');
c.label = get(line_handle, 'DisplayName');
c.line_style = get(line_handle, 'LineStyle');
c.color = get(line_handle, 'Color');
c.line_width = get(line_handle, 'LineWidth');
if isempty(panel.curves)
    panel.curves = c;
else
    panel.curves(end+1) = c; %#ok<AGROW>
end
end

function place_compact_legend(ax)
line_handles = findobj(ax, 'Type', 'Line');
if isempty(line_handles)
    return;
end

lgd = legend(ax, 'show', 'Location', 'best', 'Interpreter', 'none');
lgd.Box = 'off';
lgd.Color = 'none';

entry_count = numel(lgd.String);
if entry_count <= 4
    lgd.FontSize = 6.5;
elseif entry_count <= 8
    lgd.FontSize = 6;
else
    lgd.FontSize = 5;
end
if entry_count >= 8
    lgd.NumColumns = 2;
else
    lgd.NumColumns = 1;
end
try
    lgd.ItemTokenSize = [8 8];
catch
end

% Force legend to remain inside the parent axis bounds.
ax_units = ax.Units;
lgd_units = lgd.Units;
ax.Units = 'normalized';
lgd.Units = 'normalized';
ap = ax.Position;
lp = lgd.Position;
margin = 0.004;
lp(1) = min(max(lp(1), ap(1) + margin), ap(1) + ap(3) - lp(3) - margin);
lp(2) = min(max(lp(2), ap(2) + margin), ap(2) + ap(4) - lp(4) - margin);
lgd.Position = lp;
ax.Units = ax_units;
lgd.Units = lgd_units;
lgd.AutoUpdate = 'off';
end

function pad_axis_limits(ax)
if ~ishghandle(ax)
    return;
end

xl = xlim(ax);
yl = ylim(ax);
if any(~isfinite(xl)) || any(~isfinite(yl))
    return;
end

if strcmpi(ax.XScale, 'log')
    lo = max(realmin, xl(1) * 0.92);
    hi = xl(2) * 1.08;
    xlim(ax, [lo, hi]);
else
    dx = max(1e-9, xl(2) - xl(1));
    xlim(ax, [xl(1) - 0.03*dx, xl(2) + 0.03*dx]);
end

if strcmpi(ax.YScale, 'log')
    lo = max(realmin, yl(1) * 0.90);
    hi = yl(2) * 1.10;
    ylim(ax, [lo, hi]);
else
    dy = max(1e-9, yl(2) - yl(1));
    ylim(ax, [yl(1) - 0.06*dy, yl(2) + 0.06*dy]);
end
end

function style_complete_axis(ax)
if ~ishghandle(ax)
    return;
end
try
    ax.XRuler.TickLabelGapOffset = 4;
catch
end
try
    ax.YRuler.TickLabelGapOffset = 3;
catch
end
end

function set_compact_tick_density(ax, nx, ny)
if ~ishghandle(ax)
    return;
end
if nargin < 2 || isempty(nx)
    nx = 4;
end
if nargin < 3 || isempty(ny)
    ny = 4;
end

if strcmpi(ax.XScale, 'linear')
    xl = xlim(ax);
    if all(isfinite(xl)) && xl(2) > xl(1)
        xt = linspace(xl(1), xl(2), max(2, nx));
        if numel(unique(round(xt, 8))) >= 2
            xticks(ax, xt);
        end
    end
end

if strcmpi(ax.YScale, 'linear')
    yl = ylim(ax);
    if all(isfinite(yl)) && yl(2) > yl(1)
        yt = linspace(yl(1), yl(2), max(2, ny));
        if numel(unique(round(yt, 8))) >= 2
            yticks(ax, yt);
        end
    end
end
end

function txt = compact_legend_label(txt_in)
txt = char(txt_in);
txt = strrep(txt, 'BlueBrain', 'BB');
txt = strrep(txt, 'Canakci', 'Can');
txt = strrep(txt, 'Bezaire', 'Bez');
txt = strrep(txt, ' at Vg=', ' @');
txt = strrep(txt, ' Vg=', ' ');
txt = strrep(txt, ' mV', 'mV');
txt = strtrim(regexprep(txt, '\s+', ' '));
end

function exportgraphics_highres(fig_handle, out_path, preferred_dpi)
candidate_dpi = unique([preferred_dpi 1200 1000 900 800 700 600 500 400 300 240], 'stable');
last_err = [];
for k = 1:numel(candidate_dpi)
    dpi = candidate_dpi(k);
    try
        exportgraphics(fig_handle, out_path, 'Resolution', dpi, 'BackgroundColor', 'white');
        return;
    catch ME
        last_err = ME;
    end
end

% Final fallback via print if exportgraphics keeps failing in this environment.
try
    print(fig_handle, out_path, '-dpng', '-r300');
catch
    if ~isempty(last_err)
        rethrow(last_err);
    else
        error('Failed to export %s', out_path);
    end
end
end

function [channel_name, model_names] = figure18_panel_metadata(f, fallback_name)
channel_name = strtrim(fallback_name);
channel_name = regexprep(channel_name, '^Channel Comparison -\s*', '');
channel_name = regexprep(channel_name, '\s*\[[^\]]*\]\s*$', '');
model_names = {};

suptitle_h = findall(f, 'Type', 'Text', 'Tag', 'suptitle');
if ~isempty(suptitle_h)
    s = get(suptitle_h(1), 'String');
    if iscell(s)
        s = strjoin(s, ' ');
    end
    tok = regexp(char(s), 'Channel:\s*(.*?)\s*\|\s*Included:\s*(.*)$', 'tokens', 'once');
    if ~isempty(tok)
        channel_name = strtrim(tok{1});
        model_names = split_model_list(tok{2});
    end
end

if isempty(model_names)
    model_names = figure18_models_from_lines(f);
end

if isempty(model_names)
    model_names = {'n/a'};
end
end

function model_names = split_model_list(model_txt)
parts = strsplit(model_txt, ',');
parts = strtrim(parts);
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    model_names = {};
else
    model_names = unique(parts, 'stable');
end
end

function model_names = figure18_models_from_lines(f)
known_models = {'Canakci', 'BlueBrain', 'Bezaire'};
model_names = {};
lines = findall(f, 'Type', 'Line');
for i = 1:numel(lines)
    nm = get(lines(i), 'DisplayName');
    if isempty(nm) || ~ischar(nm)
        continue;
    end
    for m = 1:numel(known_models)
        if startsWith(strtrim(nm), known_models{m}, 'IgnoreCase', true)
            model_names{end+1} = known_models{m}; %#ok<AGROW>
        end
    end
end
if ~isempty(model_names)
    model_names = unique(model_names, 'stable');
end
end

function letter = figure18_panel_letter(idx)
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
letter = '';
n = idx;
while n > 0
    r = mod(n - 1, 26) + 1;
    letter = [alphabet(r) letter]; %#ok<AGROW>
    n = floor((n - 1) / 26);
end
end

function sec = normalize_section_name(sec_in)
sec = lower(strtrim(char(sec_in)));
if any(strcmp(sec, {'basal','basal_dend','basal-dend'}))
    sec = 'dend';
end
end

function out = sanitize_token(txt)
out = lower(strtrim(char(txt)));
out = regexprep(out, '[^a-z0-9]+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'na';
end
end

function c = model_color(name, col_can, col_bb, col_bz)
switch lower(name)
    case 'canakci'
        c = col_can;
    case 'bluebrain'
        c = col_bb;
    case 'bezaire'
        c = col_bz;
    otherwise
        c = [0 0 0];
end
end

function lbl = compact_model_label(name)
switch lower(strtrim(char(name)))
    case 'canakci'
        lbl = 'Can';
    case 'bluebrain'
        lbl = 'BB';
    case 'bezaire'
        lbl = 'Bez';
    otherwise
        lbl = strtrim(char(name));
end
end

function bb = bluebrain_extra_rates(V, celsius, cai)
% BlueBrain channels from this repo mod parameterizations.

% PC extras
[n,l,tn,tl] = bb_ka_rates(V, celsius, 11, -56, 0.05, 0.05, -1.5, 3, 0.55, 2, 0.1, -1, -40, 5, 1);
bb.pc.kap = struct('n',n,'l',l,'tau_n',tn,'tau_l',tl,'p_open',n.*l);

[m,tm] = bb_kmb_rates(V, celsius, -40, -10, -42, 0.003, 7, 0.4, 60, 0);
bb.pc.kmb = struct('m',m,'tau_m',tm,'p_open',m);

[m,tm] = bb_cal_rates(V, celsius, 5, 0.2, 0.1, 2, 4, 0.1);
h2 = cai ./ (cai + 50e-6);
bb.pc.cal = struct('m',m,'tau_m',tm,'p_open',m.^2 .* h2);

[m,tm,h,th] = bb_can_rates(V, celsius, 5, 0.2, 3, 0.03, 2, -14, 0.1);
h2 = cai ./ (cai + 50e-6);
bb.pc.can = struct('m',m,'h',h,'tau_m',tm,'tau_h',th,'p_open',m.^2 .* h .* h2);

[m,tm,h,th] = bb_cat_rates(V, celsius, 5, 0.2, 10, 0.04, 2, -28, 0.1, 0.015, 3.5, -75, 0.6);
bb.pc.cat = struct('m',m,'h',h,'tau_m',tm,'tau_h',th,'p_open',m.^2 .* h);

[m,tm] = bb_kca_rates(celsius, cai, 0.03, 0.00035, 0.5);
bb.pc.kca = struct('m',m.*ones(size(V)),'tau_m',tm.*ones(size(V)),'p_open',(m.^3).*ones(size(V)));

[o,to] = bb_cagk_rates(V, celsius, cai, 0.84, 1, 0.48e-3, 0.13e-6, 0.28, 0.48);
bb.pc.cagk = struct('o',o,'tau_o',to,'p_open',o);

% PVBC extras
bb.pv.kap = bb.pc.kap;
[n,l,tn,tl] = bb_ka_rates(V, celsius, -1, -56, 0.1, 0.05, -1.8, 3, 0.39, 2, 0.2, -1, -40, 5, 1);
bb.pv.kad = struct('n',n,'l',l,'tau_n',tn,'tau_l',tl,'p_open',n.*l);
bb.pv.kmb = bb.pc.kmb;
bb.pv.cal = bb.pc.cal;
bb.pv.can = bb.pc.can;
bb.pv.cat = bb.pc.cat;
bb.pv.kca = bb.pc.kca;
bb.pv.cagk = bb.pc.cagk;

[l,tl] = bb_hd_rates(V, celsius, -81, -8, -75, 0.011, 2.2, 0.4, 1, 4.5);
bb.pv.hd = struct('l',l,'tau_l',tl,'p_open',l);

[c,d,tc,td] = cavn_rates(V, celsius); % use Bezaire CavN equations for family-matched comparison
bb.pv.cavn = struct('c',c,'d',d,'tau_c',tc,'tau_d',td,'p_open',c.^2 .* d);

[q,tq] = kcas_rates(cai, celsius);
bb.pv.kcas = struct('q',q.*ones(size(V)),'tau_q',tq.*ones(size(V)),'p_open',(q.^2).*ones(size(V)));

[o,to] = kvcab_rates(V, cai, celsius);
bb.pv.kvcab = struct('o',o,'tau_o',to,'p_open',o);
end

function [ninf, linf, taun, taul] = bb_ka_rates(V, celsius, vhalfn, vhalfl, a0n, a0l, zetan, zetal, gmn, lmin, nmin, pw, tq, qq, qtl)
qt = 5.^((celsius-24)./10);
zeta = zetan + pw ./ (1 + exp((V - tq)./qq));
alpn = exp(1e-3 .* zeta .* (V - vhalfn) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
betn = exp(1e-3 .* zeta .* gmn .* (V - vhalfn) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
ninf = 1 ./ (1 + alpn);
taun = betn ./ (qt .* a0n .* (1 + alpn));
taun = max(taun, nmin);
alpl = exp(1e-3 .* zetal .* (V - vhalfl) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
linf = 1 ./ (1 + alpl);
taul = 0.26 .* (V + 50) ./ qtl;
taul = max(taul, lmin ./ qtl);
end

function [m_inf, tau_m] = bb_kmb_rates(V, celsius, vhalfl, kl, vhalft, a0t, zetat, gmt, b0, sh)
m_inf = 1 ./ (1 + exp((V - vhalfl - sh)./kl));
a = exp(0.0378 .* zetat .* (V - vhalft - sh));
b = exp(0.0378 .* zetat .* gmt .* (V - vhalft - sh));
tau_m = b0 + b ./ (a0t .* (1 + a));
celsius = celsius; %#ok<NASGU>
end

function [l_inf, tau_l] = bb_hd_rates(V, celsius, vhalfl, kl, vhalft, a0t, zetat, gmt, qtl, q10)
qt = q10.^((celsius-33)./10);
a = exp(0.0378 .* zetat .* (V - vhalft));
b = exp(0.0378 .* zetat .* gmt .* (V - vhalft));
l_inf = 1 ./ (1 + exp(-(V - vhalfl)./kl));
tau_l = b ./ (qtl .* qt .* a0t .* (1 + a));
end

function [m_inf, tau_m] = bb_cal_rates(V, celsius, q10, mmin, a0m, zetam, vhalfm, gmm)
qt = q10.^((celsius-25)./10);
a = 15.69 .* vtrap(-V + 81.5, 10);
b = 0.29 .* exp(-V./10.86);
m_inf = a ./ (a + b);
alpmt = exp(0.0378 .* zetam .* (V - vhalfm));
betmt = exp(0.0378 .* zetam .* gmm .* (V - vhalfm));
tau_m = betmt ./ (qt .* a0m .* (1 + alpmt));
tau_m = max(tau_m, mmin./qt);
end

function [m_inf, tau_m, h_inf, tau_h] = bb_can_rates(V, celsius, q10, mmin, hmin, a0m, zetam, vhalfm, gmm)
qt = q10.^((celsius-25)./10);
a = 0.1967 .* vtrap(-V + 19.88, 10);
b = 0.046 .* exp(-V./20.73);
m_inf = a ./ (a + b);
alpmt = exp(0.0378 .* zetam .* (V - vhalfm));
betmt = exp(0.0378 .* zetam .* gmm .* (V - vhalfm));
tau_m = betmt ./ (qt .* a0m .* (1 + alpmt));
tau_m = max(tau_m, mmin./qt);
ah = 1.6e-4 .* exp(-V./48.4);
bh = 1 ./ (exp((-V + 39.0)./10.0) + 1);
h_inf = ah ./ (ah + bh);
tau_h = max(80 .* ones(size(V)), hmin);
end

function [m_inf, tau_m, h_inf, tau_h] = bb_cat_rates(V, celsius, q10, mmin, hmin, a0m, zetam, vhalfm, gmm, a0h, zetah, vhalfh, gmh)
qt = q10.^((celsius-25)./10);
a = 0.2 .* vtrap(-V + 19.26, 10);
b = 0.009 .* exp(-V ./ 22.03);
m_inf = a ./ (a + b);
alpmt = exp(0.0378 .* zetam .* (V - vhalfm));
betmt = exp(0.0378 .* zetam .* gmm .* (V - vhalfm));
tau_m = betmt ./ (qt .* a0m .* (1 + alpmt));
tau_m = max(tau_m, mmin);
a_h = 1.0e-6 .* exp(-V ./ 16.26);
b_h = 1 ./ (exp((-V + 29.79)./10.0) + 1);
h_inf = a_h ./ (a_h + b_h);
alph_h = exp(0.0378 .* zetah .* (V - vhalfh));
beth_h = exp(0.0378 .* zetah .* gmh .* (V - vhalfh));
tau_h = beth_h ./ (a0h .* (1 + alph_h));
tau_h = max(tau_h, hmin);
end

function [m_inf, tau_m] = bb_kca_rates(celsius, cai, beta, cac, taumin)
tadj = 3.^((celsius-22)./10);
car = (cai./cac).^4;
m_inf = car ./ (1 + car);
tau_m = 1 ./ beta ./ (1 + car) ./ tadj;
tau_m = max(tau_m, taumin);
end

function [o_inf, tau_o] = bb_cagk_rates(V, celsius, cai, d1, d2, k1, k2, abar, bbar)
exp1a = k1 .* exp(-2 .* d1 .* 96485 .* V ./ (8.313424 .* (273.15 + celsius)));
exp1b = k2 .* exp(-2 .* d2 .* 96485 .* V ./ (8.313424 .* (273.15 + celsius)));
a = cai .* abar ./ (cai + exp1a);
b = bbar ./ (1 + cai ./ exp1b);
tau_o = 1 ./ (a + b);
o_inf = a .* tau_o;
end

function y = trap0(v, th, a, q)
x = v - th;
y = a .* x ./ (1 - exp(-x./q));
idx = abs(x) < 1e-6;
y(idx) = a .* q;
end

function y = vtrap(x,y0)
z = x ./ y0;
y = x ./ (exp(z)-1);
idx = abs(z) < 1e-6;
y(idx) = y0 .* (1 - z(idx)./2);
end
