function swr_plot_pc_pvbc_targeted_kinematics(V)
% Targeted PC/PVBC channel kinematics across Canakci, BlueBrain, Bezaire.
% - Voltage-dependent gate/open probability comparison
% - Ca_i sweeps only for calcium-dependent channels

if nargin < 1 || isempty(V)
    V = linspace(-100, 50, 600);
end

pc_sections = {'soma'};
pv_sections = {'soma'};

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
shared_defs{end+1} = make_channel_def('PC Na', entries); %#ok<AGROW>

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
shared_defs{end+1} = make_channel_def('PC KDR', entries); %#ok<AGROW>

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
shared_defs{end+1} = make_channel_def('PVBC Na', entries); %#ok<AGROW>

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
shared_defs{end+1} = make_channel_def('PVBC KDR', entries); %#ok<AGROW>

% PC A-type K (family-matched)
entries = struct([]);
bb_pc_kap_secs = {'soma'};
for s = 1:numel(bb_pc_kap_secs)
    sec = bb_pc_kap_secs{s};
    if isfield(bb_pc_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n','l'}, ...
            {bb_extra.pc.kap.n, bb_extra.pc.kap.l}, ...
            {bb_extra.pc.kap.tau_n, bb_extra.pc.kap.tau_l}, ...
            bb_extra.pc.kap.p_open);
    end
end
entries = add_channel_entry(entries, 'Bezaire', 'soma', ...
    {'n','l'}, {bz_pc_by_sec.soma.n_kvap, bz_pc_by_sec.soma.l_kvap}, ...
    {bz_pc_by_sec.soma.tau_n_kvap, bz_pc_by_sec.soma.tau_l_kvap}, ...
    bz_pc_by_sec.soma.p_kvap);
shared_defs{end+1} = make_channel_def('PC A-type K', entries); %#ok<AGROW>

% PVBC A-type K (family-matched)
entries = struct([]);
bb_pv_kap_secs = {'soma'};
for s = 1:numel(bb_pv_kap_secs)
    sec = bb_pv_kap_secs{s};
    if isfield(bb_pv_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'n','l'}, ...
            {bb_extra.pv.kap.n, bb_extra.pv.kap.l}, ...
            {bb_extra.pv.kap.tau_n, bb_extra.pv.kap.tau_l}, ...
            bb_extra.pv.kap.p_open);
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
shared_defs{end+1} = make_channel_def('PVBC A-type K', entries); %#ok<AGROW>

for k = 1:numel(shared_defs)
    fig_handles(end+1) = plot_channel_definition(V, shared_defs{k}, out_dir, 'shared', col_can, col_bb, col_bz); %#ok<AGROW>
end

% -------------------------------------------------------------------------
% 2) Individual model-specific channels after shared ones
% -------------------------------------------------------------------------
ind_defs = {};

% BlueBrain PC-only families
entries = struct([]);
for s = 1:numel(bb_pc_kap_secs)
    sec = bb_pc_kap_secs{s};
    if isfield(bb_pc_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m'}, {bb_extra.pc.kmb.m}, {bb_extra.pc.kmb.tau_m}, bb_extra.pc.kmb.p_open);
    end
end
ind_defs{end+1} = make_channel_def('BlueBrain PC M-type K', entries); %#ok<AGROW>

entries = struct([]);
bb_pc_cat_secs = {'soma'};
for s = 1:numel(bb_pc_cat_secs)
    sec = bb_pc_cat_secs{s};
    if isfield(bb_pc_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h'}, ...
            {bb_extra.pc.cat.m, bb_extra.pc.cat.h}, ...
            {bb_extra.pc.cat.tau_m, bb_extra.pc.cat.tau_h}, ...
            bb_extra.pc.cat.p_open);
    end
end
ind_defs{end+1} = make_channel_def('BlueBrain PC T-type Ca', entries); %#ok<AGROW>

% BlueBrain PVBC-only families
entries = struct([]);
if isfield(bb_pv_by_sec, 'soma') && section_has_gbar(bb_pv_by_sec.soma, 'gbar_Kdb')
    g = bb_pv_by_sec.soma;
    entries = add_channel_entry(entries, 'BlueBrain', 'soma', ...
        {'n'}, {g.n_kdb_inf}, {g.tau_n_kdb_ms}, g.p_open_Kdb);
end
ind_defs{end+1} = make_channel_def('BlueBrain PVBC KDB', entries); %#ok<AGROW>

entries = struct([]);
for s = 1:numel(bb_pv_kap_secs)
    sec = bb_pv_kap_secs{s};
    if isfield(bb_pv_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m'}, {bb_extra.pv.kmb.m}, {bb_extra.pv.kmb.tau_m}, bb_extra.pv.kmb.p_open);
    end
end
ind_defs{end+1} = make_channel_def('BlueBrain PVBC M-type K', entries); %#ok<AGROW>

entries = struct([]);
bb_pv_hd_secs = {'soma'};
for s = 1:numel(bb_pv_hd_secs)
    sec = bb_pv_hd_secs{s};
    if isfield(bb_pv_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'l'}, {bb_extra.pv.hd.l}, {bb_extra.pv.hd.tau_l}, bb_extra.pv.hd.p_open);
    end
end
ind_defs{end+1} = make_channel_def('BlueBrain PVBC HCN', entries); %#ok<AGROW>

entries = struct([]);
for s = 1:numel(bb_pv_hd_secs)
    sec = bb_pv_hd_secs{s};
    if isfield(bb_pv_by_sec, sec)
        entries = add_channel_entry(entries, 'BlueBrain', sec, ...
            {'m','h'}, ...
            {bb_extra.pv.cat.m, bb_extra.pv.cat.h}, ...
            {bb_extra.pv.cat.tau_m, bb_extra.pv.cat.tau_h}, ...
            bb_extra.pv.cat.p_open);
    end
end
ind_defs{end+1} = make_channel_def('BlueBrain PVBC T-type Ca', entries); %#ok<AGROW>

% Bezaire PC-only channels that are non-somatic are omitted in soma-only mode.

for k = 1:numel(ind_defs)
    fig_handles(end+1) = plot_channel_definition(V, ind_defs{k}, out_dir, 'individual', col_can, col_bb, col_bz); %#ok<AGROW>
end

% --------------------------
% Calcium-only sweeps for calcium-dependent channels
% --------------------------
cai = logspace(log10(1e-6), log10(2e-3), 350); % mM
Vfix = -40;

% Bezaire PVBC Ca-dependent channels.
[p_kcas_cai, p_kvcab_cai, p_cavl_h2_cai] = bezaire_ca_sweeps(cai, Vfix, 34);

% BlueBrain-style Ca dependent factors used in local mod set.
p_bb_kca = (cai ./ 0.00035).^4 ./ (1 + (cai ./ 0.00035).^4);     % m_inf
p_bb_kca = p_bb_kca.^3;                                           % P_open
p_bb_cah2 = cai ./ (cai + 50e-6);                                % cal/can factor

f_cai = figure('Name','Calcium-dependent channels (cai sweeps only)','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(cai*1e3, p_bb_kca, 'LineWidth', 1.6, 'DisplayName', 'BlueBrain kca P_{open}');
plot(cai*1e3, p_kcas_cai, 'LineWidth', 1.6, 'DisplayName', 'Bezaire ch\_KCaS P_{open}');
set(gca,'XScale','log');
xlabel('cai (uM)'); ylabel('P_{open}'); title('Ca-dependent K activation');
grid on; legend('Location','best');

nexttile; hold on;
plot(cai*1e3, p_bb_cah2, 'LineWidth', 1.6, 'DisplayName', 'BlueBrain cal/can Ca-factor');
plot(cai*1e3, p_cavl_h2_cai, 'LineWidth', 1.6, 'DisplayName', 'Bezaire ch\_CavL h2(cai)');
set(gca,'XScale','log');
xlabel('cai (uM)'); ylabel('Ca factor'); title('Ca-factor used by Ca channels');
grid on; legend('Location','best');

nexttile; hold on;
plot(cai*1e3, p_kvcab_cai, 'LineWidth', 1.6, 'DisplayName', sprintf('Bezaire ch\\_KvCaB, V=%g mV', Vfix));
set(gca,'XScale','log');
xlabel('cai (uM)'); ylabel('P_{open}'); title('BK-like Ca+V activation');
grid on; legend('Location','best');

set(gcf,'Position',[120 220 1250 360]);
exportgraphics(f_cai, fullfile(out_dir,'paper_compare_PC_PVBC_cai_sweeps.png'), ...
    'Resolution',180);
fig_handles(end+1) = f_cai; %#ok<AGROW>

export_target_figures_to_powerpoint(fig_handles, out_dir, ...
    'paper_compare_PC_PVBC_targeted_kinematics.pptx');
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
    catch
        % If a section is missing required parameters for this cell/model,
        % leave it out and continue plotting available sections.
    end
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

function f = plot_channel_definition(V, d, out_dir, prefix, col_can, col_bb, col_bz)
gate_styles = {'-','--',':','-.'};
popen_styles = {'-','--',':','-.'};
f = figure('Name', ['Channel Comparison - ' d.name], 'Color', 'w');
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
title([d.name ' gates (open)']); grid on; place_legend_safely(gca);

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
title([d.name ' gates (closing)']); grid on; place_legend_safely(gca);

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
title([d.name ' tau constants']); grid on; place_legend_safely(gca);

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
title([d.name ' opening probability']); grid on; place_legend_safely(gca);

included = cell(1, numel(d.entries));
for m = 1:numel(d.entries)
    included{m} = d.entries(m).model;
end
included = unique(included, 'stable');
sgtitle(sprintf('Channel: %s | Included: %s', d.name, strjoin(included, ', ')), ...
    'Interpreter', 'none');
set(f, 'Position', [130 80 1300 760]);

safe = regexprep([prefix '_' d.name], '[^A-Za-z0-9_\- ]', '_');
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
    delete(ppt_path);
end

ppt = Presentation(ppt_path);
open(ppt);

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
    png_path = fullfile(out_dir, sprintf('ppt_%02d_%s.png', k, safe_name));
    try
        exportgraphics(f, png_path, 'Resolution', 180);
    catch
        % Fallback for MATLAB versions/environments where exportgraphics
        % cannot consume this figure-handle representation.
        print(f, png_path, '-dpng', '-r180');
    end

    slide = add(ppt, 'Title and Content');
    replace(slide, 'Title', fig_name);
    replace(slide, 'Content', Picture(png_path));
    slide_count = slide_count + 1;
end

close(ppt);
fprintf('PowerPoint saved with %d slides: %s\n', slide_count, ppt_path);
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
