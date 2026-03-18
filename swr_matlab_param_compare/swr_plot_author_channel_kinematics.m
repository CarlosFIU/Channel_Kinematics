function swr_plot_author_channel_kinematics(V_shared)
% Plot explicit channel kinetics across source models with recovered gate equations.
%
% Voltage axis convention:
%   V_shared is interpreted on a shared resting-potential axis so all models
%   are aligned to the same reference rest before evaluating their native
%   gate equations.

if nargin < 1 || isempty(V_shared)
    V_shared = linspace(-100, 50, 600);
end

cai_uM = linspace(0, 500, 400);
shared_vrest_mv = shared_author_resting_potential_mv();

model_defs = {
    'BlueBrain PC (soma)',     '-',  [0.00 0.45 0.74], 'bb_pc';
    'BlueBrain PVBC (soma)',   '--', [0.85 0.33 0.10], 'bb_pv';
    'Canakci 2017 Pyr',        ':',  [0.49 0.18 0.56], 'can_pc';
    'Canakci 2017 PV',         '-.', [0.93 0.69 0.13], 'can_pv';
    'Traub 1991 CA3 Pyr',      '-',  [0.20 0.60 0.20], 'traub_ca3';
    'Traub/Miles 1995 INT',    '--', [0.64 0.08 0.18], 'traub_int'
};

schema = init_author_plot_schema(V_shared, cai_uM);
data = repmat(schema, 0, 1);
for k = 1:size(model_defs, 1)
    label = model_defs{k, 1};
    line_style = model_defs{k, 2};
    color = model_defs{k, 3};
    kind = model_defs{k, 4};
    d = build_author_model_data(kind, label, line_style, color, V_shared, cai_uM, shared_vrest_mv);
    data(end+1,1) = d; %#ok<AGROW>
end

fig_handles = gobjects(0,1);

f = figure('Name', 'Cross-paper Na_K opening probability and tau', 'Color', 'w');
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot_author_field(data, V_shared, 'p_open_na');
xlabel('V aligned to shared rest (mV)'); ylabel('P_{open}'); title('Na opening probability');
grid on; legend('Location', 'best');

nexttile; hold on;
plot_author_tau_fields(data, V_shared, {'tau_m', 'tau_h'}, {'m', 'h'});
xlabel('V aligned to shared rest (mV)'); ylabel('\tau (ms)'); title('Na tau constants');
grid on; legend('Location', 'best');

nexttile; hold on;
plot_author_field(data, V_shared, 'p_open_k');
xlabel('V aligned to shared rest (mV)'); ylabel('P_{open}'); title('K opening probability');
grid on; legend('Location', 'best');

nexttile; hold on;
plot_author_tau_fields(data, V_shared, {'tau_n'}, {'n'});
xlabel('V aligned to shared rest (mV)'); ylabel('\tau (ms)'); title('K tau constants');
grid on; legend('Location', 'best');

set(f, 'Position', [120 80 1360 760]);
fig_handles(end+1) = f; %#ok<AGROW>

traub_idx = find(contains({data.model_group}, 'traub'));
if ~isempty(traub_idx)
    f = figure('Name', 'Traub single-cell additional channel opening and tau', 'Color', 'w');
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    plot_author_field(data(traub_idx), V_shared, 'p_open_ka');
    xlabel('V aligned to shared rest (mV)'); ylabel('P_{open}'); title('K(A) opening probability');
    grid on; legend('Location', 'best');

    nexttile; hold on;
    plot_author_tau_fields(data(traub_idx), V_shared, {'tau_a', 'tau_b'}, {'a', 'b'});
    xlabel('V aligned to shared rest (mV)'); ylabel('\tau (ms)'); title('K(A) tau constants');
    grid on; legend('Location', 'best');

    nexttile; hold on;
    plot_author_field(data(traub_idx), V_shared, 'p_open_ca');
    xlabel('V aligned to shared rest (mV)'); ylabel('P_{open}'); title('Ca opening probability');
    grid on; legend('Location', 'best');

    nexttile; hold on;
    plot_author_tau_fields(data(traub_idx), V_shared, {'tau_s', 'tau_r'}, {'s', 'r'});
    xlabel('V aligned to shared rest (mV)'); ylabel('\tau (ms)'); title('Ca tau constants');
    grid on; legend('Location', 'best');

    set(f, 'Position', [100 90 1360 760]);
    fig_handles(end+1) = f; %#ok<AGROW>

    f = figure('Name', 'Traub calcium-dependent channel probability', 'Color', 'w');
    tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    plot_author_cai_field(data(traub_idx), cai_uM, 'p_open_kc_cai');
    xlabel('cai (\muM)'); ylabel('P_{open}'); title('K(C) probability @ shared V=-40 mV');
    grid on; legend('Location', 'best');

    nexttile; hold on;
    plot_author_cai_field(data(traub_idx), cai_uM, 'p_open_kahp_cai');
    xlabel('cai (\muM)'); ylabel('P_{open}'); title('K(AHP) probability');
    grid on; legend('Location', 'best');

    set(f, 'Position', [120 180 1200 360]);
    fig_handles(end+1) = f; %#ok<AGROW>
end

save_fig_exports(fig_handles);

fprintf(['\nCross-paper channel overlays now include the explicit Traub CA3 pyramidal model ' ...
         'and the Traub/Miles 1995 interneuron kinetics.\n' ...
         'For the interneuron, the gate equations are inferred from the paper statement ' ...
         'that they share the same form as the Traub CA3 pyramidal model, with different geometry and conductance densities.\n']);
end

function d = init_author_plot_schema(V_rel, cai_uM)
d = struct( ...
    'label', '', ...
    'line_style', '-', ...
    'color', [0 0 0], ...
    'model_group', '', ...
    'm_inf', nan(size(V_rel)), ...
    'h_inf', nan(size(V_rel)), ...
    'n_inf', nan(size(V_rel)), ...
    'tau_m', nan(size(V_rel)), ...
    'tau_h', nan(size(V_rel)), ...
    'tau_n', nan(size(V_rel)), ...
    'p_open_na', nan(size(V_rel)), ...
    'p_open_k', nan(size(V_rel)), ...
    'a_inf', nan(size(V_rel)), ...
    'b_inf', nan(size(V_rel)), ...
    'tau_a', nan(size(V_rel)), ...
    'tau_b', nan(size(V_rel)), ...
    's_inf', nan(size(V_rel)), ...
    'r_inf', nan(size(V_rel)), ...
    'tau_s', nan(size(V_rel)), ...
    'tau_r', nan(size(V_rel)), ...
    'c_inf', nan(size(V_rel)), ...
    'tau_c', nan(size(V_rel)), ...
    'xi_cai', nan(size(cai_uM)), ...
    'q_inf_cai', nan(size(cai_uM)), ...
    'p_open_ka', nan(size(V_rel)), ...
    'p_open_ca', nan(size(V_rel)), ...
    'p_open_kc_cai', nan(size(cai_uM)), ...
    'p_open_kahp_cai', nan(size(cai_uM)));
end

function d = build_author_model_data(kind, label, line_style, color, V_shared, cai_uM, shared_vrest_mv)
d = init_author_plot_schema(V_shared, cai_uM);
d.label = label;
d.line_style = line_style;
d.color = color;
d.model_group = kind;

switch lower(kind)
    case 'bb_pc'
        d = fill_existing_hh_model(d, author_aligned_voltage(V_shared, -70, shared_vrest_mv), @() hh_channel_kinematics_pc_params(), 'p_open_Kdr');
    case 'bb_pv'
        d = fill_existing_hh_model(d, author_aligned_voltage(V_shared, -70, shared_vrest_mv), @() hh_channel_kinematics_pvbc_params(), 'p_open_Kdrb');
    case 'can_pc'
        d = fill_existing_hh_model(d, author_aligned_voltage(V_shared, -70, shared_vrest_mv), @() hh_canakci_pyr_params(), '');
    case 'can_pv'
        d = fill_existing_hh_model(d, author_aligned_voltage(V_shared, -65, shared_vrest_mv), @() hh_canakci_pv_params(34), '');
    case {'traub_ca3','traub_int'}
        d = fill_traub_model(d, author_aligned_voltage(V_shared, 0, shared_vrest_mv), cai_uM, ...
            author_aligned_voltage(-40, 0, shared_vrest_mv));
    otherwise
        error('Unsupported author-plot model kind: %s', kind);
end
end

function V_eval = author_aligned_voltage(V_shared, model_rest_mv, shared_vrest_mv)
V_eval = V_shared + (model_rest_mv - shared_vrest_mv);
end

function vrest_mv = shared_author_resting_potential_mv()
vrest_mv = -70;
end

function d = fill_existing_hh_model(d, V_native, hh_builder, k_field)
hh = hh_builder();
g = hh_compute_gating(V_native, hh);

d.m_inf = field_or_nan(g, 'm_inf', V_native);
d.h_inf = field_or_nan(g, 'h_inf', V_native);
d.n_inf = field_or_nan(g, 'n_inf', V_native);
d.tau_m = field_or_nan(g, 'tau_m_ms', V_native);
d.tau_h = field_or_nan(g, 'tau_h_ms', V_native);
d.tau_n = field_or_nan(g, 'tau_n_ms', V_native);

if isfield(g, 'p_open_Na')
    d.p_open_na = g.p_open_Na;
elseif all(~isnan(d.m_inf)) && all(~isnan(d.h_inf))
    d.p_open_na = d.m_inf.^3 .* d.h_inf;
end

if ~isempty(k_field) && isfield(g, k_field)
    d.p_open_k = g.(k_field);
elseif all(~isnan(d.n_inf))
    d.p_open_k = d.n_inf;
end
end

function d = fill_traub_model(d, V_rel, cai_uM, Vfix_rel)
[m_inf, h_inf, tau_m, tau_h] = traub_na_rates(V_rel);
[n_inf, tau_n] = traub_kdr_rates(V_rel);
[a_inf, b_inf, tau_a, tau_b] = traub_ka_rates(V_rel);
[s_inf, r_inf, tau_s, tau_r] = traub_ca_rates(V_rel);
[c_inf, tau_c] = traub_kc_rates(V_rel);
[q_inf, ~] = traub_kahp_rates(cai_uM);
[c_fix, ~] = traub_kc_rates(Vfix_rel);

d.m_inf = m_inf;
d.h_inf = h_inf;
d.n_inf = n_inf;
d.tau_m = tau_m;
d.tau_h = tau_h;
d.tau_n = tau_n;
d.p_open_na = m_inf.^2 .* h_inf;
d.p_open_k = n_inf;

d.a_inf = a_inf;
d.b_inf = b_inf;
d.tau_a = tau_a;
d.tau_b = tau_b;
d.s_inf = s_inf;
d.r_inf = r_inf;
d.tau_s = tau_s;
d.tau_r = tau_r;
d.c_inf = c_inf;
d.tau_c = tau_c;
d.p_open_ka = a_inf .* b_inf;
d.p_open_ca = (s_inf .^ 2) .* r_inf;
d.xi_cai = min(cai_uM ./ 250, 1);
d.q_inf_cai = q_inf;
d.p_open_kc_cai = c_fix .* d.xi_cai;
d.p_open_kahp_cai = q_inf;
end

function plot_author_field(data, x, field_name)
for k = 1:numel(data)
    y = data(k).(field_name);
    if all(isnan(y))
        continue;
    end
    plot(x, y, 'LineWidth', 1.5, 'LineStyle', data(k).line_style, ...
        'Color', data(k).color, 'DisplayName', data(k).label);
end
end

function plot_author_cai_field(data, cai_uM, field_name)
for k = 1:numel(data)
    y = data(k).(field_name);
    if all(isnan(y))
        continue;
    end
    plot(cai_uM, y, 'LineWidth', 1.5, 'LineStyle', data(k).line_style, ...
        'Color', data(k).color, 'DisplayName', data(k).label);
end
end

function plot_author_tau_fields(data, x, field_names, gate_labels)
gate_styles = {'-', '--', ':', '-.'};
for k = 1:numel(data)
    for j = 1:numel(field_names)
        y = data(k).(field_names{j});
        if all(isnan(y))
            continue;
        end
        plot(x, y, 'LineWidth', 1.5, ...
            'LineStyle', gate_styles{mod(j-1, numel(gate_styles))+1}, ...
            'Color', data(k).color, ...
            'DisplayName', sprintf('%s %s', data(k).label, gate_labels{j}));
    end
end
end

function x = field_or_nan(s, f, V)
if isfield(s, f)
    x = s.(f);
else
    x = nan(size(V));
end
end

function save_fig_exports(figs)
out_dir = fullfile(pwd, 'exports_channel_kinematics');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

for k = 1:numel(figs)
    f = figs(k);
    if ~ishandle(f) || ~strcmp(get(f, 'Type'), 'figure')
        continue;
    end
    fig_name = get(f, 'Name');
    safe_name = regexprep(fig_name, '[^A-Za-z0-9_\- ]', '_');
    if isempty(strtrim(safe_name))
        safe_name = sprintf('paper_compare_%02d', k);
    end
    out_png = fullfile(out_dir, ['paper_compare_' safe_name '.png']);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(f, out_png, 'Resolution', 220);
    else
        print(f, out_png, '-dpng', '-r220');
    end
end
end

function [m_inf, h_inf, tau_m, tau_h] = traub_na_rates(V)
alpha = 0.32 .* trap_exp_ratio(13.1 - V, 4);
beta = 0.28 .* trap_exp_ratio(V - 40.1, 5);
tau_m = 1 ./ (alpha + beta);
m_inf = alpha .* tau_m;

alpha = 0.128 .* exp((17 - V) ./ 18);
beta = 4 ./ (1 + exp((40 - V) ./ 5));
tau_h = 1 ./ (alpha + beta);
h_inf = alpha .* tau_h;
end

function [n_inf, tau_n] = traub_kdr_rates(V)
alpha = 0.016 .* trap_exp_ratio(35.1 - V, 5);
beta = 0.25 .* exp((20 - V) ./ 40);
tau_n = 1 ./ (alpha + beta);
n_inf = alpha .* tau_n;
end

function [a_inf, b_inf, tau_a, tau_b] = traub_ka_rates(V)
alpha = 0.02 .* trap_exp_ratio(13.1 - V, 10);
beta = 0.0175 .* trap_exp_ratio(V - 40.1, 10);
tau_a = 1 ./ (alpha + beta);
a_inf = alpha .* tau_a;

alpha = 0.0016 .* exp((-13 - V) ./ 18);
beta = 0.05 ./ (1 + exp((10.1 - V) ./ 5));
tau_b = 1 ./ (alpha + beta);
b_inf = alpha .* tau_b;
end

function [s_inf, r_inf, tau_s, tau_r] = traub_ca_rates(V)
alpha = 1.6 ./ (1 + exp(-0.072 .* (V - 65)));
beta = 0.02 .* trap_exp_ratio(V - 51.1, 5);
tau_s = 1 ./ (alpha + beta);
s_inf = alpha .* tau_s;

alpha = 0.005 .* ones(size(V));
idx = V > 0;
alpha(idx) = exp(-V(idx) ./ 20) ./ 200;
beta = zeros(size(V));
beta(idx) = 0.005 - exp(-V(idx) ./ 20) ./ 200;
tau_r = 1 ./ (alpha + beta);
r_inf = alpha .* tau_r;
end

function [c_inf, tau_c] = traub_kc_rates(V)
alpha = zeros(size(V));
beta = zeros(size(V));
idx = V <= 50;
alpha(idx) = exp((V(idx) - 10) ./ 11 - (V(idx) - 6.5) ./ 27) ./ 18.975;
alpha(~idx) = 2 .* exp(-(V(~idx) - 6.5) ./ 27);
beta(idx) = 2 .* exp(-(V(idx) - 6.5) ./ 27) - alpha(idx);
tau_c = 1 ./ (alpha + beta);
c_inf = alpha .* tau_c;
end

function [q_inf, tau_q] = traub_kahp_rates(cai_uM)
alpha = 0.2e-4 .* cai_uM;
alpha(alpha > 0.01) = 0.01;
beta = 0.001;
tau_q = 1 ./ (alpha + beta);
q_inf = alpha .* tau_q;
end

function y = trap_exp_ratio(x, q)
z = x ./ q;
y = x ./ (exp(z) - 1);
small = abs(z) < 1e-6;
y(small) = q .* (1 - z(small) ./ 2);
end
