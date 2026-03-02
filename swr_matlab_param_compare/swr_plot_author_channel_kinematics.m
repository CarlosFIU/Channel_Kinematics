function swr_plot_author_channel_kinematics(V)
% Plot shared channel kinetics across authors/models with available gate equations.
%
% Models plotted:
%   - Blue Brain PC soma
%   - Blue Brain PVBC soma
%   - Canakci 2017 CA1 pyramidal
%   - Canakci 2017 PV basket
%
% Notes:
%   - Some papers in the comparison set use LIF neurons or do not publish full
%     channel-rate equations in the manuscript text. Those are not plotted here.

if nargin < 1 || isempty(V)
    V = linspace(-100, 50, 600);
end

model_defs = struct( ...
    'label', { ...
        'BlueBrain PC (soma)', ...
        'BlueBrain PVBC (soma)', ...
        'Canakci 2017 Pyr', ...
        'Canakci 2017 PV' ...
    }, ...
    'builder', { ...
        @() hh_channel_kinematics_pc_params(), ...
        @() hh_channel_kinematics_pvbc_params(), ...
        @() hh_canakci_pyr_params(), ...
        @() hh_canakci_pv_params(34) ...
    });

schema = struct( ...
    'label', '', ...
    'm_inf', [], ...
    'h_inf', [], ...
    'n_inf', [], ...
    'tau_m', [], ...
    'tau_h', [], ...
    'tau_n', [], ...
    'p_open_na', [], ...
    'p_open_k', []);

data = repmat(schema, 0, 1);
for k = 1:numel(model_defs)
    hh = model_defs(k).builder();
    g = hh_compute_gating(V, hh);

    d = schema;
    d.label = model_defs(k).label;
    d.m_inf = field_or_nan(g, 'm_inf', V);
    d.h_inf = field_or_nan(g, 'h_inf', V);
    d.n_inf = field_or_nan(g, 'n_inf', V);
    d.tau_m = field_or_nan(g, 'tau_m_ms', V);
    d.tau_h = field_or_nan(g, 'tau_h_ms', V);
    d.tau_n = field_or_nan(g, 'tau_n_ms', V);

    if isfield(g, 'p_open_Na')
        d.p_open_na = g.p_open_Na;
    else
        d.p_open_na = d.m_inf.^3 .* d.h_inf;
    end

    if isfield(g, 'p_open_Kdr')
        d.p_open_k = g.p_open_Kdr;
    elseif isfield(g, 'p_open_Kdrb')
        d.p_open_k = g.p_open_Kdrb;
    else
        d.p_open_k = d.n_inf;
    end

    data(end+1,1) = d; %#ok<AGROW>
end

cols = lines(numel(data));

figure('Name', 'Cross-paper Na/K gate comparison', 'Color', 'w');
tiledlayout(2,3, 'Padding','compact', 'TileSpacing','compact');

nexttile; hold on;
for k = 1:numel(data)
    plot(V, data(k).m_inf, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('m_inf'); title('Na activation (open)');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    plot(V, data(k).h_inf, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('h_inf'); title('Na inactivation (open)');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    plot(V, data(k).n_inf, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('n_inf'); title('K activation (open)');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    if all(isnan(data(k).tau_m))
        continue;
    end
    plot(V, data(k).tau_m, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('tau_m (ms)'); title('Na tau_m');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    if all(isnan(data(k).tau_h))
        continue;
    end
    plot(V, data(k).tau_h, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('tau_h (ms)'); title('Na tau_h');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    if all(isnan(data(k).tau_n))
        continue;
    end
    plot(V, data(k).tau_n, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('tau_n (ms)'); title('K tau_n');
grid on; legend('Location','best');

set(gcf, 'Position', [120 80 1400 760]);

figure('Name', 'Cross-paper channel opening probability', 'Color', 'w');
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

nexttile; hold on;
for k = 1:numel(data)
    plot(V, data(k).p_open_na, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('P_open'); title('Na opening probability');
grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(data)
    plot(V, data(k).p_open_k, 'LineWidth', 1.5, 'Color', cols(k,:), ...
        'DisplayName', data(k).label);
end
xlabel('V (mV)'); ylabel('P_open'); title('K opening probability');
grid on; legend('Location','best');

set(gcf, 'Position', [140 180 1200 360]);

save_fig_exports();

fprintf(['\nCross-paper channel overlays include models with explicit gate-rate equations.\n' ...
         'LIF-only papers and papers without full published gate equations are listed in the tables but skipped in gating plots.\n']);

end

function x = field_or_nan(s, f, V)
if isfield(s, f)
    x = s.(f);
else
    x = nan(size(V));
end
end

function save_fig_exports()
out_dir = fullfile(pwd, 'exports_channel_kinematics');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

figs = findall(groot, 'Type', 'figure');
for k = 1:numel(figs)
    f = figs(k);
    fig_name = get(f, 'Name');
    safe_name = regexprep(fig_name, '[^A-Za-z0-9_\- ]', '_');
    if isempty(strtrim(safe_name))
        safe_name = sprintf('paper_compare_%02d', k);
    end
    out_png = fullfile(out_dir, ['paper_compare_' safe_name '.png']);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(f, out_png, 'Resolution', 180);
    else
        print(f, out_png, '-dpng', '-r180');
    end
end
end
