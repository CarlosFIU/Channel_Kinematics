function hh_channel_kinematics_compartment_compare()
% Compare channel opening/closing and tau across compartments for all channels
% present in the SONATA JSON files for PC and PVBC.
%
% PC sections: soma, dend, apic
% PVBC sections: soma, dend
%
% For each mechanism found in the JSON conductance parameters, this script plots:
%   1) opening variables (gate inf)
%   2) closing variables (1 - gate inf)
%   3) tau constants
%   4) conductance-weighted open probability gbar * P_open

close all;
V = linspace(-100, 50, 600);

pc_sections = {'soma','dend','apic'};
pv_sections = {'soma','dend'};

pc = load_model_channel_data(fullfile('..','PC_dynamics_params_sonata.json'), 'PC', pc_sections, V);
pv = load_model_channel_data(fullfile('..','PVBC_dynamics_params_sonata.json'), 'PVBC', pv_sections, V);

plot_model_channels(pc, V);
plot_model_channels(pv, V);
print_channel_summary(pc);
print_channel_summary(pv);

end

function model = load_model_channel_data(json_path, model_name, sections, V)
raw = jsondecode(fileread(json_path));
model.name = model_name;
model.celsius = raw.conditions(1).celsius;
model.sections = sections;
model.V = V;

all_mechs = {};
for s = 1:numel(sections)
    sec = sections{s};
    idx = find(arrayfun(@(g) strcmp(g.section, sec) && contains(g.name, 'bar'), raw.genome));
    model.sec(s).name = sec; %#ok<AGROW>
    model.sec(s).channels = struct(); %#ok<AGROW>
    for k = 1:numel(idx)
        g = raw.genome(idx(k));
        mech = g.mechanism;
        if isempty(mech)
            continue;
        end
        all_mechs{end+1} = mech; %#ok<AGROW>
        ch.gbar = g.value;
        ch.param_name = g.name;
        kin = mechanism_kinetics(mech, V, model.celsius);
        ch.available = kin.available;
        ch.display_name = kin.display_name;
        ch.gate_names = kin.gate_names;
        ch.open_inf = kin.open_inf;
        ch.close_inf = 1 - kin.open_inf;
        ch.tau_ms = kin.tau_ms;
        ch.p_open = kin.p_open;
        ch.gbar_p_open = g.value .* kin.p_open;
        model.sec(s).channels.(mech) = ch;
    end
end
model.mechanisms = unique(all_mechs, 'stable');
end

function plot_model_channels(model, V)
for m = 1:numel(model.mechanisms)
    mech = model.mechanisms{m};

    fig_name = sprintf('%s: %s', model.name, channel_title(mech));
    figure('Name', fig_name, 'Color', 'w');
    cols = lines(numel(model.sections));

    subplot(2,2,1); hold on;
    for s = 1:numel(model.sections)
        sec_name = model.sections{s};
        if ~has_channel(model, s, mech)
            continue;
        end
        ch = model.sec(s).channels.(mech);
        for g = 1:numel(ch.gate_names)
            style = gate_line_style(g);
            plot(V, ch.open_inf(g,:), style, 'Color', cols(s,:), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s %s', sec_name, ch.gate_names{g}));
        end
    end
    xlabel('V (mV)'); ylabel('Opening variable x_\infty');
    title(sprintf('%s opening', channel_title(mech)));
    grid on; legend('Location','bestoutside');

    subplot(2,2,2); hold on;
    for s = 1:numel(model.sections)
        if ~has_channel(model, s, mech)
            continue;
        end
        sec_name = model.sections{s};
        ch = model.sec(s).channels.(mech);
        for g = 1:numel(ch.gate_names)
            style = gate_line_style(g);
            plot(V, ch.close_inf(g,:), style, 'Color', cols(s,:), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s %s', sec_name, ch.gate_names{g}));
        end
    end
    xlabel('V (mV)'); ylabel('Closing variable (1-x_\infty)');
    title(sprintf('%s closing', channel_title(mech)));
    grid on; legend('Location','bestoutside');

    subplot(2,2,3); hold on;
    for s = 1:numel(model.sections)
        if ~has_channel(model, s, mech)
            continue;
        end
        sec_name = model.sections{s};
        ch = model.sec(s).channels.(mech);
        if isempty(ch.tau_ms)
            continue;
        end
        for g = 1:size(ch.tau_ms,1)
            style = gate_line_style(g);
            plot(V, ch.tau_ms(g,:), style, 'Color', cols(s,:), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s tau_%s', sec_name, ch.gate_names{g}));
        end
    end
    xlabel('V (mV)'); ylabel('Tau (ms)');
    title(sprintf('%s tau constants', channel_title(mech)));
    grid on; legend('Location','bestoutside');

    subplot(2,2,4); hold on;
    for s = 1:numel(model.sections)
        if ~has_channel(model, s, mech)
            continue;
        end
        sec_name = model.sections{s};
        ch = model.sec(s).channels.(mech);
        plot(V, ch.gbar_p_open, '-', 'Color', cols(s,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s %s', sec_name, ch.param_name));
    end
    xlabel('V (mV)'); ylabel('gbar * P_{open} (S/cm^2)');
    title(sprintf('%s effective opening', channel_title(mech)));
    grid on; legend('Location','bestoutside');

    set(gcf, 'Position', [80 60 1400 820]);
end
end

function tf = has_channel(model, s, mech)
tf = isfield(model.sec(s).channels, mech) && model.sec(s).channels.(mech).available;
end

function print_channel_summary(model)
fprintf('\n=== %s channel summary by compartment ===\n', model.name);
for s = 1:numel(model.sections)
    fprintf('%s:\n', model.sections{s});
    mechs = fieldnames(model.sec(s).channels);
    for k = 1:numel(mechs)
        mech = mechs{k};
        ch = model.sec(s).channels.(mech);
        fprintf('  %-6s (%-24s) gbar=%-10.6g  param=%s\n', ...
            mech, channel_title(mech), ch.gbar, ch.param_name);
    end
end
fprintf('\n');
end

function kin = mechanism_kinetics(mech, V, celsius)
% Returns gate opening inf, tau (ms), and P_open for known mechanisms.
% For mechanisms without kinetic equations here, marks unavailable.

kin.available = true;
kin.display_name = channel_title(mech);

switch mech
    case 'nax'
        [m_inf, tau_m, h_inf, tau_h] = nax_like_rates(V, celsius, 0);
        kin.gate_names = {'m','h'};
        kin.open_inf = [m_inf; h_inf];
        kin.tau_ms = [tau_m; tau_h];
        kin.p_open = m_inf.^3 .* h_inf;

    case 'na3'
        [m_inf, tau_m, h_inf, tau_h, s_inf, tau_s] = na3_rates(V, celsius, 0);
        kin.gate_names = {'m','h','s'};
        kin.open_inf = [m_inf; h_inf; s_inf];
        kin.tau_ms = [tau_m; tau_h; tau_s];
        kin.p_open = m_inf.^3 .* h_inf .* s_inf;

    case {'kdr','kdrb'}
        [n_inf, tau_n] = kd_family_rates(V, celsius, 13, -3, 0.7, 0.02, 2, 0);
        kin.gate_names = {'n'};
        kin.open_inf = n_inf;
        kin.tau_ms = tau_n;
        kin.p_open = n_inf;

    case 'kdb'
        [n_inf, tau_n] = kd_family_rates(V, celsius, -33, 3, 0.7, 0.005, 2, 0);
        kin.gate_names = {'n'};
        kin.open_inf = n_inf;
        kin.tau_ms = tau_n;
        kin.p_open = n_inf;

    case 'kap'
        [n_inf, tau_n, l_inf, tau_l] = ka_rates(V, celsius, 11, -56, 0.05, 0.05, -1.5, 3, 0.55, 1, 2, 0.1, -1, -40, 5, 1);
        kin.gate_names = {'n','l'};
        kin.open_inf = [n_inf; l_inf];
        kin.tau_ms = [tau_n; tau_l];
        kin.p_open = n_inf .* l_inf;

    case 'kad'
        [n_inf, tau_n, l_inf, tau_l] = ka_rates(V, celsius, -1, -56, 0.1, 0.05, -1.8, 3, 0.39, 1, 2, 0.2, -1, -40, 5, 1);
        kin.gate_names = {'n','l'};
        kin.open_inf = [n_inf; l_inf];
        kin.tau_ms = [tau_n; tau_l];
        kin.p_open = n_inf .* l_inf;

    case 'kmb'
        [m_inf, tau_m] = kmb_rates(V, celsius, -40, -10, -42, 0.003, 7, 0.4, 60, 0);
        kin.gate_names = {'m'};
        kin.open_inf = m_inf;
        kin.tau_ms = tau_m;
        kin.p_open = m_inf;

    case 'hd'
        [l_inf, tau_l] = hd_rates(V, celsius, -81, -8, -75, 0.011, 2.2, 0.4, 1, 4.5);
        kin.gate_names = {'l'};
        kin.open_inf = l_inf;
        kin.tau_ms = tau_l;
        kin.p_open = l_inf;

    case 'cal'
        [m_inf, tau_m] = cal_rates(V, celsius, 5, 0.2, 0.1, 2, 4, 0.1);
        kin.gate_names = {'m'};
        kin.open_inf = m_inf;
        kin.tau_ms = tau_m;
        h2 = 0.001 / (0.001 + 50e-6);
        kin.p_open = m_inf.^2 .* h2;

    case 'can'
        [m_inf, tau_m, h_inf, tau_h] = can_rates(V, celsius, 5, 0.2, 3, 0.03, 2, -14, 0.1);
        kin.gate_names = {'m','h'};
        kin.open_inf = [m_inf; h_inf];
        kin.tau_ms = [tau_m; tau_h];
        h2 = 0.001 / (0.001 + 50e-6);
        kin.p_open = m_inf.^2 .* h_inf .* h2;

    case 'cat'
        [m_inf, tau_m, h_inf, tau_h] = cat_rates(V, celsius, 5, 0.2, 10, 0.04, 2, -28, 0.1, 0.015, 3.5, -75, 0.6);
        kin.gate_names = {'m','h'};
        kin.open_inf = [m_inf; h_inf];
        kin.tau_ms = [tau_m; tau_h];
        kin.p_open = m_inf.^2 .* h_inf;

    case 'kca'
        [m_inf, tau_m] = kca_rates(celsius, 50e-6, 0.03, 0.00035, 0.5);
        m_inf = m_inf .* ones(size(V));
        tau_m = tau_m .* ones(size(V));
        kin.gate_names = {'m'};
        kin.open_inf = m_inf;
        kin.tau_ms = tau_m;
        kin.p_open = m_inf.^3;

    case 'cagk'
        [o_inf, tau_o] = cagk_rates(V, celsius, 50e-6, 0.84, 1, 0.48e-3, 0.13e-6, 0.28, 0.48);
        kin.gate_names = {'o'};
        kin.open_inf = o_inf;
        kin.tau_ms = tau_o;
        kin.p_open = o_inf;

    otherwise
        kin.available = false;
        kin.gate_names = {'n/a'};
        kin.open_inf = zeros(1, numel(V));
        kin.tau_ms = nan(1, numel(V));
        kin.p_open = zeros(1, numel(V));
end
end

function t = channel_title(mech)
switch mech
    case 'nax',  t = 'nax (Na fast)';
    case 'na3',  t = 'na3 (Na fast + slow inact)';
    case 'kdr',  t = 'kdr (delayed rectifier K)';
    case 'kdrb', t = 'kdrb (delayed rectifier K)';
    case 'kdb',  t = 'kdb (D-type K)';
    case 'kap',  t = 'kap (A-type K proximal)';
    case 'kad',  t = 'kad (A-type K distal)';
    case 'kmb',  t = 'kmb (M-type K)';
    case 'hd',   t = 'hd (Ih / HCN-like)';
    case 'cal',  t = 'cal (L-type Ca)';
    case 'can',  t = 'can (N-type Ca)';
    case 'cat',  t = 'cat (T-type Ca)';
    case 'kca',  t = 'kca (Ca-dependent K)';
    case 'cagk', t = 'cagk (Ca-activated K BK-like)';
    otherwise,   t = mech;
end
end

function style = gate_line_style(i)
styles = {'-','--',':','-.'};
style = styles{mod(i-1, numel(styles))+1};
end

% ===== Channel-specific rate helpers =====
function [m_inf, tau_m, h_inf, tau_h] = nax_like_rates(V, celsius, sh)
qt = 2.^((celsius-24)./10);
tha = -30; qa = 7.2; Ra = 0.4; Rb = 0.124;
thi1 = -45; thi2 = -45; qd = 1.5; qg = 1.5;
thinf = -50; qinf = 4; mmin = 0.02; hmin = 0.5;
a_m = trap0_fun(V, tha+sh, Ra, qa);
b_m = trap0_fun(-V, -tha-sh, Rb, qa);
den_m = a_m + b_m;
m_inf = a_m ./ den_m;
tau_m = max(1 ./ den_m ./ qt, mmin);
a_h = trap0_fun(V, thi1+sh, 0.03, qd);
b_h = trap0_fun(-V, -thi2-sh, 0.01, qg);
den_h = a_h + b_h;
h_inf = 1 ./ (1 + exp((V - thinf - sh)./qinf));
tau_h = max(1 ./ den_h ./ qt, hmin);
end

function [m_inf, tau_m, h_inf, tau_h, s_inf, tau_s] = na3_rates(V, celsius, sh)
[m_inf, tau_m, h_inf, tau_h] = nax_like_rates(V, celsius, sh);
vhalfs = -60; zetas = 12; gms = 0.2; a0s = 0.0003; smin = 10;
vvh = -58; vvs = 2; ar = 1;
alpv = 1 ./ (1 + exp((V - vvh - sh)./vvs));
alps = exp(1e-3 .* zetas .* (V - vhalfs - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
bets = exp(1e-3 .* zetas .* gms .* (V - vhalfs - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
s_inf = alpv + ar .* (1 - alpv);
tau_s = max(bets ./ (a0s .* (1 + alps)), smin);
end

function [ninf, taun] = kd_family_rates(V, celsius, vhalfn, zetan, gmn, a0n, nmin, sh)
alpn = exp(1e-3 .* zetan .* (V - vhalfn - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
betn = exp(1e-3 .* zetan .* gmn .* (V - vhalfn - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
ninf = 1 ./ (1 + alpn);
taun = max(betn ./ (a0n .* (1 + alpn)), nmin);
end

function [ninf, taun, linf, taul] = ka_rates(V, celsius, vhalfn, vhalfl, a0n, a0l, zetan, zetal, gmn, gml, lmin, nmin, pw, tq, qq, qtl)
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
taul = max(taul, lmin./qtl);
end

function [m_inf, tau_m] = kmb_rates(V, celsius, vhalfl, kl, vhalft, a0t, zetat, gmt, b0, sh)
infv = 1 ./ (1 + exp((V - vhalfl - sh)./kl));
a = exp(0.0378 .* zetat .* (V - vhalft - sh));
b = exp(0.0378 .* zetat .* gmt .* (V - vhalft - sh));
tau = b0 + b ./ (a0t .* (1 + a));
% mod has no q10 on tau final equation
m_inf = infv;
tau_m = tau;
end

function [l_inf, tau_l] = hd_rates(V, celsius, vhalfl, kl, vhalft, a0t, zetat, gmt, qtl, q10)
qt = q10.^((celsius-33)./10);
a = exp(0.0378 .* zetat .* (V - vhalft));
b = exp(0.0378 .* zetat .* gmt .* (V - vhalft));
l_inf = 1 ./ (1 + exp(-(V - vhalfl)./kl));
tau_l = b ./ (qtl .* qt .* a0t .* (1 + a));
end

function [m_inf, tau_m] = cal_rates(V, celsius, q10, mmin, a0m, zetam, vhalfm, gmm)
qt = q10.^((celsius-25)./10);
a = 15.69 .* (-V + 81.5) ./ (exp((-V + 81.5)./10) - 1);
b = 0.29 .* exp(-V./10.86);
m_inf = a ./ (a + b);
alpmt = exp(0.0378 .* zetam .* (V - vhalfm));
betmt = exp(0.0378 .* zetam .* gmm .* (V - vhalfm));
tau_m = betmt ./ (qt .* a0m .* (1 + alpmt));
tau_m = max(tau_m, mmin./qt);
end

function [m_inf, tau_m, h_inf, tau_h] = can_rates(V, celsius, q10, mmin, hmin, a0m, zetam, vhalfm, gmm)
qt = q10.^((celsius-25)./10);
a = 0.1967 .* (-V + 19.88) ./ (exp((-V + 19.88)./10) - 1);
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

function [m_inf, tau_m, h_inf, tau_h] = cat_rates(V, celsius, q10, mmin, hmin, a0m, zetam, vhalfm, gmm, a0h, zetah, vhalfh, gmh)
qt = q10.^((celsius-25)./10);
a = 0.2 .* (-V + 19.26) ./ (exp((-V + 19.26)./10.0) - 1.0);
b = 0.009 .* exp(-V ./ 22.03);
m_inf = a ./ (a + b);
alpmt = exp(0.0378 .* zetam .* (V - vhalfm));
betmt = exp(0.0378 .* zetam .* gmm .* (V - vhalfm));
tau_m = betmt ./ (qt .* a0m .* (1 + alpmt));
tau_m = max(tau_m, mmin./qt);
ah = exp(0.0378 .* zetah .* (V - vhalfh));
bh = exp(0.0378 .* zetah .* gmh .* (V - vhalfh));
h_inf = 1 ./ (1 + ah);
tau_h = bh ./ (qt .* a0h .* (1 + ah));
tau_h = max(tau_h, hmin./qt);
end

function [m_inf, tau_m] = kca_rates(celsius, cai, beta, cac, taumin)
tadj = 3.^((celsius-22)./10);
car = (cai./cac).^4;
m_inf = car ./ (1 + car);
tau_m = 1 ./ beta ./ (1 + car) ./ tadj;
tau_m = max(tau_m, taumin);
end

function [o_inf, tau_o] = cagk_rates(V, celsius, cai, d1, d2, k1, k2, abar, bbar)
exp1a = k1 .* exp(-2 .* d1 .* 96485 .* V ./ (8.313424 .* (273.15 + celsius)));
exp1b = k2 .* exp(-2 .* d2 .* 96485 .* V ./ (8.313424 .* (273.15 + celsius)));
a = cai .* abar ./ (cai + exp1a);
b = bbar ./ (1 + cai ./ exp1b);
tau_o = 1 ./ (a + b);
o_inf = a .* tau_o;
end

function y = trap0_fun(v, th, a, q)
x = v - th;
y = a .* x ./ (1 - exp(-x./q));
small = abs(x) <= 1e-6;
y(small) = a .* q;
end
