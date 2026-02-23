% swr_param_compare.m
% Compare synaptic time constants across several SWR network models (from papers)
% and provide a framework to compare HH-style gating variables when kinetics are known.
%
% Papers with open-access parameter statements used here:
%  - Evangelista et al., J Neurosci 2020 (PMC7548694): synaptic tau_s (exc/inh) in ms
%  - Jahnke et al., J Neurosci 2015: alpha-kernel tau1/tau2 for exc/inh + conductance decay tau
%  - Canakci et al., PLoS ONE 2017: Exp2Syn-like rise/decay taus for AMPA/GABA + gmax
%
% NOTE: Many papers do NOT list ion-channel (HH) kinetic equations in the manuscript;
% those are typically defined in the accompanying simulator code (NEURON .mod/.hoc, Brian, etc).
% The HH comparison portion below is a scaffold: add models as you obtain kinetics/code.

clear; clc; close all;

%% -----------------------------
% 1) Synaptic time-constant database
% ------------------------------
% Convention here:
%   - For models that specify a single synaptic time constant tau_s: store as tau_decay
%   - For Exp2Syn/alpha kernel models: store tau_rise and tau_decay
% Units: milliseconds (ms)

models = struct([]);

% ------------------
% Evangelista 2020 (conductance-based LIF; time constants given as tau_s)
% Exc synapse: tau_s = 3 ms; Inh synapse: tau_s = 2 ms
% ------------------
models(end+1).name = 'Evangelista2020 (LIF)';
models(end).E_exc_mV = 0;
models(end).E_inh_mV = -70;
models(end).exc.tau_rise_ms = NaN;
models(end).exc.tau_decay_ms = 3;
models(end).inh.tau_rise_ms = NaN;
models(end).inh.tau_decay_ms = 2;
models(end).notes = 'Conductance-based LIF; manuscript reports synaptic tau_s.';

% ------------------
% Jahnke 2015
% They use an alpha kernel alpha(t) = exp(-t/tau2) - exp(-t/tau1)
% for excitatory and inhibitory inputs:
%   Exc kernel: tau1=0.3 ms, tau2=0.6 ms
%   Inh kernel: tau1=0.3 ms, tau2=2.0 ms
% They also specify a conductance decay time constant tau_ex=tau_in=10 ms for g_syn.
% Here we store kernel rise/decay as tau1/tau2 and include tau_gsyn as an extra field.
% ------------------
models(end+1).name = 'Jahnke2015 (LIF+alpha kernel)';
models(end).exc.tau_rise_ms  = 0.3;
models(end).exc.tau_decay_ms = 0.6;
models(end).inh.tau_rise_ms  = 0.3;
models(end).inh.tau_decay_ms = 2.0;
models(end).tau_gsyn_ms = 10;   % both exc/inh conductance decay constant in their g_syn ODE
models(end).notes = 'Alpha kernel + separate g_syn decay tau.';

% ------------------
% Canakci 2017 (NEURON; Exp2Syn-like kinetics)
% Exc: tau_rise = 0.5 ms, tau_decay = 3.5 ms; gmax = 0.2 nS
% Inh: tau_rise = 0.4 ms, tau_decay = 5.0 ms; gmax = 2.0 nS
% ------------------
models(end+1).name = 'Canakci2017 (NEURON)';
models(end).exc.tau_rise_ms  = 0.5;
models(end).exc.tau_decay_ms = 3.5;
models(end).exc.gmax_nS      = 0.2;
models(end).inh.tau_rise_ms  = 0.4;
models(end).inh.tau_decay_ms = 5.0;
models(end).inh.gmax_nS      = 2.0;
models(end).notes = 'Exp2Syn-like rise/decay taus and gmax reported in paper.';

%% -----------------------------
% 2) Plot synaptic tau parameters
% ------------------------------
names = string({models.name});

exc_r = arrayfun(@(m) m.exc.tau_rise_ms, models);
exc_d = arrayfun(@(m) m.exc.tau_decay_ms, models);
inh_r = arrayfun(@(m) m.inh.tau_rise_ms, models);
inh_d = arrayfun(@(m) m.inh.tau_decay_ms, models);

figure('Name','Synaptic time constants (rise/decay)','Color','w');
subplot(1,2,1);
bar([exc_r(:) exc_d(:)], 'grouped');
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',25);
ylabel('Exc tau (ms)');
legend({'rise','decay'},'Location','northwest');
title('Excitatory synapse kinetics');

subplot(1,2,2);
bar([inh_r(:) inh_d(:)], 'grouped');
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',25);
ylabel('Inh tau (ms)');
legend({'rise','decay'},'Location','northwest');
title('Inhibitory synapse kinetics');

%% -----------------------------
% 3) Optional: plot alpha/Exp2Syn kernels for visual comparison
% ------------------------------
t = linspace(0, 20, 5001); % ms

figure('Name','Normalized synaptic kernels','Color','w');
hold on;

for k = 1:numel(models)
    m = models(k);
    % Build kernel based on availability
    if ~isnan(m.exc.tau_rise_ms) && ~isnan(m.exc.tau_decay_ms)
        g = exp(-t/m.exc.tau_decay_ms) - exp(-t/m.exc.tau_rise_ms);
        g(g<0)=0;
        g = g ./ max(g);
        plot(t, g, 'DisplayName', m.name + " (exc)");
    else
        % single-exponential decay
        g = exp(-t/m.exc.tau_decay_ms);
        g = g ./ max(g);
        plot(t, g, '--', 'DisplayName', m.name + " (exc, exp)");
    end
end
xlabel('t (ms)'); ylabel('Normalized kernel');
title('Excitatory kernels (shape comparison)');
legend('Location','northeastoutside'); grid on;
hold off;

%% -----------------------------
% 4) HH gating comparison scaffold (UPDATED: Squid vs Canakci Pyr vs Canakci PV)
% ------------------------------
% Changes:
%  - Compare 3 models:
%       (1) Squid HH
%       (2) Canakci/Fink CA1 pyramidal (ca1ina + ca1ikdr)  -> hh.kind='canakci_fink'
%       (3) Canakci/Fink PV/basket (nafbwb + kdrbwb)       -> hh.kind='canakci_pv'
%  - Plot steady-states (m_inf, h_inf, n_inf) in one figure
%  - Plot time constants for Na (tau_m, tau_h) and K (tau_n) in a separate figure
%  - Optional plots for Canakci pyramidal slow Na inactivation (i_inf and tau_i)
%  - PV m is instantaneous in Nafbwb -> tau_m is NaN (skip in tau_m panel)

V = linspace(-100, 50, 600);

% ----------------------------
% Model 1: Squid HH
% ----------------------------
hh1 = hh_squid_params();                 % ensure hh1.kind is 'squid1952' or 'alphabeta'
if ~isfield(hh1,'name'); hh1.name = 'Squid HH (1952)'; end
g1  = hh_compute_gating(V, hh1);

% ----------------------------
% Model 2: Canakci/Fink CA1 pyramidal
% ----------------------------
% Your file:
%   function hh = hh_canakci_params()
%       hh.kind = 'canakci_fink';
%       hh.name = 'CA1 (Canakci/Fink)';
%       hh.vi = -60; hh.ki = 0.8;
%   end
hh2 = hh_canakci_pyr_params();               % pyramidal
if ~isfield(hh2,'name'); hh2.name = 'Canakci/Fink CA1 pyramidal'; end
g2  = hh_compute_gating(V, hh2);

% ----------------------------
% Model 3: Canakci/Fink PV basket
% ----------------------------
% Requires you created:
%   function hh = hh_canakci_pv_params(celsius)
%       hh.kind='canakci_pv'; hh.celsius=celsius; hh.name=...
%   end
hh3 = hh_canakci_pv_params(34);          % PV/basket @ 34C
if ~isfield(hh3,'name'); hh3.name = 'Canakci/Fink PV basket'; end
g3  = hh_compute_gating(V, hh3);

% ----------------------------
% Sanity checks (optional)
% ----------------------------
fprintf('Squid      m_inf range: [%g, %g]\n', min(g1.m_inf), max(g1.m_inf));
fprintf('CanakciPYR m_inf range: [%g, %g]\n', min(g2.m_inf), max(g2.m_inf));
fprintf('CanakciPV  m_inf range: [%g, %g]\n', min(g3.m_inf), max(g3.m_inf));

fprintf('Squid      h_inf range: [%g, %g]\n', min(g1.h_inf), max(g1.h_inf));
fprintf('CanakciPYR h_inf range: [%g, %g]\n', min(g2.h_inf), max(g2.h_inf));
fprintf('CanakciPV  h_inf range: [%g, %g]\n', min(g3.h_inf), max(g3.h_inf));

fprintf('Squid      n_inf range: [%g, %g]\n', min(g1.n_inf), max(g1.n_inf));
fprintf('CanakciPYR n_inf range: [%g, %g]\n', min(g2.n_inf), max(g2.n_inf));
fprintf('CanakciPV  n_inf range: [%g, %g]\n', min(g3.n_inf), max(g3.n_inf));

% ============================================================
% A) Steady-state curves
% ============================================================
figure('Name','HH gating steady-states: Squid vs Canakci Pyr vs Canakci PV','Color','w');

subplot(1,3,1);
plot(V, g1.m_inf, 'LineWidth', 1.5); hold on;
plot(V, g2.m_inf, '--', 'LineWidth', 1.5);
plot(V, g3.m_inf, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('m_\infty'); title('Na activation'); grid on;
legend({hh1.name, hh2.name, hh3.name}, 'Location','best');

subplot(1,3,2);
plot(V, g1.h_inf, 'LineWidth', 1.5); hold on;
plot(V, g2.h_inf, '--', 'LineWidth', 1.5);
plot(V, g3.h_inf, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('h_\infty'); title('Na inactivation'); grid on;
legend({hh1.name, hh2.name, hh3.name}, 'Location','best');

subplot(1,3,3);
plot(V, g1.n_inf, 'LineWidth', 1.5); hold on;
plot(V, g2.n_inf, '--', 'LineWidth', 1.5);
plot(V, g3.n_inf, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('n_\infty'); title('K activation (KDR)'); grid on;
legend({hh1.name, hh2.name, hh3.name}, 'Location','best');

set(gcf,'Position',[100 100 1200 350]);

% ============================================================
% B) Time constants (Na + K)
% ============================================================
figure('Name','HH gating time constants: Squid vs Canakci Pyr vs Canakci PV','Color','w');

% ---- tau_m: PV m is instantaneous -> tau_m is NaN -> skip PV in this panel
subplot(1,3,1);
plot(V, g1.tau_m_ms, 'LineWidth', 1.5); hold on;
plot(V, g2.tau_m_ms, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('\tau_m (ms)'); title('Na \tau_m'); grid on;
legend({hh1.name, hh2.name}, 'Location','best');

% ---- tau_h: compare all three (PV has h state)
subplot(1,3,2);
plot(V, g1.tau_h_ms, 'LineWidth', 1.5); hold on;
plot(V, g2.tau_h_ms, '--', 'LineWidth', 1.5);
plot(V, g3.tau_h_ms, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('\tau_h (ms)'); title('Na \tau_h'); grid on;
legend({hh1.name, hh2.name, hh3.name}, 'Location','best');

% ---- tau_n: compare all three
subplot(1,3,3);
plot(V, g1.tau_n_ms, 'LineWidth', 1.5); hold on;
plot(V, g2.tau_n_ms, '--', 'LineWidth', 1.5);
plot(V, g3.tau_n_ms, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('\tau_n (ms)'); title('KDR \tau_n'); grid on;
legend({hh1.name, hh2.name, hh3.name}, 'Location','best');

set(gcf,'Position',[100 520 1200 350]);

% ============================================================
% C) Optional: Canakci PYRAMIDAL slow Na inactivation i-gate
% ============================================================
if isfield(g2,'i_inf')
    figure('Name','Canakci/Fink PYR: slow Na inactivation i-gate','Color','w');

    subplot(1,2,1);
    plot(V, g2.i_inf, 'LineWidth', 1.5); grid on;
    xlabel('V (mV)'); ylabel('i_\infty'); title(['Slow Na inactivation i_\infty: ' hh2.name]);

    subplot(1,2,2);
    if isfield(g2,'tau_i_ms')
        plot(V, g2.tau_i_ms, 'LineWidth', 1.5); grid on;
        xlabel('V (mV)'); ylabel('\tau_i (ms)'); title('Slow Na inactivation \tau_i');
        % set(gca,'YScale','log'); % optional
    else
        text(0.1,0.5,'tau_i_ms not found in g2','Units','normalized'); axis off;
    end
end

% ============================================================
% D) Synaptic parameter table printout (unchanged)
% ============================================================
disp('--- Synaptic model parameter table ---');
for k = 1:numel(models)
    m = models(k);
    fprintf('%s\n', m.name);
    fprintf('  EXC: rise=%.3g ms, decay=%.3g ms\n', m.exc.tau_rise_ms, m.exc.tau_decay_ms);
    fprintf('  INH: rise=%.3g ms, decay=%.3g ms\n', m.inh.tau_rise_ms, m.inh.tau_decay_ms);
    if isfield(m,'tau_gsyn_ms')
        fprintf('  g_syn decay tau (both): %.3g ms\n', m.tau_gsyn_ms);
    end
    if isfield(m.exc,'gmax_nS')
        fprintf('  gmax exc: %.3g nS, gmax inh: %.3g nS\n', m.exc.gmax_nS, m.inh.gmax_nS);
    end
    fprintf('  Notes: %s\n\n', m.notes);
end
