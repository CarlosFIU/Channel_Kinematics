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
% 4) HH/channel dynamics comparison for THIS repository (PC + PVBC)
% ------------------------------
% Source of truth used below:
%  - Conductance values: PC_dynamics_params_sonata.json / PVBC_dynamics_params_sonata.json
%  - Kinetic equations: modfiles/naxn.mod, kdrca1.mod, na3n.mod, kdrbca1.mod, kdb.mod
%
% This section adds both:
%   (a) gate steady-states + time constants (tau)
%   (b) channel opening probabilities (P_open)

V = linspace(-100, 50, 600);

hh_pc = hh_channel_kinematics_pc_params();
hh_pv = hh_channel_kinematics_pvbc_params();

g_pc = hh_compute_gating(V, hh_pc);
g_pv = hh_compute_gating(V, hh_pv);

fprintf('PC   m_inf range:   [%g, %g]\n', min(g_pc.m_inf), max(g_pc.m_inf));
fprintf('PC   h_inf range:   [%g, %g]\n', min(g_pc.h_inf), max(g_pc.h_inf));
fprintf('PC   n_inf range:   [%g, %g]\n', min(g_pc.n_inf), max(g_pc.n_inf));
fprintf('PVBC m_inf range:   [%g, %g]\n', min(g_pv.m_inf), max(g_pv.m_inf));
fprintf('PVBC h_inf range:   [%g, %g]\n', min(g_pv.h_inf), max(g_pv.h_inf));
fprintf('PVBC n_inf (kdrb):  [%g, %g]\n', min(g_pv.n_inf), max(g_pv.n_inf));
fprintf('PVBC n_inf (kdb):   [%g, %g]\n', min(g_pv.n_kdb_inf), max(g_pv.n_kdb_inf));

% ============================================================
% A) Gate steady-state curves
% ============================================================
figure('Name','PC vs PVBC gate steady-states','Color','w');

subplot(2,3,1);
plot(V, g_pc.m_inf, 'LineWidth', 1.5); hold on;
plot(V, g_pv.m_inf, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('m_\infty'); title('Na activation m_\infty'); grid on;
legend({hh_pc.name, hh_pv.name}, 'Location','best');

subplot(2,3,2);
plot(V, g_pc.h_inf, 'LineWidth', 1.5); hold on;
plot(V, g_pv.h_inf, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('h_\infty'); title('Na inactivation h_\infty'); grid on;
legend({hh_pc.name, hh_pv.name}, 'Location','best');

subplot(2,3,3);
plot(V, g_pc.n_inf, 'LineWidth', 1.5); hold on;
plot(V, g_pv.n_inf, '--', 'LineWidth', 1.5);
plot(V, g_pv.n_kdb_inf, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('n_\infty'); title('K activation'); grid on;
legend({'PC KDR','PVBC KDRb','PVBC KDB'}, 'Location','best');

subplot(2,3,4);
plot(V, g_pc.tau_m_ms, 'LineWidth', 1.5); hold on;
plot(V, g_pv.tau_m_ms, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('\tau_m (ms)'); title('Na \tau_m'); grid on;
legend({hh_pc.name, hh_pv.name}, 'Location','best');

subplot(2,3,5);
plot(V, g_pc.tau_h_ms, 'LineWidth', 1.5); hold on;
plot(V, g_pv.tau_h_ms, '--', 'LineWidth', 1.5);
plot(V, g_pv.tau_s_ms, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('tau (ms)'); title('Na \tau_h (and PVBC \tau_s)'); grid on;
legend({'PC \tau_h','PVBC \tau_h','PVBC \tau_s'}, 'Location','best');

subplot(2,3,6);
plot(V, g_pc.tau_n_ms, 'LineWidth', 1.5); hold on;
plot(V, g_pv.tau_n_ms, '--', 'LineWidth', 1.5);
plot(V, g_pv.tau_n_kdb_ms, ':', 'LineWidth', 1.8);
xlabel('V (mV)'); ylabel('tau (ms)'); title('K \tau_n'); grid on;
legend({'PC KDR \tau_n','PVBC KDRb \tau_n','PVBC KDB \tau_n'}, 'Location','best');

set(gcf,'Position',[100 80 1300 650]);

% ============================================================
% B) Channel opening probabilities P_open(V)
% ============================================================
figure('Name','PC vs PVBC channel opening probabilities','Color','w');

subplot(1,3,1);
plot(V, g_pc.p_open_Na, 'LineWidth', 1.5); hold on;
plot(V, g_pv.p_open_Na, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('P_{open}'); title('Na P_{open}'); grid on;
legend({'PC: m^3h','PVBC: m^3hs'}, 'Location','best');

subplot(1,3,2);
plot(V, g_pc.p_open_Kdr, 'LineWidth', 1.5); hold on;
plot(V, g_pv.p_open_Kdrb, '--', 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('P_{open}'); title('KDR/KDRb P_{open}'); grid on;
legend({'PC KDR (n)','PVBC KDRb (n)'}, 'Location','best');

subplot(1,3,3);
plot(V, g_pv.p_open_Kdb, 'LineWidth', 1.5);
xlabel('V (mV)'); ylabel('P_{open}'); title('PVBC KDB P_{open}'); grid on;
legend({'PVBC KDB (n)'}, 'Location','best');

set(gcf,'Position',[120 180 1200 330]);

% ============================================================
% C) Print JSON conductance values used by the P_open curves
% ============================================================
fprintf('\n--- Soma conductances from SONATA dynamics JSON ---\n');
fprintf('PC   gbar_nax      = %.8g S/cm^2\n', hh_pc.gbar_nax);
fprintf('PC   gkdrbar_kdr   = %.8g S/cm^2\n', hh_pc.gbar_kdr);
fprintf('PVBC gbar_na3      = %.8g S/cm^2\n', hh_pv.gbar_na3);
fprintf('PVBC gkdrbar_kdrb  = %.8g S/cm^2\n', hh_pv.gbar_kdrb);
fprintf('PVBC gkdbar_kdb    = %.8g S/cm^2\n', hh_pv.gbar_kdb);

% ============================================================
% D) Synaptic parameter table printout (unchanged)
% ============================================================
disp('--- Synaptic model parameter table ---');
for k = 1:numel(models)
    m = models(k);
    fprintf('%s\n', m.name);
    fprintf('  EXC: rise=%.3g ms, decay=%.3g ms\n', m.exc.tau_rise_ms, m.exc.tau_decay_ms);
    fprintf('  INH: rise=%.3g ms, decay=%.3g ms\n', m.inh.tau_rise_ms, m.inh.tau_decay_ms);
    if isfield(m,'tau_gsyn_ms') && ~isempty(m.tau_gsyn_ms) && ~isnan(m.tau_gsyn_ms)
        fprintf('  g_syn decay tau (both): %.3g ms\n', m.tau_gsyn_ms);
    end
    if isfield(m.exc,'gmax_nS') && ~isempty(m.exc.gmax_nS) && ~isnan(m.exc.gmax_nS)
        fprintf('  gmax exc: %.3g nS, gmax inh: %.3g nS\n', m.exc.gmax_nS, m.inh.gmax_nS);
    end
    fprintf('  Notes: %s\n\n', m.notes);
end

%% -----------------------------
% 5) Paper-level channel/synapse tables (five-paper list from review doc)
% ------------------------------
[paper_tbl, channel_tbl, synapse_tbl] = swr_build_paper_tables();

out_tbl_dir = fullfile(pwd, 'exports_model_tables');
if ~exist(out_tbl_dir, 'dir')
    mkdir(out_tbl_dir);
end

writetable(paper_tbl,   fullfile(out_tbl_dir, 'paper_index.csv'));
writetable(channel_tbl, fullfile(out_tbl_dir, 'channel_kinematics_table.csv'));
writetable(synapse_tbl, fullfile(out_tbl_dir, 'synapse_connection_table.csv'));

fprintf('\n--- Exported paper-level tables ---\n');
fprintf('Paper index:      %s\n', fullfile(out_tbl_dir, 'paper_index.csv'));
fprintf('Channel table:    %s\n', fullfile(out_tbl_dir, 'channel_kinematics_table.csv'));
fprintf('Synapse table:    %s\n', fullfile(out_tbl_dir, 'synapse_connection_table.csv'));

disp(' ');
disp('Preview: channel table (first 10 rows)');
disp(channel_tbl(1:min(10,height(channel_tbl)), :));

disp(' ');
disp('Preview: synapse table (first 10 rows)');
disp(synapse_tbl(1:min(10,height(synapse_tbl)), :));


%% -----------------------------
% 7) Targeted PC/PVBC kinematics with explicit Bezaire CA1 equations
% ------------------------------
swr_plot_pc_pvbc_targeted_kinematics(V);
