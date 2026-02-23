function hh_channel_kinematics_compartment_compare()
% Compare channel kinetics across compartments for PC and PVBC SONATA JSONs.
% Compartments:
%   - PC: soma, dend, apic
%   - PVBC: soma, dend
%
% For each compartment, plots gate steady-states (m_inf/h_inf/n_inf),
% gate time constants (tau_m/tau_h/tau_n), and conductance-weighted
% channel opening proxies (gbar * P_open).

close all;
V = linspace(-100, 50, 600);

pc_sections = {'soma','dend','apic'};
pv_sections = {'soma','dend'};

pc_data = struct([]);
for i = 1:numel(pc_sections)
    hh = hh_channel_kinematics_pc_params([], pc_sections{i});
    g = hh_compute_gating(V, hh);
    pc_data(i).section = pc_sections{i};
    pc_data(i).hh = hh;
    pc_data(i).g = g;
end

pv_data = struct([]);
for i = 1:numel(pv_sections)
    hh = hh_channel_kinematics_pvbc_params([], pv_sections{i});
    g = hh_compute_gating(V, hh);
    pv_data(i).section = pv_sections{i};
    pv_data(i).hh = hh;
    pv_data(i).g = g;
    pv_data(i).has_na = ~isnan(hh.gbar_na3) && hh.gbar_na3 > 0;
    pv_data(i).has_kdrb = ~isnan(hh.gbar_kdrb) && hh.gbar_kdrb > 0;
end

plot_pc_compartments(V, pc_data);
plot_pv_compartments(V, pv_data);
print_conductance_summary(pc_data, pv_data);

end

function plot_pc_compartments(V, data)
figure('Name','PC channel kinetics across compartments','Color','w');
cols = lines(numel(data));

subplot(2,3,1); hold on;
for i = 1:numel(data), plot(V, data(i).g.m_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('m_\infty'); title('PC Na activation m_\infty'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,2); hold on;
for i = 1:numel(data), plot(V, data(i).g.h_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('h_\infty'); title('PC Na inactivation h_\infty'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,3); hold on;
for i = 1:numel(data), plot(V, data(i).g.n_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('n_\infty'); title('PC KDR activation n_\infty'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,4); hold on;
for i = 1:numel(data), plot(V, data(i).g.tau_m_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('\tau_m (ms)'); title('PC Na \tau_m'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,5); hold on;
for i = 1:numel(data), plot(V, data(i).g.tau_h_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('\tau_h (ms)'); title('PC Na \tau_h'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,6); hold on;
for i = 1:numel(data), plot(V, data(i).g.tau_n_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('\tau_n (ms)'); title('PC KDR \tau_n'); grid on;
legend({data.section}, 'Location', 'best');

set(gcf, 'Position', [90 70 1300 700]);

figure('Name','PC conductance-weighted opening across compartments','Color','w');
subplot(1,2,1); hold on;
for i = 1:numel(data)
    y = data(i).hh.gbar_nax .* data(i).g.p_open_Na;
    plot(V, y, 'LineWidth', 1.6, 'Color', cols(i,:));
end
xlabel('V (mV)'); ylabel('g_{Na}*P_{open} (S/cm^2)');
title('PC Na (nax) effective opening'); grid on;
legend({data.section}, 'Location', 'best');

subplot(1,2,2); hold on;
for i = 1:numel(data)
    y = data(i).hh.gbar_kdr .* data(i).g.p_open_Kdr;
    plot(V, y, 'LineWidth', 1.6, 'Color', cols(i,:));
end
xlabel('V (mV)'); ylabel('g_{KDR}*P_{open} (S/cm^2)');
title('PC KDR (kdrca1) effective opening'); grid on;
legend({data.section}, 'Location', 'best');
set(gcf, 'Position', [120 220 1100 380]);
end

function plot_pv_compartments(V, data)
figure('Name','PVBC channel kinetics across compartments','Color','w');
cols = lines(numel(data));

subplot(2,3,1); hold on;
for i = 1:numel(data), if data(i).has_na, plot(V, data(i).g.m_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end, end
xlabel('V (mV)'); ylabel('m_\infty'); title('PVBC Na activation m_\infty'); grid on;
legend(channel_sections(data, 'has_na'), 'Location', 'best');

subplot(2,3,2); hold on;
for i = 1:numel(data), if data(i).has_na, plot(V, data(i).g.h_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end, end
xlabel('V (mV)'); ylabel('h_\infty'); title('PVBC Na inactivation h_\infty'); grid on;
legend(channel_sections(data, 'has_na'), 'Location', 'best');

subplot(2,3,3); hold on;
for i = 1:numel(data), plot(V, data(i).g.n_inf, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('n_\infty'); title('PVBC KDRb activation n_\infty'); grid on;
legend({data.section}, 'Location', 'best');

subplot(2,3,4); hold on;
for i = 1:numel(data), if data(i).has_na, plot(V, data(i).g.tau_m_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end, end
xlabel('V (mV)'); ylabel('\tau_m (ms)'); title('PVBC Na \tau_m'); grid on;
legend(channel_sections(data, 'has_na'), 'Location', 'best');

subplot(2,3,5); hold on;
for i = 1:numel(data), if data(i).has_na, plot(V, data(i).g.tau_h_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end, end
xlabel('V (mV)'); ylabel('\tau_h (ms)'); title('PVBC Na \tau_h'); grid on;
legend(channel_sections(data, 'has_na'), 'Location', 'best');

subplot(2,3,6); hold on;
for i = 1:numel(data), plot(V, data(i).g.tau_n_ms, 'LineWidth', 1.5, 'Color', cols(i,:)); end
xlabel('V (mV)'); ylabel('\tau_n (ms)'); title('PVBC KDRb \tau_n'); grid on;
legend({data.section}, 'Location', 'best');

set(gcf, 'Position', [100 80 1300 700]);

figure('Name','PVBC conductance-weighted opening across compartments','Color','w');
subplot(1,2,1); hold on;
for i = 1:numel(data)
    if data(i).has_na
        y = data(i).hh.gbar_na3 .* data(i).g.p_open_Na;
        plot(V, y, 'LineWidth', 1.6, 'Color', cols(i,:));
    end
end
xlabel('V (mV)'); ylabel('g_{Na}*P_{open} (S/cm^2)');
title('PVBC Na (na3n) effective opening'); grid on;
legend(channel_sections(data, 'has_na'), 'Location', 'best');

subplot(1,2,2); hold on;
for i = 1:numel(data)
    y = data(i).hh.gbar_kdrb .* data(i).g.p_open_Kdrb;
    plot(V, y, 'LineWidth', 1.6, 'Color', cols(i,:));
end
xlabel('V (mV)'); ylabel('g_{KDRb}*P_{open} (S/cm^2)');
title('PVBC KDRb (kdrbca1) effective opening'); grid on;
legend({data.section}, 'Location', 'best');
set(gcf, 'Position', [140 220 1100 380]);
end

function print_conductance_summary(pc_data, pv_data)
fprintf('\n--- PC conductances by compartment (SONATA JSON) ---\n');
for i = 1:numel(pc_data)
    fprintf('%s: gbar_nax=%.8g, gkdrbar_kdr=%.8g\n', ...
        pc_data(i).section, pc_data(i).hh.gbar_nax, pc_data(i).hh.gbar_kdr);
end

fprintf('\n--- PVBC conductances by compartment (SONATA JSON) ---\n');
for i = 1:numel(pv_data)
    fprintf('%s: gbar_na3=%.8g, gkdrbar_kdrb=%.8g\n', ...
        pv_data(i).section, pv_data(i).hh.gbar_na3, pv_data(i).hh.gbar_kdrb);
end

fprintf('\nNote: PVBC gkdbar_kdb is soma-only in this JSON and is not compared across dendrite.\n');
end


function labels = channel_sections(data, flag_name)
labels = {};
for i = 1:numel(data)
    if data(i).(flag_name)
        labels{end+1} = data(i).section; %#ok<AGROW>
    end
end
if isempty(labels)
    labels = {'none'};
end
end
