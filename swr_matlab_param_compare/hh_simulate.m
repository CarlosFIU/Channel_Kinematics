function out = hh_simulate(Iinj_uAcm2, t_ms, hh, V0_mV)
% hh_simulate: simulate classic HH squid axon under a current clamp.
% Iinj_uAcm2 can be scalar or vector (same length as t_ms).
% Uses forward Euler for simplicity; switch to ode15s for stiff models.

dt = t_ms(2) - t_ms(1);
if isscalar(Iinj_uAcm2)
    Iinj_uAcm2 = Iinj_uAcm2 * ones(size(t_ms));
end

V = zeros(size(t_ms)); V(1) = V0_mV;
% Initialize gating at steady-state
g0 = hh_compute_gating(V0_mV, hh);
m = g0.m_inf; h = g0.h_inf; n = g0.n_inf;

for k = 2:numel(t_ms)
    Vk = V(k-1);

    % Gating updates
    am = hh.alpha_m(Vk); bm = hh.beta_m(Vk);
    ah = hh.alpha_h(Vk); bh = hh.beta_h(Vk);
    an = hh.alpha_n(Vk); bn = hh.beta_n(Vk);

    m = m + dt * (am*(1-m) - bm*m);
    h = h + dt * (ah*(1-h) - bh*h);
    n = n + dt * (an*(1-n) - bn*n);

    % Currents (uA/cm^2)
    INa = hh.gNa_mScm2 * (m^3) * h * (Vk - hh.ENa_mV);
    IK  = hh.gK_mScm2  * (n^4)     * (Vk - hh.EK_mV);
    IL  = hh.gL_mScm2            * (Vk - hh.EL_mV);

    dV = (Iinj_uAcm2(k-1) - INa - IK - IL) / hh.Cm_uFcm2;
    V(k) = Vk + dt*dV;
end

out.t_ms = t_ms;
out.V_mV = V;
end
