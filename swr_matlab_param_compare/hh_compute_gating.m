function g = hh_compute_gating(V, hh)
% hh_compute_gating
% Compute steady-state (x_inf) and time-constants (tau_x) for gating variables.
%
% Modes:
%   1) Alpha/Beta HH-style (classic HH) if:
%        hh.kind = 'alphabeta'  (or 'squid1952')
%      Requires fields:
%        hh.alpha_m, hh.beta_m, hh.alpha_h, hh.beta_h, hh.alpha_n, hh.beta_n
%
%   2) Canakci/Fink CA1 pyramidal kinetics if:
%        hh.kind = 'canakci_fink'
%      (Na: m,h,i; KDR: n)   ~ ca1ina.mod + ca1ikdr.mod style
%
%   3) Canakci/Fink basket/PV kinetics if:
%        hh.kind = 'canakci_pv'
%      (Na: m_inf instantaneous + h; KDR: n) ~ nafbwb.mod + kdrbwb.mod style
%      Optional: hh.celsius (default 34). Uses Q10 scaling as in the mod files.
%
% Inputs:
%   V  : vector (or scalar) membrane voltages in mV
%   hh : struct
%
% Output:
%   g : struct with fields:
%       m_inf, h_inf, n_inf, tau_m_ms, tau_h_ms, tau_n_ms
%       plus optional i_inf, tau_i_ms for canakci_fink

if ~isfield(hh,'kind') || isempty(hh.kind)
    hh.kind = 'alphabeta';
end

switch lower(hh.kind)

    % =========================================================
    % 1) Generic alpha/beta HH-style
    % =========================================================
    case {'alphabeta','squid1952'}
        am = hh.alpha_m(V); bm = hh.beta_m(V);
        ah = hh.alpha_h(V); bh = hh.beta_h(V);
        an = hh.alpha_n(V); bn = hh.beta_n(V);

        denom_m = am + bm;
        denom_h = ah + bh;
        denom_n = an + bn;

        g.m_inf = am ./ denom_m;
        g.h_inf = ah ./ denom_h;
        g.n_inf = an ./ denom_n;

        g.tau_m_ms = 1 ./ denom_m;
        g.tau_h_ms = 1 ./ denom_h;
        g.tau_n_ms = 1 ./ denom_n;

    % =========================================================
    % 2) Canakci/Fink CA1 pyramidal kinetics (stable)
    % =========================================================
    case 'canakci_fink'
        Vin = V;
        V = V(:);   % column

        % --- Na activation m (ca1ina.mod style) ---
        x = V + 30;
        a_m = 0.4   .* vtrap_1mexp(x, 7.2);     % x/(1-exp(-x/7.2))
        b_m = 0.124 .* vtrap_exp(x, 7.2);       % x/(exp(x/7.2)-1)
        denom_m = a_m + b_m;

        g.m_inf    = a_m ./ denom_m;
        g.tau_m_ms = 0.5 ./ denom_m;
        g.tau_m_ms = max(g.tau_m_ms, 0.02);

        % --- Na inactivation h ---
        x = V + 45;
        a_h = 0.03 .* vtrap_1mexp(x, 1.5);
        b_h = 0.01 .* vtrap_exp(x, 1.5);
        denom_h = a_h + b_h;

        g.h_inf    = 1 ./ (1 + exp((V + 50)./4));
        g.tau_h_ms = 0.5 ./ denom_h;
        g.tau_h_ms = max(g.tau_h_ms, 0.5);

        % --- Na slow inactivation i ---
        if ~isfield(hh,'vi'); hh.vi = -60; end
        if ~isfield(hh,'ki'); hh.ki = 0.8; end

        a_i = exp(0.45 .* (V + 66));
        b_i = exp(0.09 .* (V + 66));

        g.tau_i_ms = 3000 .* b_i ./ (1 + a_i);
        g.tau_i_ms = max(g.tau_i_ms, 10);

        g.i_inf = (1 + hh.ki .* exp((V - hh.vi)./2)) ./ (1 + exp((V - hh.vi)./2));

        % --- KDR activation n (ca1ikdr.mod style) ---
        a_n = exp(-0.11 .* (V - 13));
        b_n = exp(-0.08 .* (V - 13));

        g.n_inf    = 1 ./ (1 + a_n);
        g.tau_n_ms = 50 .* b_n ./ (1 + a_n);
        g.tau_n_ms = max(g.tau_n_ms, 2);

        % reshape back
        g.m_inf    = reshape(g.m_inf,    size(Vin));
        g.h_inf    = reshape(g.h_inf,    size(Vin));
        g.n_inf    = reshape(g.n_inf,    size(Vin));
        g.tau_m_ms = reshape(g.tau_m_ms, size(Vin));
        g.tau_h_ms = reshape(g.tau_h_ms, size(Vin));
        g.tau_n_ms = reshape(g.tau_n_ms, size(Vin));
        g.i_inf    = reshape(g.i_inf,    size(Vin));
        g.tau_i_ms = reshape(g.tau_i_ms, size(Vin));

    % =========================================================
    % 3) Canakci/Fink basket/PV kinetics (nafbwb + kdrbwb style)
    % =========================================================
    case 'canakci_pv'
        Vin = V;
        V = V(:);   % column

        % Temperature/Q10 (matches typical BWB mods: phi = 5, ref 27C)
        if ~isfield(hh,'celsius') || isempty(hh.celsius)
            hh.celsius = 34;
        end
        phih = 5; phin = 5;
        q10h = phih.^((hh.celsius - 27)./10);
        q10n = phin.^((hh.celsius - 27)./10);

        % -------- Nafbwb --------
        % m is INSTANTANEOUS in Nafbwb (no tau_m state); compute m_inf only
        am = fun3(V, -35, -0.1, -10);
        bm = fun1(V, -60,  4.0, -18);
        g.m_inf    = am ./ (am + bm);
        g.tau_m_ms = nan(size(V));   % correct for Nafbwb

        % h gate
        ah = fun1(V, -58, 0.07, -20);
        bh = fun2(V, -28, 1.00, -10);
        g.h_inf    = ah ./ (ah + bh);
        g.tau_h_ms = 1 ./ ((ah + bh) .* q10h);

        % -------- Kdrbwb --------
        an = fun3(V, -34, -0.01, -10);
        bn = fun1(V, -44,  0.125, -80);
        g.n_inf    = an ./ (an + bn);
        g.tau_n_ms = 1 ./ ((an + bn) .* q10n);

        % reshape back
        g.m_inf    = reshape(g.m_inf,    size(Vin));
        g.h_inf    = reshape(g.h_inf,    size(Vin));
        g.n_inf    = reshape(g.n_inf,    size(Vin));
        g.tau_m_ms = reshape(g.tau_m_ms, size(Vin));
        g.tau_h_ms = reshape(g.tau_h_ms, size(Vin));
        g.tau_n_ms = reshape(g.tau_n_ms, size(Vin));



    % =========================================================
    % 4) Channel_Kinematics PC soma kinetics (naxn + kdrca1)
    % =========================================================
    case 'ck_pc'
        Vin = V;
        V = V(:);

        if ~isfield(hh,'celsius') || isempty(hh.celsius)
            hh.celsius = 34;
        end
        qt = 2.^((hh.celsius-24)./10); % naxn.mod q10

        % ---- naxn.mod ----
        tha = -30; qa = 7.2; Ra = 0.4; Rb = 0.124;
        thi1 = -45; thi2 = -45; qd = 1.5; qg = 1.5;
        thinf = -50; qinf = 4;
        mmin = 0.02; hmin = 0.5;
        sh = get_field_or_default(hh,'sh_nax',0);

        a_m = trap0_fun(V, tha+sh, Ra, qa);
        b_m = trap0_fun(-V, -tha-sh, Rb, qa);
        denom_m = a_m + b_m;
        g.m_inf = a_m ./ denom_m;
        g.tau_m_ms = 1 ./ denom_m ./ qt;
        g.tau_m_ms = max(g.tau_m_ms, mmin);

        a_h = trap0_fun(V, thi1+sh, 0.03, qd);
        b_h = trap0_fun(-V, -thi2-sh, 0.01, qg);
        denom_h = a_h + b_h;
        g.h_inf = 1 ./ (1 + exp((V - thinf - sh)./qinf));
        g.tau_h_ms = 1 ./ denom_h ./ qt;
        g.tau_h_ms = max(g.tau_h_ms, hmin);

        % ---- kdrca1.mod (kdr) ----
        a0n = 0.02; zetan = -3; gmn = 0.7; vhalfn = 13; nmin = 2;
        qtn = 1.^((hh.celsius-24)./10); % q10=1 in kdrca1.mod

        alpn = exp(1e-3 .* zetan .* (V - vhalfn) .* 9.648e4 ./ (8.315 .* (273.16 + hh.celsius)));
        betn = exp(1e-3 .* zetan .* gmn .* (V - vhalfn) .* 9.648e4 ./ (8.315 .* (273.16 + hh.celsius)));
        g.n_inf = 1 ./ (1 + alpn);
        g.tau_n_ms = betn ./ (qtn .* a0n .* (1 + alpn));
        g.tau_n_ms = max(g.tau_n_ms, nmin);

        % channel opening probabilities and max conductances from JSON
        g.p_open_Na = g.m_inf.^3 .* g.h_inf;
        g.p_open_Kdr = g.n_inf;
        g.gbar_Na = get_field_or_default(hh,'gbar_nax',nan);
        g.gbar_Kdr = get_field_or_default(hh,'gbar_kdr',nan);

        g.m_inf    = reshape(g.m_inf,    size(Vin));
        g.h_inf    = reshape(g.h_inf,    size(Vin));
        g.n_inf    = reshape(g.n_inf,    size(Vin));
        g.tau_m_ms = reshape(g.tau_m_ms, size(Vin));
        g.tau_h_ms = reshape(g.tau_h_ms, size(Vin));
        g.tau_n_ms = reshape(g.tau_n_ms, size(Vin));
        g.p_open_Na = reshape(g.p_open_Na, size(Vin));
        g.p_open_Kdr = reshape(g.p_open_Kdr, size(Vin));

    % =========================================================
    % 5) Channel_Kinematics PVBC soma kinetics (na3n + kdrbca1 + kdb)
    % =========================================================
    case 'ck_pvbc'
        Vin = V;
        V = V(:);

        if ~isfield(hh,'celsius') || isempty(hh.celsius)
            hh.celsius = 34;
        end
        qt_na = 2.^((hh.celsius-24)./10); % na3n.mod q10

        % ---- na3n.mod ----
        tha = -30; qa = 7.2; Ra = 0.4; Rb = 0.124;
        thi1 = -45; thi2 = -45; qd = 1.5; qg = 1.5;
        thinf = -50; qinf = 4;
        mmin = 0.02; hmin = 0.5;
        sh = get_field_or_default(hh,'sh_na3',0);

        a_m = trap0_fun(V, tha+sh, Ra, qa);
        b_m = trap0_fun(-V, -tha-sh, Rb, qa);
        denom_m = a_m + b_m;
        g.m_inf = a_m ./ denom_m;
        g.tau_m_ms = 1 ./ denom_m ./ qt_na;
        g.tau_m_ms = max(g.tau_m_ms, mmin);

        a_h = trap0_fun(V, thi1+sh, 0.03, qd);
        b_h = trap0_fun(-V, -thi2-sh, 0.01, qg);
        denom_h = a_h + b_h;
        g.h_inf = 1 ./ (1 + exp((V - thinf - sh)./qinf));
        g.tau_h_ms = 1 ./ denom_h ./ qt_na;
        g.tau_h_ms = max(g.tau_h_ms, hmin);

        % slow inactivation s (na3n.mod)
        vhalfs = -60; zetas = 12; gms = 0.2; a0s = 0.0003; smin = 10;
        vvh = -58; vvs = 2; ar = 1; % default from na3n.mod (no slow inact reduction)
        alpv = 1 ./ (1 + exp((V - vvh - sh)./vvs));
        alps = exp(1e-3 .* zetas .* (V - vhalfs - sh) .* 9.648e4 ./ (8.315 .* (273.16 + hh.celsius)));
        bets = exp(1e-3 .* zetas .* gms .* (V - vhalfs - sh) .* 9.648e4 ./ (8.315 .* (273.16 + hh.celsius)));
        g.s_inf = alpv + ar .* (1 - alpv);
        g.tau_s_ms = bets ./ (a0s .* (1 + alps));
        g.tau_s_ms = max(g.tau_s_ms, smin);

        % ---- kdrbca1.mod ----
        sh_kdrb = get_field_or_default(hh,'sh_kdrb',0);
        [g.n_inf, g.tau_n_ms] = kd_family_rates(V, hh.celsius, 13, -3, 0.7, 0.02, 2, sh_kdrb);

        % ---- kdb.mod ----
        sh_kdb = get_field_or_default(hh,'sh_kdb',0);
        [g.n_kdb_inf, g.tau_n_kdb_ms] = kd_family_rates(V, hh.celsius, -33, 3, 0.7, 0.005, 2, sh_kdb);

        g.p_open_Na = g.m_inf.^3 .* g.h_inf .* g.s_inf;
        g.p_open_Kdrb = g.n_inf;
        g.p_open_Kdb = g.n_kdb_inf;
        g.gbar_Na = get_field_or_default(hh,'gbar_na3',nan);
        g.gbar_Kdrb = get_field_or_default(hh,'gbar_kdrb',nan);
        g.gbar_Kdb = get_field_or_default(hh,'gbar_kdb',nan);

        g.m_inf    = reshape(g.m_inf,    size(Vin));
        g.h_inf    = reshape(g.h_inf,    size(Vin));
        g.n_inf    = reshape(g.n_inf,    size(Vin));
        g.s_inf    = reshape(g.s_inf,    size(Vin));
        g.n_kdb_inf = reshape(g.n_kdb_inf, size(Vin));
        g.tau_m_ms = reshape(g.tau_m_ms, size(Vin));
        g.tau_h_ms = reshape(g.tau_h_ms, size(Vin));
        g.tau_n_ms = reshape(g.tau_n_ms, size(Vin));
        g.tau_s_ms = reshape(g.tau_s_ms, size(Vin));
        g.tau_n_kdb_ms = reshape(g.tau_n_kdb_ms, size(Vin));
        g.p_open_Na = reshape(g.p_open_Na, size(Vin));
        g.p_open_Kdrb = reshape(g.p_open_Kdrb, size(Vin));
        g.p_open_Kdb = reshape(g.p_open_Kdb, size(Vin));

    otherwise
        error('hh_compute_gating: unknown hh.kind "%s".', hh.kind);
end

end

% =========================================================
% Helpers (stable “vtrap” variants) for Canakci pyramidal
% =========================================================
function out = vtrap_exp(x, y)
% out = x / (exp(x/y) - 1), with Taylor fix near 0
x = x(:);
z = x ./ y;
out = x ./ (exp(z) - 1);
small = abs(z) < 1e-6;
out(small) = y .* (1 - z(small)./2);
end

function out = vtrap_1mexp(x, y)
% out = x / (1 - exp(-x/y)), with Taylor fix near 0
x = x(:);
z = x ./ y;
out = x ./ (1 - exp(-z));
small = abs(z) < 1e-6;
out(small) = y .* (1 + z(small)./2);
end

% =========================================================
% BWB-style helpers (fun1/fun2/fun3) for Canakci PV/basket
% =========================================================
function y = fun1(v, V0, A, B)
% A * exp((v - V0)/B)
y = A .* exp((v - V0)./B);
end

function y = fun2(v, V0, A, B)
% A / (exp((v - V0)/B) + 1)
y = A ./ (exp((v - V0)./B) + 1);
end

function y = fun3(v, V0, A, B)
% A*(v-V0)/(exp((v-V0)/B)-1) with stable expm1 handling
x = (v - V0)./B;
den = expm1(x);                 % exp(x)-1
y = A .* (v - V0) ./ den;

small = abs(x) < 1e-6;
% limit as (v-V0)->0: (v-V0)/(exp((v-V0)/B)-1) -> B
% so fun3 -> A*B
y(small) = A .* B .* (1 - x(small)./2);  % 1st-order series correction
end


function val = get_field_or_default(s, field_name, default_val)
if isfield(s, field_name)
    val = s.(field_name);
else
    val = default_val;
end
end

function y = trap0_fun(v, th, a, q)
x = v - th;
y = a .* x ./ (1 - exp(-x./q));
small = abs(x) <= 1e-6;
y(small) = a .* q;
end

function [ninf, taun] = kd_family_rates(V, celsius, vhalfn, zetan, gmn, a0n, nmin, sh)
if nargin < 8
    sh = 0;
end
q10 = 1;
qt = q10.^((celsius-24)./10);
alpn = exp(1e-3 .* zetan .* (V - vhalfn - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
betn = exp(1e-3 .* zetan .* gmn .* (V - vhalfn - sh) .* 9.648e4 ./ (8.315 .* (273.16 + celsius)));
ninf = 1 ./ (1 + alpn);
taun = betn ./ (qt .* a0n .* (1 + alpn));
taun = max(taun, nmin./qt);
end
