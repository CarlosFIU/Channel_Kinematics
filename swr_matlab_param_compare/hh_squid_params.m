function hh = hh_squid_params()
% hh_squid_params
% Classic Hodgkin-Huxley squid axon kinetics (1952 formulation).
% Returns a struct with alpha/beta function handles and conductance params.
%
% Units: V in mV, time in ms.
hh.name = 'Squid HH (1952)';
hh.kind = 'alphabeta';
hh.Cm_uFcm2  = 1.0;     % uF/cm^2
hh.gNa_mScm2 = 120.0;   % mS/cm^2
hh.gK_mScm2  = 36.0;    % mS/cm^2
hh.gL_mScm2  = 0.3;     % mS/cm^2
hh.ENa_mV    = 50.0;
hh.EK_mV     = -77.0;
hh.EL_mV     = -54.387;

% Alpha/Beta functions (avoid singularities using small eps)
epsv = 1e-9;

hh.alpha_m = @(V) (0.1*(V+40.0)) ./ (1 - exp(-(V+40.0)/10.0) + epsv);
hh.beta_m  = @(V) 4.0 * exp(-(V+65.0)/18.0);

hh.alpha_h = @(V) 0.07 * exp(-(V+65.0)/20.0);
hh.beta_h  = @(V) 1 ./ (1 + exp(-(V+35.0)/10.0));

hh.alpha_n = @(V) (0.01*(V+55.0)) ./ (1 - exp(-(V+55.0)/10.0) + epsv);
hh.beta_n  = @(V) 0.125 * exp(-(V+65.0)/80.0);
end
