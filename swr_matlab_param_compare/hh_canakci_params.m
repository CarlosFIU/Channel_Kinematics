function hh = hh_canakci_params()
% HH gating parameterization used in Canakci et al. (2017) model code,
% derived from Fink et al. (2015) eNeuro model.
%
% Source mechanisms:
%   - mod/ca1ina.mod  (Na: m, h, i)
%   - mod/ca1ikdr.mod (KDR: n)

hh.kind = 'canakci_fink';
hh.name = 'CA1 (Canakci/Fink)';

% Parameters used in ca1ina.mod for slow inactivation i_inf:
hh.vi = -60;  % mV
hh.ki = 0.8;

end
