function hh = hh_canakci_pyr_params()
% Canakci et al. (2017) / Fink et al. (2015) CA1 PYRAMIDAL kinetics
% Source mechanisms (ModelDB 230861):
%   - mod/ca1ina.mod  (Na: m, h, i)
%   - mod/ca1ikdr.mod (KDR: n)

hh.kind = 'canakci_fink';
hh.name = 'Canakci/Fink CA1 pyramidal';

% Parameters used in ca1ina.mod for slow inactivation i_inf:
hh.vi = -60;  % mV
hh.ki = 0.8;

end
