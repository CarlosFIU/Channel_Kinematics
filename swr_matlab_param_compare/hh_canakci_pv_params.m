function hh = hh_canakci_pv_params(celsius)
% Canakci et al. (2017) basket/PV kinetics (BWB mechanisms)
% Source mechanisms (ModelDB 230861):
%   - mod/nafbwb.mod  (fast Na; m is instantaneous; h is state)
%   - mod/kdrbwb.mod  (KDR; n is state)
%
% The mechanisms include Q10 scaling (phi=5) referenced to 27°C.
% We store temperature here; hh_compute_gating() uses it.

if nargin < 1 || isempty(celsius)
    celsius = 34;  % default commonly used in hippocampal modeling
end

hh.kind    = 'canakci_pv';
hh.celsius = celsius;
hh.name    = sprintf('Canakci/Fink PV basket');

end
