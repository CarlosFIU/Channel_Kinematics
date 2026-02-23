MATLAB mini-package contents
---------------------------
1) swr_param_compare.m
   - Builds a small parameter table (synaptic tau rise/decay) from SWR-model papers
   - Plots rise/decay parameters and kernel shapes
   - Includes a scaffold for HH gating comparisons

2) hh_squid_params.m + hh_compute_gating.m + hh_simulate.m
   - Classic HH (1952 squid axon) example implementation to compare gating curves.
   - Add new models by defining their alpha/beta (or m_inf/tau) functions in a similar way.

How to extend to your SWR model papers
--------------------------------------
- If a paper provides only Exp2Syn taus: add tau_rise/tau_decay.
- If a paper uses a single exponential synaptic decay: set tau_rise = NaN, tau_decay = value.
- For channel kinetics (HH gating), you usually need the simulator code:
    * NEURON: .mod files (ion channel kinetics), .hoc templates (conductances), parameter JSONs
    * Brian/Matlab: code defines alpha/beta or m_inf/tau explicitly
  Once you have those, implement the gating functions in MATLAB and run the same plots.

3) hh_channel_kinematics_pc_params.m + hh_channel_kinematics_pvbc_params.m
   - Load PC/PVBC soma conductances from repository SONATA JSON parameter files.
   - Used by swr_param_compare.m to compare gate tau and channel opening probabilities
     using the repository mod-file equations (naxn/kdrca1 for PC; na3n/kdrbca1/kdb for PVBC).
