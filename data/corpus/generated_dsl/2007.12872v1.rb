# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# e+e- -> mu+mu- Born cross-section scan, sqrt(s) ~ 3.80 - 4.60 GeV.
# The DSL enumerates 45 representative BOSS 703/705/706/707 energy points.
data_names = %w[
  703_3810 703_3872 703_3900 703_4009 703_4090
  703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4245 703_4246 703_4260
  703_4270 703_4280 703_4310 703_4360 703_4390
  703_4420 703_4470 703_4530 703_4575 703_4600
  705_4130 705_4160 705_4290 705_4315 705_4340
  705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
data_points = data_names.map { |n| DatasetManager.real_data.find(n) }

# Inclusive MC is available for the subset of those points listed in the BESIII MC table.
incmc_names = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4237 703_4246 703_4260
  703_4270 703_4280 703_4360 703_4420 703_4600
  705_4130 705_4160
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
incmc_points = incmc_names.map { |n| DatasetManager.inclusive_mc.find(n) }

# Decay card for the signal process (EvtGen format).
# Continuum e+e- -> mu+mu- approximated with KKMC beam radiation + VLL (paper uses BABAYAGA).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 mu+ mu- VLL;
  Enddecay
  End
DECAYCARD

# Generate one 200k-event exclusive signal MC per scan energy point
# (same decay card / cross-section, differing only by the related real dataset).
exMCs_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_mumu_scan"   # auto-suffixed per energy point, e.g. exmc_mumu_scan_703_4180
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "MuMuScan"
mumu_alg = Algorithm.new(alg_name)
mumu_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 4.260]})            # nominal CMS energy for the 4C fit (per-run energy from conditions DB)
        .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain.
event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta 0.8                # |cos(theta)| < 0.8
                  Vz        10.0               # |Vz| < 10 cm
                  Vr        1.0                # Vr < 1 cm in the transverse plane
                  nChrp     "==1"              # exactly one positively charged track
                  nChrn     "==1"              # exactly one negatively charged track
                  nNet      "==0"              # net charge zero
                }
               .pid(method: :probability) {    # Muon identification (probability method)
                  prob_cut 0.001               # PID probability > 0.001
                  identify :muon, against: [:electron, :pion, :kaon, :proton]  # mu+ and mu- vs e, pi, K, p hypotheses
                  nmup "==1"                   # exactly one mu+
                  nmum "==1"                   # exactly one mu-
                }
               .kinematic_fit([:mup, :mum]) {  # 4C kinematic fit to the mu+mu- hypothesis
                  nominal
                  constrain_four_momentum      # constrain total four-momentum to the CMS energy
                  chi2_cut 60                  # chi^2 < 60
                }

# Capture BOSS-side procedures that have no DSL representation.
mumu_alg
  .note(:pid_correction_method, "paper-level tighter muon-ID E/p window 0.05 < E/p < 0.40 is not applied in the BOSS selection; deferred to ROOT-level selection")
  .note(:efficiency_curve, "energy-dependent total-momentum cut sum(|p(mu+)|,|p(mu-)|) > 0.9*sqrt(s) cannot be expressed with a single ECMS constant in BOSS; deferred to per-energy ROOT-level selection")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on all real data, inclusive MC and the per-point exclusive signal MC.
root_files = mumu_alg.execute_on(data_points + incmc_points + exMCs_mumu)