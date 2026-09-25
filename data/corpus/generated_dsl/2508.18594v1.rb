# ============================================================================
# Search for the Λc–Σc bound state Hc in e+e- -> π+ Hc-, Hc- -> π- Λc+ Λ̄c-
# Partial reconstruction of the Λc+ side; the Λ̄c- is treated as missing.
# ============================================================================

### ------------------------------- Datasets -------------------------------
data_4914  = DatasetManager.real_data.find("707_4914")      # √s = 4918.02 MeV
data_4946  = DatasetManager.real_data.find("707_4946")      # √s = 4950.93 MeV
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")   # corresponding inclusive MC
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

### ------------------------------ Decay card ------------------------------
# e+e- -> π+π- Λc+ Λ̄c- (KKMC top mother ψ(4260)),
# Λc+ -> p K- π+ and the charge-conjugate decay (PHSP).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ K- pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0 anti-p- K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------- Exclusive MC samples -------------------------
# 100k-event exclusive MC at each energy point (same card / cross section,
# one sample per data point).
exMC_signal = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_hc_pipi_lambdac"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### ------------------------- Event selection (BOSS) -----------------------
alg_name = "HcSearch"
hc_alg   = Algorithm.new(alg_name)
hc_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 4.91802]})   # c.m. energy constant 4.91802 GeV

event_selection = Selection.new
  # charged tracks
  .select_track {
      cos_theta 0.93     # |cosθ| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr         1.0     # Vr < 1 cm
      nChrp     ">=2"    # at least two positive tracks
      nChrn     ">=2"    # at least two negative tracks
  }
  # photons
  .select_photon {
      tdc_emc_start     0      # EMC timing window 0 …
      tdc_emc_end      14      # … 14
      energyThreshold_b 0.025  # barrel E > 25 MeV
      energyThreshold_e 0.050  # endcap E > 50 MeV
      angle_to_track   10.0    # > 10° from the nearest charged track
  }
  # particle identification: probability method, 0.001 probability cut
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]    # p / p̄
      identify :kaon,   against: [:pion, :proton]  # K+ / K-
      identify :pion,   against: [:kaon, :proton]  # π+ / π-
      nprp ">=1"   # at least one p
      nkm  ">=1"   # at least one K-
      npip ">=1"   # at least one π+
  }
  # Λc+ candidate: secondary-vertex fit of p K- π+
  .secondary_vertex_fit([:prp, :km, :pip]) {
      build_virtual_particle(:Lambda_c).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Partial reconstruction of the Λc+ side; the Λ̄c- (recID 4, with its whole
  # decay subtree) is left unreconstructed and enters as the missing particle.
  .partial_miss([4]) {
      best_combination_by_mass :Lambda_c, 2.28646   # keep the Λc+ closest to its nominal mass
  }

# BOSS-side procedures that cannot be expressed in the DSL constructs above
hc_alg
  .note(:background_veto, "the Λc+ candidate is additionally required to lie in the
    invariant-mass window 2.270–2.300 GeV with fit χ² < 200 (nominal); the qq̄
    background is modelled from the Λc+ sidebands [2.190,2.250] and [2.320,2.380] GeV,
    together with the e+e- -> ΣcΣc, ΛcΣcπ, ΛcΛc(2595) and ΛcΛc(2625) channels, in the
    unbinned maximum-likelihood fit to the RM distributions.")
  .note(:remaining_pion_vertex_fit, "the remaining π+ and π- tracks are vertex-fitted
    together with the Λc+ (χ² < 200, combination with the minimum χ² kept); this step
    is not expressible inside the partial-reconstruction block.")
  .note(:efficiency_curve, "signal MC is generated for 15 Hc mass–width combinations
    (m_Hc = 4715–4735 MeV, Γ = 5/10/20 MeV) and used to model the RM signal shape; the
    90% C.L. Bayesian upper limit is obtained from the unbinned fit with multiplicative
    systematics via likelihood convolution and additive systematics taken as the most
    conservative limit across fit variations.")

# Generate the full algorithm for the process defined in the decay card
hc_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on both data points, their inclusive MC and the exclusive MC samples
root_files = hc_alg.execute_on([data_4914, data_4946, incMC_4914, incMC_4946] + exMC_signal)