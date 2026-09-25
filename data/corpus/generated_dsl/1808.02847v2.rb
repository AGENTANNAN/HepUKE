# frozen_string_literal: true

### Dataset preparation — e+e- -> pi+ D0 D*- ###
# Five high-luminosity points between 4.05 and 4.60 GeV (BOSS 703 release)
data_points = [
  DatasetManager.real_data.find("703_4230"),  # 4.2263 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.2580 GeV
  DatasetManager.real_data.find("703_4360"),  # 4.3583 GeV
  DatasetManager.real_data.find("703_4420"),  # 4.4156 GeV
  DatasetManager.real_data.find("703_4600")   # 4.5995 GeV
]

# Matched inclusive MC samples at the same energies
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600")
]

# Decay card: pi+ D0 D*- produced in phase space.
# D0 -> K- pi+ is fully reconstructed; D*- -> anti-D0 pi- is NOT reconstructed
# (it is inferred from the recoil mass against the D0 pi+ system).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ D0 D*-
    PHSP;
    Enddecay

    Decay D0
    1.0 K- pi+
    PHSP;
    Enddecay

    Decay D*-
    1.0 anti-D0 pi-
    PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi-
    PHSP;
    Enddecay

    End
DECAYCARD

# Signal MC: 500k events of the pi+ D0 D*- phase-space mode at each energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_piD0Dst_phsp"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PiD0Dst"
pi_d0_dst = Algorithm.new(alg_name)
pi_d0_dst.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.258]})

# Full selection chain up to (and substituting for) the kinematic fit.
event_selection = Selection.new
  .select_track {                                   # charged-track selection
    cos_theta 0.93                                  # |cos(theta)| < 0.93
    Vz        10.0                                  # |Vz| < 10 cm
    Vr        1.0                                   # Vr < 1 cm
    nTot      ">=4"                                 # at least 4 charged tracks
  }
  .pid(method: :probability) {                      # PID by the probability method
    prob_cut 0.001                                  # probability cut = 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion, :proton]       # K- (and K+) vs pi, p, and leptons
    identify :pion, against: [:kaon, :proton]       # pi+ (and pi-) vs K, p, and leptons
    nkm  ">=1"                                      # at least one K-
    npip ">=2"                                      # at least two pi+
  }
  # Partial reconstruction: tag the D0 (built from K- pi+) and the bachelor pi+;
  # the D*- is inferred from the recoil mass against the D0 pi+ system.
  # No kinematic fit is applied.
  .partial_rec([1, 2]) {
    best_combination_by_mass :D0, 1.8648            # D0 combination closest to m(D0)
    require_recoil_mass 1.85, 2.15                  # D*- region (loose)
  }

# BOSS-side procedures that have no formal DSL construct
pi_d0_dst
  .note(:d0_mass_window, "the D0 candidate is formed from the K- pi+ pair whose invariant mass is closest to the nominal m(D0); only K- pi+ combinations satisfying |M(K- pi+) - m(D0)| < 15 MeV/c^2 are kept")
  .note(:rm_cor_correction, "the corrected recoil mass RM_cor(D0 pi+) = RM(D0 pi+) + M(K- pi+) - m(D0) defines the D*- signal variable; the tight signal window |RM_cor - DeltaM - m(D*-)| < 20 MeV/c^2 and the sideband 1.91 < RM_cor < 1.95 GeV/c^2 are applied in ROOT, where the signal yield is extracted from an unbinned maximum-likelihood fit to RM_cor")
  .note(:background_veto, "e+e- -> D*D* background suppressed by rejecting events with M(D0 pi+) < 2.03 GeV/c^2")

pi_d0_dst.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, matched inclusive MC, and signal MC at all five points
root_files = pi_d0_dst.execute_on(data_points + incMC_points + exMCs_signal)