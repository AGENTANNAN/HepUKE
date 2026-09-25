# =====================================================================
# Dataset preparation
# =====================================================================
# psi(2S) data @ 3.686 GeV (106 M psi' events, 156.4 pb^-1) + inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Off-resonance data @ 3.65 GeV (44.1 pb^-1) + inclusive MC
off_data   = DatasetManager.real_data.find("709_3650")
off_incMC  = DatasetManager.inclusive_mc.find("709_3650")

# psi(3770) data @ 3.773 GeV (921.8 pb^-1) + inclusive MC
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# =====================================================================
# Decay cards (EvtGen format)
# =====================================================================
# Signal: psi(2S) -> gamma chi_c0, chi_c0 -> gamma gamma
decay_card_c0_gg = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 P2GC0;
  Enddecay

  Decay chi_c0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Signal: psi(2S) -> gamma chi_c2, chi_c2 -> gamma gamma
decay_card_c2_gg = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 P2GC2;
  Enddecay

  Decay chi_c2
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Control: psi(2S) -> gamma chi_c0, chi_c0 -> K+ K-
decay_card_c0_kk = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 P2GC0;
  Enddecay

  Decay chi_c0
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# Control: psi(2S) -> gamma chi_c2, chi_c2 -> K+ K-
decay_card_c2_kk = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 P2GC2;
  Enddecay

  Decay chi_c2
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# Peaking background: chi_c0 -> pi0 pi0
decay_card_c0_pi0pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 P2GC0;
  Enddecay

  Decay chi_c0
  1.000 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background: chi_c2 -> eta eta
decay_card_c2_etaeta = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 P2GC2;
  Enddecay

  Decay chi_c2
  1.000 eta eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Continuum (off-resonance) background e+ e- -> gamma gamma (gamma).
# Top mother psi(4260) is the KKMC convention for a QED continuum final state.
decay_card_babayaga = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# =====================================================================
# Exclusive MC samples (100k events each)
# =====================================================================
exMC_c0_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic0_ggg"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c0_gg
  config.cross_section   = :default
end

exMC_c2_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_ggg"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c2_gg
  config.cross_section   = :default
end

exMC_c0_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic0_kk"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c0_kk
  config.cross_section   = :default
end

exMC_c2_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_kk"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c2_kk
  config.cross_section   = :default
end

exMC_c0_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_chic0_pi0pi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c0_pi0pi0
  config.cross_section   = :default
end

exMC_c2_etaeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_chic2_etaeta"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_c2_etaeta
  config.cross_section   = :default
end

# Shared Babayaga e+e- -> gamma gamma (gamma) sample used for BOTH off-resonance
# points (3.65 GeV and 3.773 GeV) -> one signal MC per energy point.
exMC_babayaga = DatasetManager.create_exclusive_mc_for([off_data, psi3770_data]) do |config|
  config.sample_name   = "exmc_babayaga_gg"
  config.events        = 100000
  config.decay_card    = decay_card_babayaga
  config.cross_section = :default
end

# =====================================================================
# Event selection (BOSS)
# =====================================================================
alg_name = "ChiCGammaGG"   # psi(2S) -> gamma chi_c(J) -> gamma gamma gamma
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) c.m. energy
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common three-photon chain: chi_c0 and chi_c2 share the identical final state
# (gamma gamma gamma) and therefore one algorithm / one selection chain.
event_selection = Selection.new
event_selection.select_track {          # No charged tracks in gamma gamma gamma
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     "==0"        # Zero positive charged tracks
                  nChrn     "==0"        # Zero negative charged tracks
                  nNet      "==0"        # Net charge zero
                }
               .select_photon {          # Exactly three photon candidates
                  tdc_emc_start     0    # EMC timing window [0, 14] (x 700 ns)
                  tdc_emc_end       14
                  angle_to_track    10.0 # >=10 deg to the nearest charged track
                  energyThreshold_b 0.025 # E > 25 MeV (barrel)
                  energyThreshold_e 0.050 # E > 50 MeV (endcap)
                  nGam              "==3" # Exactly 3 photons
                }
               .kinematic_fit([:gamma, :gamma, :gamma]) {  # 4C fit over the three photons
                  nominal                    # Nominal fit (chi2 and 4-momenta from here)
                  constrain_four_momentum    # 4C energy-momentum conservation
                  chi2_cut 200               # Loose cut in BOSS; chi2 <= 80 applied in ROOT analysis
                }

# BOSS-side procedures that have no dedicated DSL primitive
my_algorithm
  .note(:photon_fiducial, "each photon is required to satisfy a fiducial |cos(theta)| < 0.75 " \
        "with respect to the e+ beam direction; applied after select_photon because the photon " \
        "selection block exposes no polar-angle cut")
  .note(:photon_assignment, "in each gamma gamma gamma candidate the lowest-energy photon is " \
        "assigned to the radiative gamma1 and the two higher-energy photons to the chi_c0/chi_c2 " \
        "decay, used to form the chi_c invariant mass")
  .note(:babayaga_generator, "the continuum e+e- -> gamma gamma (gamma) sample for the two " \
        "off-resonance points is produced with the Babayaga generator including ISR and is shared " \
        "between the 3.65 GeV and 3.773 GeV data sets")

# Generate the complete algorithm (kinematic variables from the chi_c0 -> gamma gamma card;
# the chi_c2 channel has the identical gamma gamma gamma final state)
my_algorithm.with_decay_card(decay_card_c0_gg).apply(event_selection)

# Execute on data, inclusive MC, signal MC, peaking-background MC and the
# off-resonance continuum MC. (The gamma K+ K- control samples use a dedicated
# K+ K- gamma selection and are reconstructed in a companion spec.)
root_files = my_algorithm.execute_on([
  psip_data, psip_incMC,
  off_data, off_incMC,
  psi3770_data, psi3770_incMC,
  exMC_c0_gg, exMC_c2_gg,
  exMC_c0_pi0pi0, exMC_c2_etaeta,
  *exMC_babayaga
])