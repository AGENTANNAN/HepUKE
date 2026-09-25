### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")    # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686") # corresponding inclusive MC

# --- Decay cards: three chi_cJ modes, all with the common final state gamma + 4 pi0 ---
decay_card_chi_c0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 P2GC0;
  Enddecay

  Decay chi_c0
  1.000 pi0 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chi_c1 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 pi0 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chi_c2 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 P2GC2;
  Enddecay

  Decay chi_c2
  1.000 pi0 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC: 500k events per chi_cJ mode ---
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chi_c0_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chi_c1_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chi_c2_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaChiC4Pi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})   # ECMS = 3.686 GeV

# The three chi_cJ modes share an identical final state (gamma + 4 pi0) and an
# identical selection chain, so a single Algorithm instance serves all of them.
event_selection = Selection.new
  .select_track {                          # charged-track requirement: zero tracks + quality cuts
    cos_theta 0.93                         # |cos(theta)| < 0.93
    Vz        10.0                         # |Vz| < 10 cm
    Vr        1.0                          # Vr < 1 cm
    nChrp     "==0"                        # no positive tracks
    nChrn     "==0"                        # no negative tracks
  }
  .select_photon {                         # exactly nine good photons
    tdc_emc_start     0                    # EMC timing window 0-14
    tdc_emc_end       14
    energyThreshold_b 0.025                # > 25 MeV in the barrel
    energyThreshold_e 0.050                # > 50 MeV in the endcap
    nGam              "==9"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                # pi0 reconstruction from photon pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # mass-constrain to pi0
    invariant_mass_of(:gamma, :gamma).within(0.110, 0.150)                # pi0 mass window 110-150 MeV/c^2
    chi2_cut 25                                                          # chi2 < 25 for the pi0 fit
    npi0     ">=4"                                                       # at least four pi0 candidates
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pi0]) {       # nominal 4C fit of gamma + 4 pi0
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# Inexpressible BOSS-side selection criteria preserved as notes
alg.note(:chi2_4pi_selection, "best combination of one radiative gamma plus four pi0 " \
        "chosen by minimizing chi2_4pi = sum_i (m_gg,i - m_pi0)^2 / sigma^2 with " \
        "sigma = 6.5 MeV/c^2; require chi2_4pi < 15")
   .note(:total_energy_window, "total energy of the nine photons required to lie in " \
        "[3.45, 3.80] GeV to suppress background")
   .note(:background_veto, "J/psi veto: reject events with any di-pi0 recoil mass within " \
        "100 MeV/c^2 of m_J/psi")
   .note(:background_veto, "K_S K_S veto for chi_c0 / chi_c2: reject events with " \
        "sqrt((m12 - m_K_S)^2 + (m34 - m_K_S)^2) < 100 MeV/c^2 for any assignment of the " \
        "four pi0 into two di-pi0 pairs")
   .with_decay_card(decay_card_chi_c0)
   .apply(event_selection)

# Execute on real data, inclusive MC and the three signal exclusive MC samples
root_files = alg.execute_on([psip_data, psip_incMC, exMC_chi_c0, exMC_chi_c1, exMC_chi_c2])