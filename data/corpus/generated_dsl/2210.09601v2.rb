# ================================================================
# Dataset preparation
# ================================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi   (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi   inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(2S) inclusive MC

# ----------------------------------------------------------------
# Decay cards (EvtGen format)
# ----------------------------------------------------------------
# Signal at J/psi:  J/psi -> eta Sigma+ anti-Sigma-,
#   eta -> gamma gamma, Sigma+ -> p pi0, anti-Sigma- -> pbar pi0, pi0 -> gamma gamma
# (the negatively-charged Sigma in the reaction is the charge-conjugate anti-Sigma-,
#  whose pbar pi0 decay gives the required p / pbar final state)
decay_card_jpsi_signal = <<~DECAYCARD
  Decay J/psi
  1.000 eta Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Signal at psi(2S):  psi(2S) -> eta Sigma+ anti-Sigma- (same subsequent decays)
decay_card_psip_signal = <<~DECAYCARD
  Decay psi(2S)
  1.000 eta Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background at J/psi:  J/psi -> pi0 Sigma+ anti-Sigma-
decay_card_jpsi_bkg = <<~DECAYCARD
  Decay J/psi
  1.000 pi0 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background at psi(2S):  psi(2S) -> gamma chi_c0, chi_c0 -> pi0 Sigma+ anti-Sigma-
decay_card_psip_chi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi0 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background at psi(2S):  psi(2S) -> gamma chi_c1, chi_c1 -> pi0 Sigma+ anti-Sigma-
decay_card_psip_chi1 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.000 pi0 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background at psi(2S):  psi(2S) -> gamma chi_c2, chi_c2 -> pi0 Sigma+ anti-Sigma-
decay_card_psip_chi2 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.000 pi0 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ----------------------------------------------------------------
# Exclusive MC samples
# ----------------------------------------------------------------
# 500k-event phase-space signal MC at each resonance
exMC_jpsi_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_etaSigmaSigma"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_signal
  config.cross_section   = :default
end

exMC_psip_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_etaSigmaSigma"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_psip_signal
  config.cross_section   = :default
end

# Exclusive MC for the peaking backgrounds (used to model their shapes)
exMC_jpsi_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_pi0SigmaSigma"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_jpsi_bkg
  config.cross_section   = :default
end

exMC_psip_chi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachi0_pi0SigmaSigma"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_psip_chi0
  config.cross_section   = :default
end

exMC_psip_chi1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachi1_pi0SigmaSigma"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_psip_chi1
  config.cross_section   = :default
end

exMC_psip_chi2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachi2_pi0SigmaSigma"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_psip_chi2
  config.cross_section   = :default
end

# ================================================================
# Event selection (BOSS)
# ================================================================

# ----------------------------------------------------------------
# J/psi (3.097 GeV):  e+e- -> eta Sigma+ anti-Sigma-
# ----------------------------------------------------------------
alg_name_jpsi = "EtaSigmaSigmaJpsi"
alg_jpsi = Algorithm.new(alg_name_jpsi)
alg_jpsi.set_header(["#{alg_name_jpsi}Alg/#{alg_name_jpsi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = 3.097 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_jpsi = Selection.new
sel_jpsi
  .select_track {                    # Charged track selection
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        20.0                   # |Vz| < 20 cm
    nChrp     "==1"                  # exactly one positive track
    nChrn     "==1"                  # exactly one negative track
    nNet      "==0"                  # net charge zero
  }
  .select_photon {                   # Photon selection
    tdc_emc_start 0                  # EMC timing window [0, 14]
    tdc_emc_end   14
    energyThreshold_b 0.025          # 25 MeV in the barrel
    energyThreshold_e 0.050          # 50 MeV in the end-cap
    nGam ">=6"                       # at least six photons
  }
  .pid(method: :probability) {       # PID with the probability method
    prob_cut 0.001                   # probability > 0.001
    identify :proton, against: [:kaon, :pion]   # separate p / pbar from pi, K
    nprp "==1"                       # exactly one proton
    nprm "==1"                       # exactly one anti-proton
  }
  # Reconstruct two pi0 (gamma gamma, 1-C mass constraint) for the 7C fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # Reconstruct eta (gamma gamma, 1-C mass constraint) for the 7C fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # Nominal 7C kinematic fit: 4C energy-momentum conservation + eta + 2 pi0 mass constraints
  .kinematic_fit([:prp, :prm, :eta, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 30                     # chi^2 < 30 for J/psi
  }

# Background subtraction of the peaking J/psi -> pi0 Sigma+ anti-Sigma- component is
# performed on the NTuple with control samples and the dedicated exclusive-MC shapes.
alg_jpsi.note(:peaking_background_subtraction,
              "peaking J/psi -> pi0 Sigma+ anti-Sigma- background subtracted using control " \
              "samples and the shapes of the dedicated exclusive MC generated for this channel")

alg_jpsi.with_decay_card(decay_card_jpsi_signal).apply(sel_jpsi)

# ----------------------------------------------------------------
# psi(2S) (3.686 GeV):  e+e- -> eta Sigma+ anti-Sigma-
# ----------------------------------------------------------------
alg_name_psip = "EtaSigmaSigmaPsip"
alg_psip = Algorithm.new(alg_name_psip)
alg_psip.set_header(["#{alg_name_psip}Alg/#{alg_name_psip}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy = 3.686 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_psip = Selection.new
sel_psip
  .select_track {                    # Charged track selection
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        20.0                   # |Vz| < 20 cm
    nChrp     "==1"                  # exactly one positive track
    nChrn     "==1"                  # exactly one negative track
    nNet      "==0"                  # net charge zero
  }
  .select_photon {                   # Photon selection
    tdc_emc_start 0                  # EMC timing window [0, 14]
    tdc_emc_end   14
    energyThreshold_b 0.025          # 25 MeV in the barrel
    energyThreshold_e 0.050          # 50 MeV in the end-cap
    nGam ">=6"                       # at least six photons
  }
  .pid(method: :probability) {       # PID with the probability method
    prob_cut 0.001                   # probability > 0.001
    identify :proton, against: [:kaon, :pion]   # separate p / pbar from pi, K
    nprp "==1"                       # exactly one proton
    nprm "==1"                       # exactly one anti-proton
  }
  # Reconstruct two pi0 (gamma gamma, 1-C mass constraint) for the 7C fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # Reconstruct eta (gamma gamma, 1-C mass constraint) for the 7C fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # Nominal 7C kinematic fit: 4C energy-momentum conservation + eta + 2 pi0 mass constraints
  .kinematic_fit([:prp, :prm, :eta, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 25                     # chi^2 < 25 for psi(2S)
  }

# Background subtraction of the peaking psi(2S) -> gamma chi_c0,1,2 -> gamma pi0 Sigma+ anti-Sigma-
# components is performed on the NTuple with control samples and the exclusive-MC shapes.
alg_psip.note(:peaking_background_subtraction,
              "peaking psi(2S) -> gamma chi_c0,1,2, chi_c0,1,2 -> pi0 Sigma+ anti-Sigma- " \
              "backgrounds subtracted using control samples and the shapes of the dedicated " \
              "exclusive MC generated for these channels")

alg_psip.with_decay_card(decay_card_psip_signal).apply(sel_psip)

# ================================================================
# Execute on real data, inclusive MC and exclusive MC
# ================================================================
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC,
                                       exMC_jpsi_signal, exMC_jpsi_bkg])

root_files_psip = alg_psip.execute_on([psip_data, psip_incMC,
                                       exMC_psip_signal,
                                       exMC_psip_chi0, exMC_psip_chi1, exMC_psip_chi2])