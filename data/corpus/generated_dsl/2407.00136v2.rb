# =============================================================================
#  h_c -> e+e- eta_c (EM Dalitz) / h_c -> gamma eta_c  with eta_c UNDETECTED
#  R = B(h_c -> e+e- eta_c) / B(h_c -> gamma eta_c)
#
#  Mode I  : psi(2S) real data + inclusive MC, 500k exclusive MC each for
#            psi(2S) -> pi0 h_c, h_c -> e+e- eta_c  and  h_c -> gamma eta_c
#  Mode II : 17 XYZ scan points (4.130-4.780 GeV) real data, 100k exclusive MC
#            per point for e+e- -> pi+pi- h_c, h_c -> e+e- eta_c / gamma eta_c
#
#  The eta_c is not reconstructed -> the h_c is built with the eta_c missing by
#  partial reconstruction (best mass combination at 3.525 GeV + recoil-mass
#  window 2.92-3.08 GeV on the eta_c).  partial_rec replaces the 4C kinematic
#  fit entirely.
# =============================================================================

### Datasets ###
# Mode I: psi(2S) data + inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Mode II: 17 XYZ scan points spanning 4.130 - 4.780 GeV (BOSS 705/706/707)
xyz_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780")
]

### Decay cards (EvtGen syntax) ###
# Mode I signal: psi(2S) -> pi0 h_c, h_c -> e+e- eta_c, pi0 -> gamma gamma
decay_card_modeI_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 e+ e- eta_c PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode I reference (normalisation): psi(2S) -> pi0 h_c, h_c -> gamma eta_c
decay_card_modeI_reference = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II signal: e+e- -> pi+pi- h_c (top mother psi(4260)), h_c -> e+e- eta_c
decay_card_modeII_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 e+ e- eta_c PHSP;
  Enddecay

  End
DECAYCARD

# Mode II reference (normalisation): e+e- -> pi+pi- h_c, h_c -> gamma eta_c
decay_card_modeII_reference = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
# Mode I: 500k events each, associated to the psi(2S) real data sample
exMC_modeI_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "modeI_psip_pi0hc_hc_ee_etac"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI_signal
  config.cross_section   = :default
end

exMC_modeI_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "modeI_psip_pi0hc_hc_gamma_etac"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI_reference
  config.cross_section   = :default
end

# Mode II: 100k events per energy point, one signal and one reference MC per point
exMCs_modeII_sig = DatasetManager.create_exclusive_mc_for(xyz_points) do |config|
  config.sample_name   = "modeII_pipihc_hc_ee_etac"
  config.events        = 100000
  config.decay_card    = decay_card_modeII_signal
  config.cross_section = :default
end

exMCs_modeII_ref = DatasetManager.create_exclusive_mc_for(xyz_points) do |config|
  config.sample_name   = "modeII_pipihc_hc_gamma_etac"
  config.events        = 100000
  config.decay_card    = decay_card_modeII_reference
  config.cross_section = :default
end

### Event selection (BOSS) ###

# -----------------------------------------------------------------------------
# Mode I signal:  psi(2S) -> pi0 h_c ,  h_c -> e+e- eta_c
#   decay card recIDs: 0=psi(2S) 1=pi0 2=h_c 3=gamma 4=gamma 5=e+ 6=e- 7=eta_c
# -----------------------------------------------------------------------------
alg_modeI_sig = Algorithm.new("HcEeEtacModeI")
alg_modeI_sig.set_header(["HcEeEtacModeIAlg/HcEeEtacModeI.h"])
             .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI_sig = Selection.new
  .select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nChrp     ">=1"       # at least 1 positive track
      nChrn     ">=1"       # at least 1 negative track
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025   # 25 MeV in the barrel
      energyThreshold_e 0.050   # 50 MeV in the endcap
      nGam              ">=2"   # at least 2 photons (pi0 -> gamma gamma)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> electron
      nlp "==1"                 # 1 e+
      nlm "==1"                 # 1 e-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # pi0 mass constraint (1C)
      chi2_cut 25
      npi0 ">=1"
  }
  # partial reconstruction: h_c built with the eta_c missing
  .partial_rec([1, 2, 3, 4, 5, 6]) {   # pi0, h_c, 2 gamma, e+, e-  (eta_c left out)
      best_combination_by_mass :h_c, 3.525   # best h_c mass combination
      require_recoil_mass 2.92, 3.08         # recoil (eta_c) mass window
  }

alg_modeI_sig
  .note(:ep_ratio_cut, "E/p of the electron candidates required in [0.5, 1.2]; applied at cut level on the reconstructed electron tracks")
  .note(:pi0_mass_window, "pi0 candidates required to lie in the mass window [0.120, 0.145] GeV/c2")
  .note(:pi0_photon_energy, "each photon from the pi0 required to have energy > 40 MeV")
  .note(:background_veto, "J/psi and pi0 -> e+e- gamma (Dalitz) backgrounds vetoed")
  .note(:photon_conversion_veto, "photon-conversion background suppressed by requiring Delta_xy > 2 cm")
  .note(:ee_energy_window, "E(e+e-) required in [0.470, 0.540] GeV")
  .with_decay_card(decay_card_modeI_signal)
  .apply(sel_modeI_sig)

# -----------------------------------------------------------------------------
# Mode I reference:  psi(2S) -> pi0 h_c ,  h_c -> gamma eta_c
#   decay card recIDs: 0=psi(2S) 1=pi0 2=h_c 3=gamma 4=gamma 5=gamma 6=eta_c
# -----------------------------------------------------------------------------
alg_modeI_ref = Algorithm.new("HcGammaEtacModeI")
alg_modeI_ref.set_header(["HcGammaEtacModeIAlg/HcGammaEtacModeI.h"])
             .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI_ref = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=3"   # pi0 -> gamma gamma plus the h_c -> gamma
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .partial_rec([1, 2, 3, 4, 5]) {   # pi0, h_c, 3 gamma  (eta_c left out)
      best_combination_by_mass :h_c, 3.525
      require_recoil_mass 2.92, 3.08
  }

alg_modeI_ref
  .note(:pi0_mass_window, "pi0 candidates required to lie in the mass window [0.120, 0.145] GeV/c2")
  .note(:pi0_photon_energy, "each photon from the pi0 required to have energy > 40 MeV")
  .note(:background_veto, "J/psi and pi0 -> e+e- gamma (Dalitz) backgrounds vetoed")
  .note(:photon_conversion_veto, "photon-conversion background suppressed by requiring Delta_xy > 2 cm")
  .with_decay_card(decay_card_modeI_reference)
  .apply(sel_modeI_ref)

# -----------------------------------------------------------------------------
# Mode II signal:  e+e- -> pi+pi- h_c ,  h_c -> e+e- eta_c
#   17-point XYZ scan (4.130-4.780 GeV)
#   decay card recIDs: 0=psi(4260) 1=pi+ 2=pi- 3=h_c 4=e+ 5=e- 6=eta_c
# -----------------------------------------------------------------------------
alg_modeII_sig = Algorithm.new("HcEeEtacModeII")
alg_modeII_sig.set_header(["HcEeEtacModeIIAlg/HcEeEtacModeII.h"])
              .set_constant({"ECMS" => [:double, 4.26]})   # placeholder, scan energy set per point
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII_sig = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"       # at least 2 positive tracks (pi+ e+)
      nChrn     ">=2"       # at least 2 negative tracks (pi- e-)
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      nlp   "==1"           # 1 e+
      nlm   "==1"           # 1 e-
      npip  "==1"           # 1 pi+
      npim  "==1"           # 1 pi-
  }
  .partial_rec([1, 2, 3, 4, 5]) {   # pi+, pi-, h_c, e+, e-  (eta_c left out)
      best_combination_by_mass :h_c, 3.525
      require_recoil_mass 2.92, 3.08
  }

alg_modeII_sig
  .note(:ep_ratio_cut, "E/p of the electron candidates required in [0.5, 1.2]; applied at cut level on the reconstructed electron tracks")
  .note(:background_veto, "J/psi and pi0 -> e+e- gamma (Dalitz) backgrounds vetoed")
  .note(:photon_conversion_veto, "photon-conversion background suppressed by requiring Delta_xy > 2 cm")
  .note(:ee_energy_window, "E(e+e-) required in [0.470, 0.540] GeV")
  .note(:ecms_scan, "Mode II is a 17-point XYZ scan (4.130-4.780 GeV); the ECMS constant above is a placeholder and the CMS energy is set per scan point at generation time")
  .with_decay_card(decay_card_modeII_signal)
  .apply(sel_modeII_sig)

# -----------------------------------------------------------------------------
# Mode II reference:  e+e- -> pi+pi- h_c ,  h_c -> gamma eta_c
#   decay card recIDs: 0=psi(4260) 1=pi+ 2=pi- 3=h_c 4=gamma 5=eta_c
# -----------------------------------------------------------------------------
alg_modeII_ref = Algorithm.new("HcGammaEtacModeII")
alg_modeII_ref.set_header(["HcGammaEtacModeIIAlg/HcGammaEtacModeII.h"])
              .set_constant({"ECMS" => [:double, 4.26]})   # placeholder, scan energy set per point
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII_ref = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=1"       # at least 1 positive track
      nChrn     ">=1"       # at least 1 negative track
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=1"   # at least 1 photon (h_c -> gamma)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip  "==1"           # 1 pi+
      npim  "==1"           # 1 pi-
  }
  .partial_rec([1, 2, 3, 4]) {   # pi+, pi-, h_c, gamma  (eta_c left out)
      best_combination_by_mass :h_c, 3.525
      require_recoil_mass 2.92, 3.08
  }

alg_modeII_ref
  .note(:background_veto, "J/psi and pi0 -> e+e- gamma (Dalitz) backgrounds vetoed")
  .note(:photon_conversion_veto, "photon-conversion background suppressed by requiring Delta_xy > 2 cm")
  .note(:ecms_scan, "Mode II is a 17-point XYZ scan (4.130-4.780 GeV); the ECMS constant above is a placeholder and the CMS energy is set per scan point at generation time")
  .with_decay_card(decay_card_modeII_reference)
  .apply(sel_modeII_ref)

### Execute on datasets ###
# Mode I: psi(2S) data + inclusive MC + the corresponding signal/reference exclusive MC
root_files_modeI_sig = alg_modeI_sig.execute_on([psip_data, psip_incMC, exMC_modeI_sig])
root_files_modeI_ref = alg_modeI_ref.execute_on([psip_data, psip_incMC, exMC_modeI_ref])

# Mode II: the 17 scan points plus the 100k/point signal/reference exclusive MC sets
root_files_modeII_sig = alg_modeII_sig.execute_on(xyz_points + exMCs_modeII_sig)
root_files_modeII_ref = alg_modeII_ref.execute_on(xyz_points + exMCs_modeII_ref)