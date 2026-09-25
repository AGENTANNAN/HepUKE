# BESIII Analysis: First observation of hc -> e+e- eta_c (EM Dalitz)
# Paper: 2407.00136v2
# Mode I: psi(3686) -> pi0 hc,  Mode II: e+e- -> pi+pi- hc
# Ratio R = B(hc -> e+e- eta_c) / B(hc -> gamma eta_c); eta_c undetected

### Dataset preparation ###
# Mode I: psi(2S) data
data_psip = DatasetManager.real_data.find("709_3686")
incMC_psip = DatasetManager.inclusive_mc.find("709_3686")

# Mode II: XYZ scan data
data_4130 = DatasetManager.real_data.find("705_4130")    # 4.130 GeV
data_4160 = DatasetManager.real_data.find("705_4160")    # 4.160 GeV
data_4210 = DatasetManager.real_data.find("703_4210")    # 4.210 GeV
data_4230 = DatasetManager.real_data.find("703_4230")    # 4.230 GeV
data_4237 = DatasetManager.real_data.find("703_4237")    # 4.237 GeV
data_4246 = DatasetManager.real_data.find("703_4246")    # 4.246 GeV
data_4260 = DatasetManager.real_data.find("703_4260")    # 4.260 GeV
data_4290 = DatasetManager.real_data.find("705_4290")    # 4.290 GeV
data_4315 = DatasetManager.real_data.find("705_4315")    # 4.315 GeV
data_4340 = DatasetManager.real_data.find("705_4340")    # 4.340 GeV
data_4360 = DatasetManager.real_data.find("703_4360")    # 4.360 GeV
data_4380 = DatasetManager.real_data.find("705_4380")    # 4.380 GeV
data_4400 = DatasetManager.real_data.find("705_4400")    # 4.400 GeV
data_4420 = DatasetManager.real_data.find("703_4420")    # 4.420 GeV
data_4440 = DatasetManager.real_data.find("705_4440")    # 4.440 GeV
data_4750 = DatasetManager.real_data.find("707_4750")    # 4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")    # 4.780 GeV

scan_points = [
  data_4130, data_4160, data_4210, data_4230,
  data_4237, data_4246, data_4260, data_4290,
  data_4315, data_4340, data_4360, data_4380,
  data_4400, data_4420, data_4440, data_4750, data_4780
]

# ============================================================
# Mode I: psi(3686) -> pi0 hc
# ============================================================

# Decay card: Mode I signal — hc -> e+e- eta_c
decay_card_modeI_sig = <<~DECAYCARD
  Decay psi(2S)
  1.0 pi0 hc HELAMP 0 0 1 0 1 0 0 0;
  Enddecay

  Decay hc
  1.0 e+ e- eta_c PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: Mode I reference — hc -> gamma eta_c
decay_card_modeI_ref = <<~DECAYCARD
  Decay psi(2S)
  1.0 pi0 hc HELAMP 0 0 1 0 1 0 0 0;
  Enddecay

  Decay hc
  1.0 gamma eta_c E1RAD;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for Mode I signal
exMC_modeI_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc2ee_etac_modeI"
  config.related_dataset = data_psip
  config.events          = 500_000
  config.decay_card      = decay_card_modeI_sig
  config.cross_section   = :default
end

# Exclusive MC for Mode I reference
exMC_modeI_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc2gamma_etac_modeI"
  config.related_dataset = data_psip
  config.events          = 500_000
  config.decay_card      = decay_card_modeI_ref
  config.cross_section   = :default
end

# --- Mode I Signal Algorithm: psi(3686) -> pi0 hc, hc -> e+e- eta_c ---
alg_modeI_sig = Algorithm.new("HcEEetacModeI")
alg_modeI_sig.set_header(["HcEEetacModeIAlg/HcEEetacModeI.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })

sel_modeI_sig = Selection.new

# Charged tracks: e+e- pair
sel_modeI_sig.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
end

# Photon selection: at least 2 for pi0
sel_modeI_sig.select_photon do
  tdc_emc_start   0
  tdc_emc_end     700
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track  10.0
  nGam            ">=2"
end

# Electron identification for e+e- pair
sel_modeI_sig.pid(method: :probability) do
  prob_cut   0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  nep        "==1"
  nem        "==1"
end

# Reconstruct pi0 from gamma pairs
sel_modeI_sig.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0       ">=1"
end

# Partial reconstruction: reconstruct pi0 + e+e-, miss eta_c
# RecID mapping: 0=psi(2S)[skip], 1=pi0, 2=hc, 3=gamma1(pi0), 4=gamma2(pi0), 5=ep, 6=em, 7=eta_c[miss]
sel_modeI_sig.partial_miss([7]) do
  best_combination_by_mass :hc, 3.525    # constrain to hc PDG mass
  require_recoil_mass 2.92, 3.08        # eta_c mass region
end

alg_modeI_sig.with_decay_card(decay_card_modeI_sig).apply(sel_modeI_sig)

alg_modeI_sig
  .note(:ep_over_p_cut,
    "E/p cut 0.5–1.2 applied to higher-momentum e+/e- track for improved eID; " \
    "applied in BOSS-level selection on electron candidates")
  .note(:pi0_gamma_energy,
    "Photons from pi0 candidate must be in EMC barrel with E > 40 MeV; " \
    "applied in BOSS pi0 reconstruction")
  .note(:pi0_mass_window,
    "M(gamma gamma) in [120, 145] MeV/c2 for pi0 candidates; " \
    "applied before Kalman kinematic fit")
  .note(:jpsi_veto,
    "Recoil mass of each pi+pi- (gammagamma) pair vetoed in J/psi +/- 7 (30) MeV/c2; " \
    "suppresses psi(2S) -> pi+pi- J/psi and gamma gamma J/psi backgrounds")
  .note(:pi0_to_ee_gamma_veto,
    "M(e+e- gamma') vetoed in pi0 +/- 15 MeV/c2 to suppress pi0 -> gamma e+e-; " \
    "gamma' is either photon from pi0 candidate")
  .note(:photon_conversion_veto,
    "Photon conversion finder applied; delta_xy < 2 cm rejects conversions at beam pipe " \
    "and MDC inner wall")
  .note(:e1_gamma_pi0_veto,
    "For hc -> gamma eta_c reference: E1 gamma combined with any other photon " \
    "vetoed if M(gammagamma) in pi0 +/- 15 MeV/c2")
  .note(:ee_energy_window,
    "E(e+e-) in [470, 540] MeV for eta_c selection in hc -> e+e- eta_c; " \
    "E(gamma) in same range for hc -> gamma eta_c reference")

# --- Mode I Reference Algorithm: psi(3686) -> pi0 hc, hc -> gamma eta_c ---
alg_modeI_ref = Algorithm.new("HcGammaEtacModeI")
alg_modeI_ref.set_header(["HcGammaEtacModeIAlg/HcGammaEtacModeI.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })

sel_modeI_ref = Selection.new

# Photon selection: at least 3 (2 for pi0 + 1 E1 gamma)
sel_modeI_ref.select_photon do
  tdc_emc_start   0
  tdc_emc_end     700
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track  10.0
  nGam            ">=3"
end

# Reconstruct pi0 from gamma pairs
sel_modeI_ref.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0       ">=1"
end

# Partial reconstruction: reconstruct pi0 + gamma, miss eta_c
# RecID: 0=psi(2S)[skip], 1=pi0, 2=hc, 3=gamma1(pi0), 4=gamma2(pi0), 5=gamma(E1), 6=eta_c[miss]
sel_modeI_ref.partial_miss([6]) do
  best_combination_by_mass :hc, 3.525
  require_recoil_mass 2.92, 3.08
end

alg_modeI_ref.with_decay_card(decay_card_modeI_ref).apply(sel_modeI_ref)

alg_modeI_ref
  .note(:pi0_mass_window,
    "M(gamma gamma) in [120, 145] MeV/c2 for pi0 candidates")
  .note(:pi0_gamma_energy,
    "Photons from pi0 candidate must be in EMC barrel with E > 40 MeV")
  .note(:jpsi_veto,
    "Recoil mass veto for J/psi background suppression")
  .note(:e1_gamma_pi0_veto,
    "E1 gamma vetoed if M(gammagamma) in pi0 window with any other photon")
  .note(:e1_gamma_energy_window,
    "E(gamma) in [470, 540] MeV for eta_c selection")

# ============================================================
# Mode II: e+e- -> pi+pi- hc at XYZ scan points
# ============================================================

# Decay card: Mode II signal — hc -> e+e- eta_c
decay_card_modeII_sig = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- hc PHSP;
  Enddecay

  Decay hc
  1.0 e+ e- eta_c PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: Mode II reference — hc -> gamma eta_c
decay_card_modeII_ref = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- hc PHSP;
  Enddecay

  Decay hc
  1.0 gamma eta_c E1RAD;
  Enddecay

  End
DECAYCARD

# Exclusive MC for Mode II signal (scan)
exMC_modeII_sig = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_hc2ee_etac_modeII"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII_sig
  config.cross_section = :default
end

# Exclusive MC for Mode II reference (scan)
exMC_modeII_ref = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_hc2gamma_etac_modeII"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII_ref
  config.cross_section = :default
end

# --- Mode II Signal Algorithm: e+e- -> pi+pi- hc, hc -> e+e- eta_c ---
alg_modeII_sig = Algorithm.new("HcEEetacModeII")
alg_modeII_sig.set_header(["HcEEetacModeIIAlg/HcEEetacModeII.h"])

sel_modeII_sig = Selection.new

# Charged tracks: pi+pi- from hc production + e+e- from hc decay
sel_modeII_sig.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=2"
end

# Electron and pion identification
sel_modeII_sig.pid(method: :probability) do
  prob_cut   0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  nep        "==1"
  nem        "==1"
  npip       "==1"    # production pi+
  npim       "==1"    # production pi-
end

# Partial reconstruction: reconstruct pi+pi- + e+e-, miss eta_c
# RecID: 0=psi(4260)[skip], 1=pip, 2=pim, 3=hc, 4=ep, 5=em, 6=eta_c[miss]
sel_modeII_sig.partial_miss([6]) do
  best_combination_by_mass :hc, 3.525
  require_recoil_mass 2.92, 3.08
end

alg_modeII_sig.with_decay_card(decay_card_modeII_sig).apply(sel_modeII_sig)

alg_modeII_sig
  .note(:ep_over_p_cut,
    "E/p cut 0.5–1.2 applied to higher-momentum e+/e- track for improved eID")
  .note(:photon_conversion_veto,
    "Photon conversion finder: delta_xy < 2 cm")
  .note(:ee_energy_window,
    "E(e+e-) in [470, 540] MeV for eta_c selection")

# --- Mode II Reference Algorithm: e+e- -> pi+pi- hc, hc -> gamma eta_c ---
alg_modeII_ref = Algorithm.new("HcGammaEtacModeII")
alg_modeII_ref.set_header(["HcGammaEtacModeIIAlg/HcGammaEtacModeII.h"])

sel_modeII_ref = Selection.new

# Charged tracks: pi+pi- from hc production
sel_modeII_ref.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
end

# Photon: E1 gamma
sel_modeII_ref.select_photon do
  tdc_emc_start   0
  tdc_emc_end     700
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track  10.0
  nGam            ">=1"
end

# Pion identification for production pi+pi-
sel_modeII_ref.pid(method: :probability) do
  prob_cut   0.001
  identify :pion, against: [:kaon]
  npip       "==1"
  npim       "==1"
end

# Partial reconstruction: reconstruct pi+pi- + gamma, miss eta_c
# RecID: 0=psi(4260)[skip], 1=pip, 2=pim, 3=hc, 4=gamma, 5=eta_c[miss]
sel_modeII_ref.partial_miss([5]) do
  best_combination_by_mass :hc, 3.525
  require_recoil_mass 2.92, 3.08
end

alg_modeII_ref.with_decay_card(decay_card_modeII_ref).apply(sel_modeII_ref)

alg_modeII_ref
  .note(:e1_gamma_energy_window,
    "E(gamma) in [470, 540] MeV for eta_c selection")
  .note(:e1_gamma_pi0_veto,
    "E1 gamma combined with any other photon vetoed if M(gammagamma) in pi0 window")

# Execute Mode I algorithms on psi(2S) data
alg_modeI_sig.execute_on([data_psip, incMC_psip, exMC_modeI_sig])
alg_modeI_ref.execute_on([data_psip, incMC_psip, exMC_modeI_ref])

# Execute Mode II algorithms on scan points
all_scan_datasets = scan_points + [exMC_modeII_sig].flatten + [exMC_modeII_ref].flatten
alg_modeII_sig.execute_on(all_scan_datasets)
alg_modeII_ref.execute_on(all_scan_datasets)