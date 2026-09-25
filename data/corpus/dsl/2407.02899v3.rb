# BESIII Analysis: J/psi -> p pbar eta branching fraction measurement
# Paper: 2407.02899v3
# Two eta channels: eta -> gamma gamma (I) and eta -> pi+ pi- pi0 (II)
# Both use 4C kinematic fit; eta mass unconstrained in both

### Dataset preparation ###
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Channel I: J/psi -> p pbar eta, eta -> gamma gamma
# ============================================================

decay_card_etagg = <<~DECAYCARD
  Decay J/psi
  1.0 p+ anti-p- eta PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_etagg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ppbar_etagg"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_card_etagg
  config.cross_section   = :default
end

alg_etagg = Algorithm.new("JpsiPpbarEtaGG")
alg_etagg.set_header(["JpsiPpbarEtaGGAlg/JpsiPpbarEtaGG.h"])
         .set_constant({ "ECMS" => [:double, 3.097] })

sel_etagg = Selection.new

# Charged tracks: p and anti-p (exactly 2 tracks)
sel_etagg.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
end

# Photons: at least 2 for eta -> gamma gamma
sel_etagg.select_photon do
  tdc_emc_start   0
  tdc_emc_end     700
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track  20.0
  nGam            ">=2"
end

# Proton identification — mass hypothesis only; no PID discrimination in real analysis
sel_etagg.pid(method: :probability) do
  prob_cut   0.001
  identify :proton, against: [:kaon]
  nprp       "==1"
  nprm       "==1"
end

# Reconstruct eta from gamma-gamma pairs with wide mass window
sel_etagg.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).between(0.200, 0.900)
  chi2_cut 200
  neta       ">=1"
end

# 4C kinematic fit — eta mass left unconstrained
sel_etagg.kinematic_fit([:prp, :prm, :eta]) do
  constrain_four_momentum
  chi2_cut 200
  nominal
end

alg_etagg.with_decay_card(decay_card_etagg).apply(sel_etagg)

alg_etagg
  .note(:no_pid_discrimination,
    "In the eta->gamma gamma channel, no PID discrimination is applied to charged tracks; " \
    "background suppression relies on the kinematic fit. Proton mass hypothesis is set " \
    "for kinematic fit purposes only.")
  .note(:photon_endcap_veto,
    "Combinations where both photons are in the EMC endcaps (|cos(theta)| > 0.8 / " \
    "0.86 < |cos(theta)| < 0.92) are rejected; applied in BOSS photon pairing")
  .note(:photon_angle_cut,
    "Photon candidate angle to nearest charged track > 20 deg " \
    "to reject bremsstrahlung and hadronic split-offs")
  .note(:vertex_fit,
    "Vertex fit applied to all charged tracks before kinematic fit to ensure common origin")
  .note(:best_candidate,
    "When multiple candidates per event, the one with minimum chi2 of the kinematic fit is selected")
  .note(:amplitude_model_mc,
    "Signal MC uses amplitude model (ComPWA helicity formalism, 7 N* resonances) " \
    "rather than pure PHSP; reweighting applied at ROOT analysis level")
  .note(:qed_background,
    "Continuum data at 3.080 GeV analyzed with same selection (adjusted CMS energy) " \
    "to estimate QED background; applied at ROOT level")
  .note(:eta_signal_region,
    "Eta signal region: M(gamma gamma) in [492, 587] MeV/c^2; " \
    "sideband regions [350, 462] and [632, 700] MeV/c^2; ROOT-level analysis")

# ============================================================
# Channel II: J/psi -> p pbar eta, eta -> pi+ pi- pi0
# ============================================================

decay_card_eta3pi = <<~DECAYCARD
  Decay J/psi
  1.0 p+ anti-p- eta PHSP;
  Enddecay

  Decay eta
  1.0 pi+ pi- pi0 ETA_DALITZ;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_eta3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ppbar_eta3pi"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta3pi
  config.cross_section   = :default
end

alg_eta3pi = Algorithm.new("JpsiPpbarEta3Pi")
alg_eta3pi.set_header(["JpsiPpbarEta3PiAlg/JpsiPpbarEta3Pi.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })

sel_eta3pi = Selection.new

# Charged tracks: p, anti-p, pi+, pi- (at least 4 total)
sel_eta3pi.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=2"
end

# Photons: at least 2 for pi0 -> gamma gamma
sel_eta3pi.select_photon do
  tdc_emc_start   0
  tdc_emc_end     700
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track  20.0
  nGam            ">=2"
end

# PID: proton vs pion discrimination; no kaon requirement
sel_eta3pi.pid(method: :probability) do
  prob_cut   0.001
  identify :proton, against: [:pion]
  identify :pion,   against: [:proton]
  nprp       "==1"
  nprm       "==1"
  npip       ">=1"
  npim       ">=1"
end

# Reconstruct pi0 from gamma-gamma pairs with mass constraint
sel_eta3pi.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0       ">=1"
end

# 4C + pi0 mass constraint kinematic fit — eta mass unconstrained
sel_eta3pi.kinematic_fit([:prp, :prm, :pip, :pim, :pi0]) do
  constrain_four_momentum
  invariant_mass_of(:pi0).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  nominal
end

alg_eta3pi.with_decay_card(decay_card_eta3pi).apply(sel_eta3pi)

alg_eta3pi
  .note(:eta_3pi_mass_window,
    "M(pi+ pi- pi0) in [200, 900] MeV/c^2 for eta candidate selection; " \
    "applied in BOSS-level selection before kinematic fit")
  .note(:photon_endcap_veto,
    "Combinations where both photons are in the EMC endcaps are rejected")
  .note(:photon_angle_cut,
    "Photon candidate angle to nearest charged track > 20 deg")
  .note(:vertex_fit,
    "Vertex fit applied to all charged tracks before kinematic fit")
  .note(:best_candidate,
    "When multiple candidates per event, the one with minimum chi2 of the kinematic fit is selected")
  .note(:amplitude_model_mc,
    "Signal MC uses amplitude model (ComPWA, 7 N* resonances); " \
    "reweighting applied at ROOT analysis level")
  .note(:qed_background,
    "Continuum data at 3.080 GeV analyzed with same selection to estimate QED background")
  .note(:eta_signal_region,
    "Eta signal region: M(pi+ pi- pi0) in [502, 602] MeV/c^2; " \
    "sideband regions [407, 492] and [622, 725] MeV/c^2; ROOT-level analysis")
  .note(:omega_veto,
    "M(pi+ pi- pi0) vetoed in omega(782) mass region to suppress " \
    "J/psi -> p pbar omega background")
  .note(:eta_pi_pi_gamma_peaking,
    "Peaking background from J/psi -> p pbar eta(eta -> pi+ pi- gamma) estimated " \
    "from inclusive MC and subtracted at ROOT level (~1.5% of signal)")
  .note(:pid_no_kaon,
    "No kaon likelihood requirement; protons: L(p) > L(pi); pions: L(pi) > L(p)")

# Execute both channels on J/psi data
alg_etagg.execute_on([data_jpsi, incMC_jpsi, exMC_etagg])
alg_eta3pi.execute_on([data_jpsi, incMC_jpsi, exMC_eta3pi])