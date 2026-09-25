# DSL: 2503.17165v1 — CP test with Σ+→pπ0, Σ̄-→p̄π0 at J/ψ and ψ(3686)
# BESIII: (1.0087±0.0044)×10^10 J/ψ + (2.7124±0.0143)×10^9 ψ(3686)
# Ordinary analysis — Selection + kinematic fit

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ── Datasets ──────────────────────────────────────────────────────────
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ── Decay cards ───────────────────────────────────────────────────────
# J/ψ → Σ+ Σ̄-; Σ+ → p π0; Σ̄- → p̄ π0
decay_card_jpsi = <<~DECAYCARD
  Decay J/psi
  1.000 Sigma+ Sigma_bar- HELAMP 1.0 0.0 1.0 0.0;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay Sigma_bar-
  1.000 p_bar- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
DECAYCARD

# ψ(3686) → Σ+ Σ̄-; same subsequent decays
decay_card_psip = <<~DECAYCARD
  Decay psi(2S)
  1.000 Sigma+ Sigma_bar- HELAMP 1.0 0.0 1.0 0.0;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay Sigma_bar-
  1.000 p_bar- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
DECAYCARD

# ── Exclusive MC ──────────────────────────────────────────────────────
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Sigma_p_pi0_jpsi"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Sigma_p_pi0_psip"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

# ── Algorithm: J/ψ ────────────────────────────────────────────────────
alg_jpsi = Algorithm.new("SigmaP2PPi0_Jpsi")
alg_jpsi.set_header(["SigmaP2PPi0_JpsiAlg/SigmaP2PPi0_Jpsi.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .with_decay_card(decay_card_jpsi)

event_selection_jpsi = Selection.new
  # Track selection: 2 charged tracks, net charge zero
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        2.0       # paper: < 2 cm in transverse plane
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  # Photon selection: at least 2 photons for π0→γγ
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # Isolated photon selection (anti-proton can fake photons in EMC)
  .select_isolated_photon do
    angle_to_prm_track 20.0
    angle_to_prp_track 20.0
    nGam               ">=2"
  end
  # PID: identify protons and anti-protons
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  end
  # Remove identified protons from generic charged lists, assign remainder as pions
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({ chrgp: :pip, chrgn: :pim })
  # Kalman kinematic fit: reconstruct π0 → γγ
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200         # 1C fit χ² < 200 (paper: χ²_1C < 200)
    npi0 ">=1"
  end
  # 2C kinematic fit: e+e- → p π0 p̄ π0 with four-momentum conservation + π0 mass
  # Paper: 2C fit constraining M(γγ) to π0 mass and four-momentum conservation
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
    constrain_four_momentum
    chi2_cut 30          # paper: χ²_2C < 30
    nominal
  end

alg_jpsi.with_decay_card(decay_card_jpsi).apply(event_selection_jpsi)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

# ── Algorithm: ψ(3686) ────────────────────────────────────────────────
alg_psip = Algorithm.new("SigmaP2PPi0_Psip")
alg_psip.set_header(["SigmaP2PPi0_PsipAlg/SigmaP2PPi0_Psip.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .with_decay_card(decay_card_psip)

event_selection_psip = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        2.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .select_isolated_photon do
    angle_to_prm_track 20.0
    angle_to_prp_track 20.0
    nGam               ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({ chrgp: :pip, chrgn: :pim })
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  end
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
    constrain_four_momentum
    chi2_cut 30
    nominal
  end

alg_psip.with_decay_card(decay_card_psip).apply(event_selection_psip)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])

# ── Notes ─────────────────────────────────────────────────────────────
# - Proton DCA cut (>0.34 cm from IP to suppress prompt tracks) not expressible
#   in select_track; captured as note.
# - Secondary vertex fit for Σ+/Σ̄- reconstruction happens in ROOT analysis stage
#   (angular analysis), not in BOSS event selection.
# - Sideband background subtraction, 5D angular fit, CP asymmetry extraction
#   are ROOT-stage procedures, out of scope.
alg_jpsi.note(:proton_dca_cut, "Proton/anti-proton DCA > 0.34 cm from IP")
alg_psip.note(:proton_dca_cut, "Proton/anti-proton DCA > 0.34 cm from IP")
alg_jpsi.note(:sigma_mass_windows, "Signal region: M(pπ0) ∈ [1.172, 1.200], M(pbarπ0) ∈ [1.167, 1.212] GeV/c²")
alg_psip.note(:sigma_mass_windows, "Signal region: M(pπ0) ∈ [1.172, 1.200], M(pbarπ0) ∈ [1.167, 1.212] GeV/c²")