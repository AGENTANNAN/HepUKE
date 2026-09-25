# ============================================================
# Datasets
# ============================================================
# ψ(2S) real data at 3.686 GeV (106.4M ψ(2S) events) and the corresponding inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
# 3.650 GeV continuum (off-resonance) real data used for QED / continuum background
psip_cont  = DatasetManager.real_data.find("709_3650")

# ============================================================
# Decay cards (EvtGen) — one per signal mode
# ============================================================
# ψ(2S) → π0 J/ψ, π0 → γγ, J/ψ → e+e-
decay_card_pi0_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 J/psi PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ψ(2S) → π0 J/ψ, π0 → γγ, J/ψ → μ+μ-
decay_card_pi0_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 J/psi PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ψ(2S) → η J/ψ, η → γγ, J/ψ → e+e-
decay_card_eta_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 eta J/psi PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ψ(2S) → η J/ψ, η → γγ, J/ψ → μ+μ-
decay_card_eta_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 eta J/psi PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ============================================================
# Exclusive MC — 1,000,000 events for each of the four modes
# ============================================================
exMC_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_jpsi_ee"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pi0_ee
  config.cross_section   = :default
end

exMC_pi0_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_jpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pi0_mumu
  config.cross_section   = :default
end

exMC_eta_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_eta_jpsi_ee"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_ee
  config.cross_section   = :default
end

exMC_eta_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_eta_jpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_mumu
  config.cross_section   = :default
end

# ============================================================
# Shared event selection — the π0 and η channels use one common chain
# ============================================================
event_selection = Selection.new
  .select_track {              # Exactly two charged tracks: one positive, one negative, net charge zero
    cos_theta 0.93             # |cosθ| < 0.93
    Vz        10.0             # |Vz| < 10 cm (beam direction)
    Vr        1.0              # Vr < 1 cm (transverse plane)
    nChrp     "==1"            # Exactly 1 positively charged track
    nChrn     "==1"            # Exactly 1 negatively charged track
    nNet      "==0"            # Net charge zero
  }
  .select_photon {             # At least two photons (π0/η → γγ)
    tdc_emc_start     0        # EMC timing window start
    tdc_emc_end       14       # EMC timing window end
    angle_to_track    10.0     # Angle to nearest charged track > 10°
    energyThreshold_b 0.025    # 25 MeV (barrel)
    energyThreshold_e 0.050    # 50 MeV (endcap)
    nGam              ">=2"    # At least 2 photons
  }
  .pid(method: :probability) { # Lepton identification via the probability method
    # High-momentum tracks (p > 1.0 GeV) are treated as leptons; an electron if the EMC
    # energy exceeds 0.6 GeV, otherwise a muon
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"                  # Exactly 1 l+
    nlm "==1"                  # Exactly 1 l-
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) { # 4C fit to γγ l+l- with beam four-momentum constraint
    nominal                    # Nominal fit — corrected four-momenta are saved
    constrain_four_momentum    # Constrain total four-momentum to the beam four-momentum
    chi2_cut 200               # Loose χ² < 200 cut (tight cut performed in ROOT)
  }

# ============================================================
# Algorithms — one per meson channel (π0, η), each covering both lepton flavours
# ============================================================
# --- π0 channel: ψ(2S) → π0 J/ψ →
# ---   π0 → γγ with J/ψ → e+e- (exMC_pi0_ee) or J/ψ → μ+μ- (exMC_pi0_mumu)
alg_name_pi0 = "PsipPi0Jpsi"
alg_pi0 = Algorithm.new(alg_name_pi0)
alg_pi0.set_header(["#{alg_name_pi0}Alg/#{alg_name_pi0}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:pid_correction_method,
             "Lepton identification applies an E/p criterion in addition to the " \
             "probability-method high-momentum lepton selection: muons require " \
             "0.08 < E/p < 0.22 and electrons require E/p > 0.8, where E is the " \
             "associated EMC cluster energy and p is the MDC track momentum.")
alg_pi0.with_decay_card(decay_card_pi0_ee).apply(event_selection.dup)
root_files_pi0 = alg_pi0.execute_on([psip_data, psip_cont, psip_incMC, exMC_pi0_ee, exMC_pi0_mumu])

# --- η channel: ψ(2S) → η J/ψ →
# ---   η → γγ with J/ψ → e+e- (exMC_eta_ee) or J/ψ → μ+μ- (exMC_eta_mumu)
alg_name_eta = "PsipEtaJpsi"
alg_eta = Algorithm.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:pid_correction_method,
             "Lepton identification applies an E/p criterion in addition to the " \
             "probability-method high-momentum lepton selection: muons require " \
             "0.08 < E/p < 0.22 and electrons require E/p > 0.8, where E is the " \
             "associated EMC cluster energy and p is the MDC track momentum.")
alg_eta.with_decay_card(decay_card_eta_ee).apply(event_selection.dup)
root_files_eta = alg_eta.execute_on([psip_data, psip_cont, psip_incMC, exMC_eta_ee, exMC_eta_mumu])