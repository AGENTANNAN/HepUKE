# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding ψ(3686) inclusive MC

# ============================================================
# Decay cards (EvtGen format)
# ============================================================

# ---- Signal: ψ(3686) → γ χ_cJ → γγ J/ψ (J = 0,1,2) ----
# J/ψ → e+e-
decay_card_chicJ_ee = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3333 gamma chi_c2 P2GC2;
    Enddecay
    Decay chi_c0
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# J/ψ → μ+μ-
decay_card_chicJ_mumu = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3333 gamma chi_c2 P2GC2;
    Enddecay
    Decay chi_c0
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# ---- Signal: ψ(3686) → γ η_c(2S) → γγ J/ψ ----
decay_card_etac2S_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay
    Decay eta_c(2S)
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_etac2S_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay
    Decay eta_c(2S)
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# ---- Background: ψ(3686) → π0 J/ψ ----
decay_card_pi0Jpsi_ee = <<~DECAYCARD
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

decay_card_pi0Jpsi_mumu = <<~DECAYCARD
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

# ---- Background: ψ(3686) → π0 π0 J/ψ ----
decay_card_2pi0Jpsi_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 pi0 J/psi PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_2pi0Jpsi_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 pi0 J/psi PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# ---- Background: nonresonant ψ(3686) → γγ J/ψ ----
decay_card_nrJpsi_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_nrJpsi_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# ============================================================
# Exclusive MC samples
# ============================================================

# γχ_cJ → γγ J/ψ : 1,000,000 events
exMC_chicJ_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_chicJ_gammagammaJpsi_ee"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_chicJ_ee
  config.cross_section   = :default
end
exMC_chicJ_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_chicJ_gammagammaJpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_chicJ_mumu
  config.cross_section   = :default
end

# γη_c(2S) → γγ J/ψ : 500,000 events
exMC_etac2S_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_etac2S_gammagammaJpsi_ee"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_etac2S_ee
  config.cross_section   = :default
end
exMC_etac2S_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_etac2S_gammagammaJpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_etac2S_mumu
  config.cross_section   = :default
end

# ψ(3686) → π0 J/ψ : 500,000 events
exMC_pi0Jpsi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_pi0Jpsi_ee"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_pi0Jpsi_ee
  config.cross_section   = :default
end
exMC_pi0Jpsi_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_pi0Jpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_pi0Jpsi_mumu
  config.cross_section   = :default
end

# ψ(3686) → π0π0 J/ψ : 500,000 events
exMC_2pi0Jpsi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_2pi0Jpsi_ee"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_2pi0Jpsi_ee
  config.cross_section   = :default
end
exMC_2pi0Jpsi_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_2pi0Jpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_2pi0Jpsi_mumu
  config.cross_section   = :default
end

# nonresonant ψ(3686) → γγ J/ψ : 500,000 events
exMC_nrJpsi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_nonres_gammagammaJpsi_ee"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_nrJpsi_ee
  config.cross_section   = :default
end
exMC_nrJpsi_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_nonres_gammagammaJpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_nrJpsi_mumu
  config.cross_section   = :default
end

# ============================================================
# Event selection (BOSS)  — common γγ l+l- selection
# ============================================================
alg_name = "GammaGammaLepJpsi"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})          # ψ(3686) center-of-mass energy (GeV)
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                       # charged-track selection
      cos_theta 0.93                    # |cosθ| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm
      momentum  "> 1.0"                 # p > 1.0 GeV/c
      nChrp     "==1"                   # exactly one positive track
      nChrn     "==1"                   # exactly one negative track
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # photon selection
      tdc_emc_start     0               # EMC TDC window start
      tdc_emc_end       14              # EMC TDC window end
      angle_to_track    10.0            # ≥ 10° from any charged track
      energyThreshold_b 0.025           # barrel energy threshold 25 MeV
      energyThreshold_e 0.025           # endcap energy threshold 25 MeV
      nGam              ">=2"           # 2–4 photons (upper bound applied at ROOT level)
  }
  .pid(method: :probability) {          # particle identification
      prob_cut 0.001                    # PID probability > 0.001
      # tracks with p > 1.0 GeV/c → lepton; lepton with EMC E > 0.6 GeV → e, otherwise μ
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlp ">=1"                         # at least one positive lepton
      nlm ">=1"                         # at least one negative lepton
  }
  # 4C kinematic fit on γγ l+l- ; DSL automatically picks the two-photon
  # combination with the smallest χ². χ² < 60, flagged nominal.
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
      nominal
      constrain_four_momentum                         # 4C constraint to ψ(3686) four-momentum
      invariant_mass_of(:lp, :lm).within(3.08, 3.12)  # J/ψ window 3.08–3.12 GeV/c²
      invariant_mass_of(:gamma, :gamma).out_of(0.11, 0.15)  # veto π0 → γγ (0.11–0.15 GeV/c²)
      invariant_mass_of(:gamma, :gamma).out_of(0.51, 0.57)  # veto η  → γγ (0.51–0.57 GeV/c²)
      chi2_cut 60
  }

# Radiative Bhabha suppression (e+e- mode only): cosθ(e+) < 0.3 and cosθ(e-) > -0.3.
# Not expressible on the combined lepton lists inside this BOSS selection.
alg.note(:background_veto,
  "radiative Bhabha suppression for the e+e- mode: require cos(theta_e+) < 0.3 and " \
  "cos(theta_e-) > -0.3; applied on the lepton candidates (e+e- channel only)")

# Render the algorithm from the signal decay card and execute on all datasets.
alg.with_decay_card(decay_card_chicJ_ee).apply(event_selection)

root_files = alg.execute_on([
  psip_data, psip_incMC,
  exMC_chicJ_ee,  exMC_chicJ_mumu,
  exMC_etac2S_ee, exMC_etac2S_mumu,
  exMC_pi0Jpsi_ee, exMC_pi0Jpsi_mumu,
  exMC_2pi0Jpsi_ee, exMC_2pi0Jpsi_mumu,
  exMC_nrJpsi_ee, exMC_nrJpsi_mumu
])