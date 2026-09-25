# Study of η → π+π- l+l- (l=e, μ) via J/ψ → γη
# J/ψ, 1.0087×10^10 events, arXiv:2501.10130v1
# BF + TFF + CP asymmetry + ALP search

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_ee = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- e+ e- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- mu+ mu- PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_eta_pipiee"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card_ee
  config.cross_section = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_eta_pipimumu"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card_mumu
  config.cross_section = :default
end

# Algorithm: J/ψ → γη, η → π⁺π⁻ℓ⁺ℓ⁻
# Hybrid PID: high-momentum leptons for ℓ⁺ℓ⁻, pions for π⁺π⁻
# e/μ separation done by selecting best combined χ² = χ²_4C + Σχ²_PID
alg = Algorithm.new("JpsiGammaEtaPiPiLL")
alg.set_header(["JpsiGammaEtaPiPiLLAlg/JpsiGammaEtaPiPiLL.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

# ── Selection ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==4"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  # Hybrid PID: high-momentum leptons for J/ψ→ℓ⁺ℓ⁻ type,
  # pions for η→π⁺π⁻. e/μ separation by combined χ² = χ²_4C + Σχ²_PID in ROOT.
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
    nlp  "==1"
    nlm  "==1"
  }
  # 4C kinematic fit: γ π⁺π⁻ ℓ⁺ℓ⁻
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg.with_decay_card(decay_card_ee).apply(event_selection)

# Additional cuts applied in ROOT:
# - χ²_sum = χ²_4C + Σχ²_PID < 50 (combined PID+kinematic fit quality)
# - Photon conversion veto: Φ_ee vs R_xy
# - e/μ separation: select hypothesis with smallest χ²_sum
# - Signal region e⁺e⁻: M(π⁺π⁻e⁺e⁻) in [0.53, 0.57] GeV/c²
# - Signal region μ⁺μ⁻: M(π⁺π⁻μ⁺μ⁻) in [0.531, 0.567] GeV/c²
# - ALP search: additional selection on M(l⁺l⁻) spectrum
alg.note(:root_cuts, "ROOT-level: χ²_sum(χ²_4C+Σχ²_PID)<50, photon conversion veto, e/μ separation by min χ²_sum, signal region ee [0.53,0.57] GeV/c², μμ [0.531,0.567] GeV/c²")
alg.note(:analysis_goals, "BF measurement, TFF extraction, CP asymmetry test, ALP search in M(l+l-) spectrum")

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal_ee, exMC_signal_mumu])