# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# Decay card for the signal chain: ψ(3686) → γ η_c, η_c → p p̄ η, η → γγ
# (final state p p̄ γγγ)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c       PHSP;
    Enddecay

    Decay eta_c
    1.000 p+ anti-p- eta    PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal chain: 500k events, default cross section
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gamma_etac_ppbar_eta"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaEtaCToPPbarEta"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  .select_track {
      cos_theta 0.93     # |cosθ| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr        1.0      # Vr < 1 cm
      nChrp     "==1"    # exactly one positive track
      nChrn     "==1"    # exactly one negative track
      nNet      "==0"    # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0    # more than 10° from any charged track
      energyThreshold_b 0.025   # 25 MeV (EMC barrel)
      energyThreshold_e 0.050   # 50 MeV (EMC endcap)
      nGam              ">=3"   # at least three photons
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]   # identify p and p̄ at once
      nprp "==1"                                  # exactly one proton
      nprm "==1"                                  # exactly one anti-proton
  }
  .select_isolated_photon {
      angle_to_prp_track 20.0   # > 20° from the proton track
      angle_to_prm_track 20.0   # > 20° from the anti-proton track
      nGam               ">=3"  # at least three photons remain
  }
  # Nominal 4C kinematic fit to p p̄ γγγ
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma]) {
      nominal                  # only this fit defines the corrected four-momenta
      constrain_four_momentum
      chi2_cut 200             # loose cut; tight χ²(3γ) < 65 applied in ROOT
  }
  # Competing 2γ hypothesis — no chi2_cut / no nominal, χ² stored for background suppression
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
      constrain_four_momentum
  }
  # Competing 4γ hypothesis — no chi2_cut / no nominal, χ² stored for background suppression
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma]) {
      constrain_four_momentum
  }

# Background-suppression procedure that has no dedicated DSL method:
# η candidate choice and the associated mass vetoes/windows.
alg.note(:background_veto, "η candidate formed from the two lowest-energy photons with 0.4 < M(γγ) < 0.7 GeV/c² (otherwise the next photon pair); π⁰ veto window 0.12–0.15 GeV/c²; ψ(3686) → η J/ψ recoil-mass veto window 3.072–3.125 GeV/c²; η signal window (0.520, 0.570) GeV/c² with sidebands (0.475–0.500) ∪ (0.600–0.625) GeV/c²")

alg.note(:background_veto, "χ²(3γ) < 65 and χ²(3γ) required below both the 2γ and 4γ alternative-hypothesis χ²; applied at analysis level using the stored chi2 values")

alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on([psip_data, psip_incMC, exMC_signal])