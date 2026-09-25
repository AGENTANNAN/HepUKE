### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding ψ(3686) inclusive MC

# Decay card for the signal process ψ(3686) → γ χ_c0, χ_c0 → φ φ η, φ → K+K−, η → γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0        PHSP;
    Enddecay

    Decay chi_c0
    1.0000 phi phi eta         PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-               VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma         PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for ψ(3686) → γ χ_c0, χ_c0 → φ φ η → γ K+K− K+K− γγ
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic0_phiphieta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "ChiCJPhiPhiEta"
chi_cJ_alg = Algorithm.new(alg_name)
chi_cJ_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                        # exactly 2 positive and 2 negative charged tracks
    cos_theta 0.93                       # |cosθ| < 0.93
    Vz 10.0                              # |Vz| < 10 cm
    Vr 1.0                               # Vr < 1 cm
    nChrp "==2"
    nChrn "==2"
    nNet "==0"                           # net charge zero
  }
  .select_photon {                       # at least 3 photons
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0                  # each ≥ 10° away from any charged track
    energyThreshold_b 0.025              # E > 25 MeV (barrel)
    energyThreshold_e 0.050              # E > 50 MeV (endcap)
    nGam ">=3"
  }
  .pid(method: :probability) {           # kaon identification, 0.001 probability cut
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # all charged tracks treated as kaons
    nkp "==2"
    nkm "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # η → γγ reconstruction (mass-constrained)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :kp, :km, :eta]) {   # 4C fit to γ K+K−K+K−η
    nominal
    constrain_four_momentum                              # total four-momentum = CMS
    invariant_mass_of(:kp, :km).constrain_to_nominal_mass_of(:phi)   # two φ candidates
    chi2_cut 200                                         # loose cut; optimal cut applied in ROOT
  }

# Procedures that cannot be expressed as DSL constructs are preserved as notes.
chi_cJ_alg
  .note(:combinatorial_selection, "the photon not used in the η is identified as the radiative photon; the two φ candidates are chosen from the K+K− pairs by minimizing ΔM²=(M(K+K−)_1−m_φ)²+(M(K+K−)_2−m_φ)², each within M(K+K−) ∈ [1.005, 1.035] GeV/c²; the η is the γγ pair closest to the nominal η mass within [0.52, 0.58] GeV/c²")
  .note(:background_veto, "π⁰ veto: reject events with any γγ combination having M(γγ) ∈ [0.115, 0.150] GeV/c²; M(γη) veto: reject events with M(γη) ∈ [1.00, 1.04] GeV/c²; η-recoil-mass cut: M_recoil(η) < 3.05 GeV/c²")
  .note(:signal_regions, "χ_cJ signal regions on the recoil mass against the radiative photon: χ_c0 ∈ [3.38, 3.45], χ_c1 ∈ [3.48, 3.54], χ_c2 ∈ [3.54, 3.60] GeV/c² (applied at ROOT level)")
  .note(:helix_correction, "helix-parameter track corrections derived from a K_S⁰ control sample applied to charged tracks before the kinematic fit")
  .note(:efficiency_curve, "E1 radiative-transition weighting factor (Eγ/Eγ0)³ applied for efficiency correction")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
chi_cJ_alg.execute_on([psip_data, psip_incMC, exMC_signal])