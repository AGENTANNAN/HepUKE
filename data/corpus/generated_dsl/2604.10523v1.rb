# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")        # ψ(3686) real data (BOSS 709, 3.686 GeV)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # ψ(3686) inclusive MC

# ---- Decay cards: ψ(2S) → γ χ_cJ ; χ_cJ → ρ+ρ− (VSS) ; ρ± → π±π0 (VSS) ; π0 → γγ (PHSP) ----
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 rho+ rho- VSS;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 VSS;
    Enddecay

    Decay rho-
    1.000 pi- pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 rho+ rho- VSS;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 VSS;
    Enddecay

    Decay rho-
    1.000 pi- pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.000 rho+ rho- VSS;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 VSS;
    Enddecay

    Decay rho-
    1.000 pi- pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC: 500k events for each of the three χ_cJ modes ----
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic0_2rho"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic1_2rho"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic2_2rho"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# The three χ_cJ modes share identical final state (γ π+π− π0π0) and identical
# selection criteria, so a single Algorithm / Selection chain serves all of them.
alg_name = "GammaChiCJTo2Rho"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {
      cos_theta 0.93        # |cosθ| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nChrp     "==1"       # exactly one positive track
      nChrn     "==1"       # exactly one negative track
      nNet      "==0"       # net charge zero
  }
  .select_photon {
      tdc_emc_start     0      # EMC timing start 0 ns
      tdc_emc_end       14     # EMC timing end 700 ns (14 × 50 ns)
      angle_to_track    10.0   # > 10° to nearest charged track
      energyThreshold_b 0.025  # 25 MeV (barrel)
      energyThreshold_e 0.050  # 50 MeV (endcap)
      nGam              ">=5"  # radiative γ + four γ from the two π0
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # π+ / π− separated from K and p
      npip "==1"
      npim "==1"
  }
  # Reconstruct the two π0 from photon pairs: 1-C mass fit each, χ² < 25
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=2"
  }
  # 6C fit: four-momentum constraint plus the two π0 mass constraints (π0 pre-fitted above)
  .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) {
      nominal
      constrain_four_momentum
      chi2_cut 200          # loose BOSS-level cut; tighter χ² < 60 applied later in ROOT
  }
  # Competing hypothesis with one extra photon (χ²_more) — four-momentum only, χ² stored
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0, :pi0]) {
      constrain_four_momentum
  }
  # Competing hypothesis with one photon less (χ²_less) — four-momentum only, χ² stored
  .kinematic_fit([:pip, :pim, :pi0, :pi0]) {
      constrain_four_momentum
  }

# BOSS-side procedures that cannot be expressed as formal DSL constructs
algorithm
  .note(:helix_correction, "charged-track helix-parameter correction applied to all
    tracks before the 6C kinematic fit; the efficiency difference with/without the
    correction is evaluated by re-running the BOSS selection")
  .note(:background_veto, "additional vetoes applied in the analysis: π0π0 J/ψ recoil
    region M_recoil(π0π0) > 3.07 GeV/c²; π+π− J/ψ recoil region
    M_recoil(π+π−) ∈ [3.09, 3.10] GeV/c²; η/J/ψ region
    3.09 < M(π+π−π0) < 3.18 GeV/c²; and a 2D K_S K_S window")
  .with_decay_card(decay_card_chic1)   # header defines the shared γ π+π− π0π0 kinematics
  .apply(event_selection)

# Execute on real data, inclusive MC and the three signal-MC samples
root_files = algorithm.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])