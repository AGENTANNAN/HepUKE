# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(2S) real data at 3.686 GeV (106 M events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC sample

# Decay cards for ψ(2S) → γ χ_cJ, χ_cJ → Λ Λ̄ π+π− (phase space),
# with Λ → p π− and Λ̄ → p̄ π+. One card per χ_cJ mother.
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0                          PHSP;
    Enddecay

    Decay chi_c0
    1.0000 Lambda0 anti-Lambda0 pi+ pi-          PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                           HypWK;
    Enddecay

    End
DECAYCARD

decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1                          PHSP;
    Enddecay

    Decay chi_c1
    1.0000 Lambda0 anti-Lambda0 pi+ pi-          PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                           HypWK;
    Enddecay

    End
DECAYCARD

decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 Lambda0 anti-Lambda0 pi+ pi-          PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                           HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for the χ_c1 mode (phase-space χ_c1 → Λ Λ̄ π+π−),
# with the card cloned (appropriate χ_cJ mother) for the χ_c0 and χ_c2 modes.
exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c1_LLbarpipi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c0_LLbarpipi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c2_LLbarpipi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# χ_c0, χ_c1 and χ_c2 share the identical final state γ Λ Λ̄ π+π− and identical
# selection, so a single Algorithm instance covers all three modes.
alg_name = "PsipGammaChiCJToLLbarPipi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                       # charged-track selection
    cos_theta 0.93                      # |cosθ| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # |Vr| < 1 cm
    nChrp     ">=2"                     # at least 2 positive tracks
    nChrn     ">=2"                     # at least 2 negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025             # 25 MeV in the EMC barrel
    energyThreshold_e 0.050             # 50 MeV in the EMC endcap
    nGam              ">=1"             # at least one photon
  }
  # No charged PID: each positive track is a candidate for both p and π+,
  # each negative track for both p̄ and π−.
  .assign({:chrgp => :prp, :chrgn => :prm})   # positives as p candidates, negatives as p̄ candidates
  .assign({:chrgp => :pip, :chrgn => :pim})   # positives as π+ candidates, negatives as π− candidates
  .secondary_vertex_fit([:prp, :pim]) {       # Λ → p π−, pick combination minimising mass difference
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # Λ̄ → p̄ π+, pick combination minimising mass difference
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :pip, :pim]) {   # 4C fit to γ Λ Λ̄ π+π−
    nominal                                   # nominal fit — corrected four-momenta are the ones stored
    constrain_four_momentum                   # 4-momentum conservation against the CMS energy
    chi2_cut 200                              # loose BOSS cut; tighter χ² < 80 applied in ROOT
  }

# BOSS-side notes for procedures that cannot be expressed in the DSL.
# These are the ROOT-side background vetoes quoted in the description.
my_algorithm
  .note(:background_veto_jpsi,
        "ROOT-side J/psi veto: reject events whose pi+pi- recoil mass lies inside 3.088-3.108 GeV/c2")
  .note(:background_veto_sigma0,
        "ROOT-side Sigma0 veto: reject events with M(gamma Lambda) or M(gamma Lambda_bar) inside 1.183-1.202 GeV/c2")
  .note(:background_veto_xi,
        "ROOT-side Xi veto: reject events with M(Lambda pi-) or M(Lambda_bar pi+) inside 1.312-1.331 GeV/c2")

# Attach the decay card (identical final state for all three χ_cJ modes) and generate the algorithm.
my_algorithm.with_decay_card(decay_card_chi_c1).apply(event_selection)

# Execute on real data, inclusive MC and the three exclusive signal MC samples.
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_chi_c0, exMC_chi_c1, exMC_chi_c2])