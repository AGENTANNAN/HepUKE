# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data (3.686 GeV)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the main tag mode:
#   Xi(1530)- -> pi- Xi0 ; Xi0 -> pi0 Lambda ; Lambda -> p pi- ; pi0 -> gamma gamma
decay_card_pim_xi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi(1530)- anti-Xi(1530)+ PHSP;
    Enddecay

    Decay Xi(1530)-
    1.0000 pi- Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 pi0 Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the alternative tag mode (requires its own card):
#   Xi(1530)- -> pi0 Xi- ; Xi- -> pi- Lambda ; Lambda -> p pi- ; pi0 -> gamma gamma
decay_card_pi0_xim = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi(1530)- anti-Xi(1530)+ PHSP;
    Enddecay

    Decay Xi(1530)-
    1.0000 pi0 Xi- PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for psi(3686) -> Xi(1530)- anti-Xi(1530)+ (full tag decay chain)
exMC_pim_xi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_Xi1530_Xi1530bar_pimXi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_pim_xi0
  config.cross_section   = :default
end

# Exclusive MC for the pi0 Xi- tag mode
exMC_pi0_xim = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_Xi1530_Xi1530bar_pi0Xim"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_pi0_xim
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Selection common to both tag modes (identical final state: p, pi-, pi-, pi0 -> gamma gamma)
common_selection = Selection.new
  .select_track {                    # charged track selection
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        100.0                  # |Vz| < 100 cm
    Vr        10.0                   # Vr < 10 cm
    nChrp     ">=1"                  # at least one positive track
    nChrn     ">=2"                  # at least two negative tracks
  }
  .select_photon {                   # photon selection
    tdc_emc_start     0              # TDC window [0, 14] (= 0-700 ns)
    tdc_emc_end       14
    energyThreshold_b 0.025          # 25 MeV (barrel)
    energyThreshold_e 0.050          # 50 MeV (endcap)
    nGam              ">=2"          # at least two photons
  }
  .pid(method: :probability) {       # PID by the probability method
    prob_cut 0.001                   # probability > 0.001
    identify :proton, against: [:kaon, :pion]     # p / anti-p vs K, pi
    identify :pion,   against: [:kaon, :proton]   # pi+ / pi- vs K, p
    nprp ">=1"                       # at least one proton
    npim ">=2"                       # at least two pi-
  }
  .remove([:prp <= :chrgp])          # remove protons from the positive charged list
  .assign({:chrgp => :pip})          # remaining positive tracks assigned as pi+
  .kalman_kinematic_fit([:gamma, :gamma]) {      # 1C fit of the two photons to the pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20                      # chi^2 < 20
    npi0 ">=1"                       # at least one pi0 candidate
  }
  .secondary_vertex_fit([:prp, :pim]) {           # build Lambda from p pi- (secondary vertex)
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }

# --- Tag mode pi- Xi0: reconstruct the tag chain, infer the anti-baryon as recoil ---
# recIDs: 0 psi(2S), 1 Xi(1530)-, 2 anti-Xi(1530)+, 3 pi-, 4 Xi0, 5 pi0, 6 Lambda0,
#         7 gamma, 8 gamma, 9 p, 10 pi-
selection_pim_xi0 = common_selection.dup
  .partial_rec([1, 3, 4, 5, 6, 7, 8, 9, 10]) {
    best_combination_by_mass :Xi0,    1.31486     # Xi0 nominal mass (GeV/c^2)
    best_combination_by_mass :Xi1530, 1.535       # Xi(1530)- nominal mass (GeV/c^2)
    require_recoil_mass 1.20, 1.70                # recoil (anti-baryon) mass window (GeV/c^2)
  }

# --- Tag mode pi0 Xi- ---
# recIDs: 0 psi(2S), 1 Xi(1530)-, 2 anti-Xi(1530)+, 3 pi0, 4 Xi-, 5 gamma, 6 gamma,
#         7 pi-, 8 Lambda0, 9 p, 10 pi-
selection_pi0_xim = common_selection.dup
  .partial_rec([1, 3, 4, 5, 6, 7, 8, 9, 10]) {
    best_combination_by_mass :Xi_minus, 1.32171   # Xi- nominal mass (GeV/c^2)
    best_combination_by_mass :Xi1530,   1.535     # Xi(1530)- nominal mass (GeV/c^2)
    require_recoil_mass 1.20, 1.70                # recoil (anti-baryon) mass window (GeV/c^2)
  }

# Algorithm for the pi- Xi0 tag mode
alg_pim_xi0 = Algorithm.new("Xi1530TagPimXi0")
alg_pim_xi0.set_header(["Xi1530TagPimXi0Alg/Xi1530TagPimXi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:mass_windows, "Xi0/Xi- mass window 10 MeV/c^2, Lambda mass window 5 MeV/c^2, Xi(1530)- mass window 15 MeV/c^2; applied in the ROOT-level analysis")
           .note(:tag_mode_pi0_xim, "the pi0 Xi- tag mode requires a separate decay card; both tag modes are combined in the full analysis")
           .with_decay_card(decay_card_pim_xi0).apply(selection_pim_xi0)

# Algorithm for the pi0 Xi- tag mode (separate decay card)
alg_pi0_xim = Algorithm.new("Xi1530TagPi0Xim")
alg_pi0_xim.set_header(["Xi1530TagPi0XimAlg/Xi1530TagPi0Xim.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:mass_windows, "Xi0/Xi- mass window 10 MeV/c^2, Lambda mass window 5 MeV/c^2, Xi(1530)- mass window 15 MeV/c^2; applied in the ROOT-level analysis")
          .with_decay_card(decay_card_pi0_xim).apply(selection_pi0_xim)

# Execute on real data, inclusive MC and the corresponding signal MC
root_files_pim_xi0 = alg_pim_xi0.execute_on([psip_data, psip_incMC, exMC_pim_xi0])
root_files_pi0_xim = alg_pi0_xim.execute_on([psip_data, psip_incMC, exMC_pi0_xim])