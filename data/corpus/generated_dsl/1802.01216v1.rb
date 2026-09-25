# =============================================================================
# Dataset preparation
# =============================================================================
# Cross-section scan e+e- -> (K+K- , K_S0 K_S0 , pi+pi-) J/psi at sqrt(s) = 4.189-4.600 GeV
# (14 energy points). 703_4260 (4.260 GeV) is used here as the representative point;
# the identical selection chain is applied at every energy point.
data_4260  = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# -----------------------------------------------------------------------------
# Decay cards (EvtGen format) - J/psi -> l+ l-
# -----------------------------------------------------------------------------
# Mode I: e+e- -> K+ K- J/psi
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K+ K- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II: e+e- -> K_S0 K_S0 J/psi
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K_S0 K_S0 J/psi    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode III: e+e- -> pi+ pi- J/psi
decay_card_modeIII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# -----------------------------------------------------------------------------
# Exclusive signal MC: 100k events per mode
# -----------------------------------------------------------------------------
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_KKJpsi_ll"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_KsKsJpsi_ll"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_pipiJpsi_ll"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

# =============================================================================
# Event selection (BOSS)
# =============================================================================
# -----------------------------------------------------------------------------
# Mode I: e+e- -> K+ K- J/psi  (J/psi -> l+ l-)
# -----------------------------------------------------------------------------
alg_name_I = "KKJpsi"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:background_veto, "radiative-Bhabha veto: reject events in which any oppositely charged track pair has cos(opening angle) > 0.98; the same veto is applied identically at every energy point")

sel_modeI = Selection.new
sel_modeI.select_track {            # charged-track quality cuts (common to all modes)
    cos_theta 0.93                  # |cos(theta)| < 0.93
    Vz 10.0                         # |Vz| < 10 cm
    Vr 1.0                          # Vr < 1 cm
    nChrp ">=2"                     # >= 2 positive tracks (K+, l+)
    nChrn ">=2"                     # >= 2 negative tracks (K-, l-)
    nNet  "==0"                     # net charge zero
  }
  .pid(method: :probability) {      # PID: high-momentum leptons + kaons
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion, :proton]   # K+ / K-, separated from pi and p
    nkp ">=1"                       # >= 1 K+
    nkm ">=1"                       # >= 1 K-
    nlp ">=1"                       # >= 1 l+
    nlm ">=1"                       # >= 1 l-
  }
  .kinematic_fit([:kp, :km, :lp, :lm]) {        # 4C kinematic fit K+K-l+l-
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# -----------------------------------------------------------------------------
# Mode II: e+e- -> K_S0 K_S0 J/psi  (J/psi -> l+ l-)
# -----------------------------------------------------------------------------
alg_name_II = "KsKsJpsi"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:background_veto, "radiative-Bhabha veto: reject events in which any oppositely charged track pair has cos(opening angle) > 0.98; the same veto is applied identically at every energy point")
          .note(:ks0_vertex_quality, "each K_S0 candidate required to satisfy secondary-vertex fit chi2 < 100 and flight significance L/sigma > 4")
          .note(:ks0_mass_window, "K_S0 candidates required to have M(pi+pi-) within [471, 524] MeV/c^2")

sel_modeII = Selection.new
sel_modeII.select_track {           # charged-track quality cuts
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"                     # >= 3 positive tracks (2 pi+ from K_S0, l+)
    nChrn ">=3"                     # >= 3 negative tracks (2 pi- from K_S0, l-)
    nNet  "==0"
  }
  .pid(method: :probability) {      # PID: high-momentum leptons only
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"                       # >= 1 l+
    nlm ">=1"                       # >= 1 l-
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])      # keep only leftover tracks as pion candidates
  .assign({:chrgp => :pip, :chrgn => :pim})    # remaining tracks assumed pions for K_S0
  .secondary_vertex_fit([:pip, :pim]) {        # first K_S0 (pi+ pi-)
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {        # second K_S0 (pi+ pi-)
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) {   # 6C fit: 4C + J/psi + K_S0 mass constraints
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    chi2_cut 200
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# -----------------------------------------------------------------------------
# Mode III: e+e- -> pi+ pi- J/psi  (J/psi -> l+ l-)
# -----------------------------------------------------------------------------
alg_name_III = "pipiJpsi"
alg_modeIII = Algorithm.new(alg_name_III)
alg_modeIII.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
           .set_constant({"ECMS" => [:double, 4.260]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:background_veto, "radiative-Bhabha veto: reject events in which any oppositely charged track pair has cos(opening angle) > 0.98; the same veto is applied identically at every energy point")

sel_modeIII = Selection.new
sel_modeIII.select_track {          # charged-track quality cuts (common to all modes)
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"                     # >= 2 positive tracks (pi+, l+)
    nChrn ">=2"                     # >= 2 negative tracks (pi-, l-)
    nNet  "==0"
  }
  .pid(method: :probability) {      # PID: high-momentum leptons + pions
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]  # pi+ / pi-, separated from K and p
    npip ">=1"                      # >= 1 pi+
    npim ">=1"                      # >= 1 pi-
    nlp ">=1"                       # >= 1 l+
    nlm ">=1"                       # >= 1 l-
  }
  .kinematic_fit([:pip, :pim, :lp, :lm]) {     # 4C kinematic fit pi+pi-l+l-
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)

# =============================================================================
# Execute on datasets
# =============================================================================
alg_modeI.execute_on([data_4260, incMC_4260, exMC_modeI])
alg_modeII.execute_on([data_4260, incMC_4260, exMC_modeII])
alg_modeIII.execute_on([data_4260, incMC_4260, exMC_modeIII])