# =============================================================================
# BESIII  psi(3686) -> gamma eta_c(2S)
#   eta_c(2S) -> K_S0 K^± pi^∓   (charge-conjugate modes)
#   eta_c(2S) -> K^+ K^- pi^0
# BOSS part only: dataset preparation + event selection up to the kinematic fit.
# =============================================================================

### Dataset description ###
psip_data       = DatasetManager.real_data.find("709_3686")     # psi(3686) data @ 3.686 GeV (~106 M events, 156 pb^-1)
continuum_data  = DatasetManager.real_data.find("709_3650")     # 3.65 GeV continuum data (~42 pb^-1)
psip_incMC      = DatasetManager.inclusive_mc.find("709_3686")  # inclusive MC @ 3.686 GeV
continuum_incMC = DatasetManager.inclusive_mc.find("709_3650")  # inclusive MC @ 3.65 GeV

### Decay cards (EvtGen) ###
# ---- Signal modes ----
# Mode I : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+ pi-
decay_card_kskp = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S) PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K- pi+  (charge conjugate)
decay_card_kskm = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S) PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 K_S0 K- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode III : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K+ K- pi0
decay_card_kkpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S) PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Background modes (exclusive MC) ----
# Fake-photon K_S0 K pi : psi(3686) -> K_S0 K+ pi- (a hadron/shower fakes the signal photon)
decay_card_bg_fake_kskpi = <<~DECAYCARD
  Decay psi(2S)
  1.000 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Fake-photon K K pi0 : psi(3686) -> K+ K- pi0
decay_card_bg_fake_kkpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# FSR K_S0 K pi : psi(3686) -> K_S0 K+ pi- with final-state radiation (PHOTOS)
decay_card_bg_fsr_kskpi = <<~DECAYCARD
  Decay psi(2S)
  1.000 K_S0 K+ pi- PHOTOS PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# omega K+ K-
decay_card_bg_omega_kk = <<~DECAYCARD
  Decay psi(2S)
  1.000 omega K+ K- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# pi0 K_S0 K pi
decay_card_bg_pi0_kskpi = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# pi0 K+ K- pi0
decay_card_bg_pi0_kkpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# pi+ pi- J/psi , J/psi -> K+ K-
decay_card_bg_pipi_jpsi_kk = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 K+ K- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# pi+ pi- J/psi , J/psi -> leptons
decay_card_bg_pipi_jpsi_ll = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# eta J/psi , eta -> 3 pi (pi+ pi- pi0)
decay_card_bg_eta_jpsi_3pi = <<~DECAYCARD
  Decay psi(2S)
  1.000 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# eta J/psi , eta -> gamma gamma
decay_card_bg_eta_jpsi_2g = <<~DECAYCARD
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

# Continuum gamma* -> K_S0 K pi (top mother set to psi(4260) per BESIII KKMC convention)
decay_card_bg_continuum = <<~DECAYCARD
  Decay psi(4260)
  1.000 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples (100k events each) ###
# ---- Signal ----
exMC_sig_kskp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2s_kskp"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kskp
  config.cross_section   = :default
end

exMC_sig_kskm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2s_kskm"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kskm
  config.cross_section   = :default
end

exMC_sig_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2s_kkpi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kkpi0
  config.cross_section   = :default
end

# ---- Backgrounds (resonant, at the psi(3686)) ----
exMC_bg_fake_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_fake_kskpi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_fake_kskpi
  config.cross_section   = :default
end

exMC_bg_fake_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_fake_kkpi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_fake_kkpi0
  config.cross_section   = :default
end

exMC_bg_fsr_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_fsr_kskpi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_fsr_kskpi
  config.cross_section   = :default
end

exMC_bg_omega_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_omega_kk"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_omega_kk
  config.cross_section   = :default
end

exMC_bg_pi0_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_pi0_kskpi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_pi0_kskpi
  config.cross_section   = :default
end

exMC_bg_pi0_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_pi0_kkpi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_pi0_kkpi0
  config.cross_section   = :default
end

exMC_bg_pipi_jpsi_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_pipi_jpsi_kk"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_pipi_jpsi_kk
  config.cross_section   = :default
end

exMC_bg_pipi_jpsi_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_pipi_jpsi_ll"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_pipi_jpsi_ll
  config.cross_section   = :default
end

exMC_bg_eta_jpsi_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_eta_jpsi_3pi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_eta_jpsi_3pi
  config.cross_section   = :default
end

exMC_bg_eta_jpsi_2g = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_eta_jpsi_2g"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_eta_jpsi_2g
  config.cross_section   = :default
end

# ---- Background (continuum, at 3.65 GeV) ----
exMC_bg_continuum = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_bg_continuum_kskpi"
  config.related_dataset = continuum_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_continuum
  config.cross_section   = :default
end

# =============================================================================
### Event selection (BOSS) ###
# -----------------------------------------------------------------------------
# Algorithm A1 : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+ pi-
# -----------------------------------------------------------------------------
alg_name_kskp = "EtaC2SToKsKPiP"
alg_kskp = Algorithm.new(alg_name_kskp)
alg_kskp.set_header(["#{alg_name_kskp}Alg/#{alg_name_kskp}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

selection_kskp = Selection.new
  .select_track {                       # charged track quality + multiplicity
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==2"                     # exactly two positive tracks (pi+ from K_S0, K+)
    nChrn     "==2"                     # exactly two negative tracks (pi- from K_S0, pi-)
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon (shower) selection
    tdc_emc_start     0
    tdc_emc_end       14                # EMC timing in [0, 14]
    energyThreshold_b 0.025             # E > 25 MeV in the barrel
    energyThreshold_e 0.050             # E > 50 MeV in the endcap
    angle_to_track    10.0              # at least 10 deg from any charged track
    nGam              ">=1"             # at least one photon (the signal gamma)
  }
  .pid(method: :probability) {          # PID: kaon identification against pi/p
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp      ">=1"                      # the K+ of this mode must be identified
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])          # take the identified kaon out
  .assign({:chrgp => :pip, :chrgn => :pim})        # remaining tracks are pions
  .secondary_vertex_fit([:pip, :pim]) {            # build K_S0 from pi+pi- pairs
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list       # used pions left the two tracks of "K pi"
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {     # 4C fit of gamma K_S0 K+ pi-
    nominal
    constrain_four_momentum
    chi2_cut 200                                   # loose; paper requires chi2_4C < 50 in ROOT
  }

alg_kskp
  .note(:ks0_mass_window, "K_S0 candidates required to satisfy |M(pi+pi-) - m(K_S0)| < 7 MeV/c^2; " \
        "candidate chosen among all opposite-charge pairs by the best secondary-vertex fit")
  .note(:ks0_vertex_distance, "K_S0 decay vertex required to be at least 0.5 cm from the interaction point")
  .note(:ks0_exclusivity, "the remaining tracks must not form another good K_S0 candidate")
  .note(:recoil_mass_cut, "recoil mass of all pi+pi- pairs required < 3.05 GeV/c^2 (suppresses " \
        "psi(3686) and continuum backgrounds; recoil-mass helper not expressible in the current DSL)")
  .note(:chi2_combination, "best candidate combination selected by chi2_4C + chi2_PID(K) + chi2_PID(pi); " \
        "default kinematic_fit candidate ordering uses chi2_4C only")
  .note(:background_veto, "invariant mass of the two charged tracks (K pi) evaluated under the muon mass " \
        "hypothesis required < 2.9 GeV/c^2 to veto J/psi and eta J/psi backgrounds")
  .with_decay_card(decay_card_kskp)
  .apply(selection_kskp)

# -----------------------------------------------------------------------------
# Algorithm A2 : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K- pi+  (charge conjugate)
# -----------------------------------------------------------------------------
alg_name_kskm = "EtaC2SToKsKPiM"
alg_kskm = Algorithm.new(alg_name_kskm)
alg_kskm.set_header(["#{alg_name_kskm}Alg/#{alg_name_kskm}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

selection_kskm = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkm      ">=1"                      # the K- of this mode must be identified
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :km, :pip]) {     # 4C fit of gamma K_S0 K- pi+
    nominal
    constrain_four_momentum
    chi2_cut 200                                   # loose; paper requires chi2_4C < 50 in ROOT
  }

alg_kskm
  .note(:ks0_mass_window, "K_S0 candidates required to satisfy |M(pi+pi-) - m(K_S0)| < 7 MeV/c^2; " \
        "candidate chosen among all opposite-charge pairs by the best secondary-vertex fit")
  .note(:ks0_vertex_distance, "K_S0 decay vertex required to be at least 0.5 cm from the interaction point")
  .note(:ks0_exclusivity, "the remaining tracks must not form another good K_S0 candidate")
  .note(:recoil_mass_cut, "recoil mass of all pi+pi- pairs required < 3.05 GeV/c^2 (suppresses " \
        "psi(3686) and continuum backgrounds; recoil-mass helper not expressible in the current DSL)")
  .note(:chi2_combination, "best candidate combination selected by chi2_4C + chi2_PID(K) + chi2_PID(pi); " \
        "default kinematic_fit candidate ordering uses chi2_4C only")
  .note(:background_veto, "invariant mass of the two charged tracks (K pi) evaluated under the muon mass " \
        "hypothesis required < 2.9 GeV/c^2 to veto J/psi and eta J/psi backgrounds")
  .with_decay_card(decay_card_kskm)
  .apply(selection_kskm)

# -----------------------------------------------------------------------------
# Algorithm B : psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K+ K- pi0 (--> gamma gamma)
# -----------------------------------------------------------------------------
alg_name_kkpi0 = "EtaC2SToKKPi0"
alg_kkpi0 = Algorithm.new(alg_name_kkpi0)
alg_kkpi0.set_header(["#{alg_name_kkpi0}Alg/#{alg_name_kkpi0}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

selection_kkpi0 = Selection.new
  .select_track {
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"                     # exactly one positive track (K+)
    nChrn     "==1"                     # exactly one negative track (K-)
    nNet      "==0"                     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"             # two photons from pi0 + one signal gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp      "==1"                      # exactly one K+
    nkm      "==1"                      # exactly one K-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # 1C mass-constrained pi0 reconstruction
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pi0]) {         # 4C (total 5C with the pi0 mass constraint)
    nominal
    constrain_four_momentum
    chi2_cut 200                                     # loose; paper requires chi2_5C < 30 in ROOT
  }

alg_kkpi0
  .note(:photon_pi0_selection, "the best recoil photon and the best pi0 candidate are chosen by the " \
        "smallest chi2_5C of the gamma K+ K- pi0 fit")
  .note(:background_veto, "invariant mass of the two charged tracks (K+ K-) evaluated under the muon mass " \
        "hypothesis required < 2.9 GeV/c^2 to veto J/psi and eta J/psi backgrounds")
  .with_decay_card(decay_card_kkpi0)
  .apply(selection_kkpi0)

# =============================================================================
### Execute on datasets (real data, inclusive MC, signal & background exclusive MC) ###
### =============================================================================
boss_datasets = [
  psip_data, continuum_data,
  psip_incMC, continuum_incMC,
  exMC_sig_kskp, exMC_sig_kskm, exMC_sig_kkpi0,
  exMC_bg_fake_kskpi, exMC_bg_fake_kkpi0, exMC_bg_fsr_kskpi, exMC_bg_omega_kk,
  exMC_bg_pi0_kskpi, exMC_bg_pi0_kkpi0,
  exMC_bg_pipi_jpsi_kk, exMC_bg_pipi_jpsi_ll,
  exMC_bg_eta_jpsi_3pi, exMC_bg_eta_jpsi_2g,
  exMC_bg_continuum
]

root_files_kskp  = alg_kskp.execute_on(boss_datasets)
root_files_kskm  = alg_kskm.execute_on(boss_datasets)
root_files_kkpi0 = alg_kkpi0.execute_on(boss_datasets)