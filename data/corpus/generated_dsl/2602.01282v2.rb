# =====================================================================
# Dataset preparation
# =====================================================================
# J/psi (3.097 GeV) real data and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---------------------------------------------------------------------
# Decay cards (EvtGen)
# Signal: J/psi -> Lambda (tag, Lambda -> p pi-)
#                   + [ Lambda-bar p -> K+ pi+ pi- + k pi0 ]
# (the Lambda-bar annihilates on a proton at rest in the beam-pipe
#  cooling-oil layer; the observable final state is written directly)
# ---------------------------------------------------------------------

# mode A : k = 1 -> K+ pi+ pi- pi0
decay_card_k1 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 K+ pi+ pi-             PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# mode B : k = 2 -> K+ pi+ pi- 2 pi0
decay_card_k2 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 K+ pi+ pi- pi0 pi0     PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# mode C : k = 3 -> K+ pi+ pi- 3 pi0
decay_card_k3 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 K+ pi+ pi- pi0 pi0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# mode D : K*(892)+ sub-process -> K*+ pi+ pi-, with K*+ -> K+ pi0
decay_card_kstar = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 K*+ pi+ pi-            PHSP;
    Enddecay

    Decay K*+
    1.0000 K+ pi0                         VSS;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------
# Exclusive MC: 200k events per signal mode
# ---------------------------------------------------------------------
exMC_k1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_annih_kpipipim_pi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_k1
  config.cross_section   = :default
end

exMC_k2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_annih_kpipipim_2pi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_k2
  config.cross_section   = :default
end

exMC_k3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_annih_kpipipim_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_k3
  config.cross_section   = :default
end

exMC_kstar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_annih_kstar_pip_pim"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_kstar
  config.cross_section   = :default
end

# =====================================================================
# Event selection (BOSS)
# =====================================================================

# ---------------------------------------------------------------------
# mode A : k = 1  (Lambda K+ pi+ pi- pi0)
# ---------------------------------------------------------------------
alg_k1 = Algorithm.new("JpsiAnnihKpipipimPi0")
alg_k1.set_header(["JpsiAnnihKpipipimPi0Alg/JpsiAnnihKpipipimPi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:lambda_vertex_quality,
            "tag-side Lambda secondary-vertex fit: among p pi- combinations passing the fit
             only the candidate with the best vertex chi2 (ndf = 1) is kept; the track above
             500 MeV/c is taken as the proton candidate and the track below 500 MeV/c as the
             pion candidate (implemented as momentum conditions on the prp / pim lists).")
      .note(:recoil_mass_window,
            "single-tag recoil mass against the tagged Lambda required in
             RM_Lambda in [1.071, 1.160] GeV/c2; stored and applied at ROOT level
             (no BOSS-DSL construct for a recoil-mass window).")
      .note(:annihilation_vertex,
            "the Lambda-bar p annihilation vertex is required to lie in the beam-pipe
             cooling-oil layer, R_xy in [2.9, 3.6] cm, reconstructed from the
             K+ pi+ pi- (+ pi0) final state.")
      .note(:oil_proton_momentum,
            "inferred momentum of the proton at rest in the cooling oil required
             |p_oil| < 0.05 GeV/c.")

sel_k1 = Selection.new
  .select_track {                       # charged tracks
      cos_theta   0.93                  # |cos(theta)| < 0.93
      nChrp       ">=2"                 # at least two positive tracks
      nChrn       ">=2"                 # at least two negative tracks
  }
  .select_photon {                      # photons
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025           # 25 MeV (barrel)
      energyThreshold_e 0.050           # 50 MeV (endcap)
      nGam       ">=2"                  # k = 1 : at least 2 photons
  }
  .pid(method: :probability) {          # PID, probability method
      prob_cut   0.001
      identify :proton, against: [:kaon, :pion]   # p+/p- vs K and pi
      identify :kaon,   against: [:pion]         # K+/K- vs pi
      identify :pion,   against: [:kaon]         # pi+/pi- vs K
  }
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" }  # p > 500 MeV/c  -> proton candidate
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.5" }  # p < 500 MeV/c  -> pion candidate
  .kalman_kinematic_fit([:gamma, :gamma]) {       # pi0 from gamma gamma (1C mass constraint)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0  ">=1"                       # k = 1 : at least one pi0
  }
  .secondary_vertex_fit([:prp, :pim]) {           # tag-side Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      invariant_mass_of(:prp, :pim).between(1.110, 1.121)  # M(p pi-) window
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0]) {        # nominal 4C fit + Lambda mass constraint
      nominal
      constrain_four_momentum           # 4-momentum conservation to CMS
      invariant_mass_of(:Lambda).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
  }

alg_k1.with_decay_card(decay_card_k1).apply(sel_k1)

# ---------------------------------------------------------------------
# mode B : k = 2  (Lambda K+ pi+ pi- 2 pi0)
# ---------------------------------------------------------------------
alg_k2 = Algorithm.new("JpsiAnnihKpipipim2Pi0")
alg_k2.set_header(["JpsiAnnihKpipipim2Pi0Alg/JpsiAnnihKpipipim2Pi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:lambda_vertex_quality,
            "tag-side Lambda secondary-vertex fit: best vertex chi2 (ndf = 1) retained;
             p > 500 MeV/c taken as proton candidate, p < 500 MeV/c as pion candidate.")
      .note(:recoil_mass_window,
            "single-tag recoil mass against the tagged Lambda required in
             RM_Lambda in [1.071, 1.160] GeV/c2 (applied in ROOT).")
      .note(:annihilation_vertex,
            "annihilation vertex required in the cooling-oil layer, R_xy in [2.9, 3.6] cm.")
      .note(:oil_proton_momentum,
            "inferred oil-proton momentum required |p_oil| < 0.05 GeV/c.")

sel_k2 = Selection.new
  .select_track {
      cos_theta   0.93
      nChrp       ">=2"
      nChrn       ">=2"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam       ">=4"                  # k = 2 : at least 4 photons
  }
  .pid(method: :probability) {
      prob_cut   0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
  }
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" }
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.5" }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0  ">=2"                       # k = 2 : at least two pi0
  }
  .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      invariant_mass_of(:prp, :pim).between(1.110, 1.121)
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0, :pi0]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:Lambda).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
  }

alg_k2.with_decay_card(decay_card_k2).apply(sel_k2)

# ---------------------------------------------------------------------
# mode C : k = 3  (Lambda K+ pi+ pi- 3 pi0)
# ---------------------------------------------------------------------
alg_k3 = Algorithm.new("JpsiAnnihKpipipim3Pi0")
alg_k3.set_header(["JpsiAnnihKpipipim3Pi0Alg/JpsiAnnihKpipipim3Pi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:lambda_vertex_quality,
            "tag-side Lambda secondary-vertex fit: best vertex chi2 (ndf = 1) retained;
             500 MeV/c momentum assignment for the proton / pion candidates.")
      .note(:recoil_mass_window,
            "single-tag recoil mass against the tagged Lambda required in
             RM_Lambda in [1.071, 1.160] GeV/c2 (applied in ROOT).")
      .note(:annihilation_vertex,
            "annihilation vertex required in the cooling-oil layer, R_xy in [2.9, 3.6] cm.")
      .note(:oil_proton_momentum,
            "inferred oil-proton momentum required |p_oil| < 0.05 GeV/c.")

sel_k3 = Selection.new
  .select_track {
      cos_theta   0.93
      nChrp       ">=2"
      nChrn       ">=2"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam       ">=6"                  # k = 3 : at least 6 photons
  }
  .pid(method: :probability) {
      prob_cut   0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
  }
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" }
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.5" }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0  ">=3"                       # k = 3 : at least three pi0
  }
  .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      invariant_mass_of(:prp, :pim).between(1.110, 1.121)
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0, :pi0, :pi0]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:Lambda).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
  }

alg_k3.with_decay_card(decay_card_k3).apply(sel_k3)

# ---------------------------------------------------------------------
# mode D : K*(892)+ sub-process  (Lambda K*+ pi+ pi-, K*+ -> K+ pi0)
# ---------------------------------------------------------------------
alg_kstar = Algorithm.new("JpsiAnnihKstarPipPim")
alg_kstar.set_header(["JpsiAnnihKstarPipPimAlg/JpsiAnnihKstarPipPim.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:lambda_vertex_quality,
              "tag-side Lambda secondary-vertex fit: best vertex chi2 (ndf = 1) retained;
               500 MeV/c momentum assignment for the proton / pion candidates.")
        .note(:recoil_mass_window,
              "single-tag recoil mass against the tagged Lambda required in
               RM_Lambda in [1.071, 1.160] GeV/c2 (applied in ROOT).")
        .note(:annihilation_vertex,
              "annihilation vertex required in the cooling-oil layer, R_xy in [2.9, 3.6] cm.")
        .note(:oil_proton_momentum,
              "inferred oil-proton momentum required |p_oil| < 0.05 GeV/c.")

sel_kstar = Selection.new
  .select_track {
      cos_theta   0.93
      nChrp       ">=2"
      nChrn       ">=2"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam       ">=2"                  # K*+ pi+ pi- : at least 2 photons (K*+ -> K+ pi0)
  }
  .pid(method: :probability) {
      prob_cut   0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
  }
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" }
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.5" }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # pi0 of the K*+ -> K+ pi0 decay
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0  ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      invariant_mass_of(:prp, :pim).between(1.110, 1.121)
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:Lambda).constrain_to_nominal_mass_of(:Lambda)
      invariant_mass_of(:kp, :pi0).constrain_to_nominal_mass_of(:k_starp)  # K*(892)+ -> K+ pi0
      chi2_cut 200
  }

alg_kstar.with_decay_card(decay_card_kstar).apply(sel_kstar)

# =====================================================================
# Execute on real data, inclusive MC and the four exclusive-MC samples
# =====================================================================
root_files_k1    = alg_k1.execute_on([jpsi_data, jpsi_incMC, exMC_k1])
root_files_k2    = alg_k2.execute_on([jpsi_data, jpsi_incMC, exMC_k2])
root_files_k3    = alg_k3.execute_on([jpsi_data, jpsi_incMC, exMC_k3])
root_files_kstar = alg_kstar.execute_on([jpsi_data, jpsi_incMC, exMC_kstar])