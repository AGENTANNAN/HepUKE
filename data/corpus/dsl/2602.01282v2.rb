# BESIII: Λbar p → K+ π+ π- + k π0 (k=1,2,3) using J/ψ → Λ Λbar
# Λbar annihilates on protons at rest in the beam-pipe cooling oil layer.
# Λ is tagged via Λ → p π-; the signal side is reconstructed from remaining
# charged tracks (K+ π+ π-) and π0 candidates.

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### Decay cards ###
decay_card_k1 = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi-  pi0                     PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma gamma                           PHSP;
  Enddecay

  End
DECAYCARD

decay_card_k2 = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi-  pi0  pi0                PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma gamma                           PHSP;
  Enddecay

  End
DECAYCARD

decay_card_k3 = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi-  pi0  pi0  pi0           PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma gamma                           PHSP;
  Enddecay

  End
DECAYCARD

decay_card_kstar = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K*+  pi+  pi-                         PHSP;
  Enddecay

  Decay K*+
  1.0000  K+  pi0                               VSS;
  Enddecay

  Decay pi0
  1.0000  gamma gamma                           PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples ###
exMC_k1 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kpipi_pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_k1
  c.cross_section   = :default
end

exMC_k2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kpipi_2pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_k2
  c.cross_section   = :default
end

exMC_k3 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kpipi_3pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_k3
  c.cross_section   = :default
end

exMC_kstar = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kstarp_pipi"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_kstar
  c.cross_section   = :default
end

######################################################################
# Common note strings used for procedures not expressible in the DSL
######################################################################
NOTE_ANNIH_VERTEX = "Λbar p annihilation vertex is required to lie in the " \
  "beam-pipe cooling-oil layer: transverse distance of the K+ π+ π- vertex " \
  "from the beam axis R_xy ∈ [2.9, 3.6] cm."
NOTE_POIL = "Proton momentum inferred from p_p = p_final - p_Λbar (Λbar taken " \
  "opposite to the reconstructed Λ momentum); require |p_oil| < 0.05 GeV/c to " \
  "select at-rest hydrogen targets in the cooling oil."
NOTE_TAG_TRACKS = "On tag side: charged tracks with momentum > 500 MeV/c are " \
  "taken as proton candidates; those with momentum < 500 MeV/c as pion candidates."
NOTE_ST_REGION = "Single-tag Λ signal region defined on RM_Λ (recoil mass " \
  "against the reconstructed p π-): RM_Λ ∈ [1.071, 1.160] GeV/c^2."

######################################################################
# Algorithm 1: k = 1  (Λbar p → K+ π+ π- π0)
######################################################################
alg_k1 = Algorithm.new("LbarpKpipiPi0")
alg_k1.set_header(["LbarpKpipiPi0Alg/LbarpKpipiPi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_k1 = Selection.new
sel_k1.select_track {
        cos_theta 0.93
        nChrp    ">=2"
        nChrn    ">=2"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion]
        identify :pion,   against: [:kaon]
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=1"
      }
      .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
        chi2_cut 200
      }

alg_k1.note(:tag_track_momentum_split, NOTE_TAG_TRACKS)
      .note(:lambda_tag_mass_window, "Tagged Λ candidates required to have " \
            "M(p π-) ∈ [1.110, 1.121] GeV/c^2; multiple candidates resolved " \
            "by keeping the p π- combination with smallest vertex-fit χ² " \
            "(secondary vertex χ² < 200, ndf = 1).")
      .note(:st_signal_region, NOTE_ST_REGION)
      .note(:annihilation_vertex_region, NOTE_ANNIH_VERTEX)
      .note(:poil_cut, NOTE_POIL)
      .note(:best_pi0_selection, "Best π0 chosen by minimizing χ²_1C(π0) " \
            "over all π0 candidates in the k=1 topology.")
      .note(:kstar892_subprocess, "Sub-process Λbar p → K*(892)+ π+ π- with " \
            "K*(892)+ → K+ π0 studied in the M(K+ π0) spectrum for k=1; " \
            "K*(892)+ mass = 891.67 MeV/c^2, width = 51.4 MeV/c^2.")

alg_k1.with_decay_card(decay_card_k1).apply(sel_k1)
alg_k1.execute_on([jpsi_data, jpsi_incMC, exMC_k1, exMC_kstar])

######################################################################
# Algorithm 2: k = 2  (Λbar p → K+ π+ π- π0 π0)
######################################################################
alg_k2 = Algorithm.new("LbarpKpipi2Pi0")
alg_k2.set_header(["LbarpKpipi2Pi0Alg/LbarpKpipi2Pi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_k2 = Selection.new
sel_k2.select_track {
        cos_theta 0.93
        nChrp    ">=2"
        nChrn    ">=2"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion]
        identify :pion,   against: [:kaon]
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=2"
      }
      .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0, :pi0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
        chi2_cut 200
      }

alg_k2.note(:tag_track_momentum_split, NOTE_TAG_TRACKS)
      .note(:lambda_tag_mass_window, "Tagged Λ candidates required to have " \
            "M(p π-) ∈ [1.110, 1.121] GeV/c^2; multiple candidates resolved " \
            "by keeping the p π- combination with smallest vertex-fit χ².")
      .note(:st_signal_region, NOTE_ST_REGION)
      .note(:annihilation_vertex_region, NOTE_ANNIH_VERTEX)
      .note(:poil_cut, NOTE_POIL)
      .note(:best_pi0_selection, "Best π0 pair chosen by minimizing " \
            "χ²_1C(π0_1) + χ²_1C(π0_2) among all π0 combinations.")

alg_k2.with_decay_card(decay_card_k2).apply(sel_k2)
alg_k2.execute_on([jpsi_data, jpsi_incMC, exMC_k2])

######################################################################
# Algorithm 3: k = 3  (Λbar p → K+ π+ π- π0 π0 π0)
######################################################################
alg_k3 = Algorithm.new("LbarpKpipi3Pi0")
alg_k3.set_header(["LbarpKpipi3Pi0Alg/LbarpKpipi3Pi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_k3 = Selection.new
sel_k3.select_track {
        cos_theta 0.93
        nChrp    ">=2"
        nChrn    ">=2"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=6"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion]
        identify :pion,   against: [:kaon]
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=3"
      }
      .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:Lambda, :kp, :pip, :pim, :pi0, :pi0, :pi0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
        chi2_cut 200
      }

alg_k3.note(:tag_track_momentum_split, NOTE_TAG_TRACKS)
      .note(:lambda_tag_mass_window, "Tagged Λ candidates required to have " \
            "M(p π-) ∈ [1.110, 1.121] GeV/c^2; multiple candidates resolved " \
            "by keeping the p π- combination with smallest vertex-fit χ².")
      .note(:st_signal_region, NOTE_ST_REGION)
      .note(:annihilation_vertex_region, NOTE_ANNIH_VERTEX)
      .note(:poil_cut, NOTE_POIL)
      .note(:best_pi0_selection, "Best π0 triplet chosen by minimizing " \
            "χ²_1C(π0_1) + χ²_1C(π0_2) + χ²_1C(π0_3) among all combinations.")

alg_k3.with_decay_card(decay_card_k3).apply(sel_k3)
alg_k3.execute_on([jpsi_data, jpsi_incMC, exMC_k3])
