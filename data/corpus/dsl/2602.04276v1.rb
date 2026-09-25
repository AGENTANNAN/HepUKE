# BESIII: Λbar p → light mesons (six channels)
# Λbar comes from J/ψ → Λ Λbar and annihilates on protons at rest in the
# cooling oil of the beam pipe. Λ is tagged via Λ → p π-.
# Signal reactions:
#   Λbar p → K_S0 π+
#   Λbar p → K+ π+ π-
#   Λbar p → K_S0 2π+ π-
#   Λbar p → K+ 2π+ 2π-
#   Λbar p → K+ 3π+ 3π-
#   Λbar p → 2K+ K-

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

######################################################################
# Common shared notes and helpers
######################################################################
NOTE_LAM_TAG = "Tagged Λ from Λ → p π-: best p π- combination chosen by " \
  "smallest |M(p π-) - m_Λ|; Λ signal region |M(p π-) - m_Λ| < 0.003 GeV/c^2 " \
  "(sideband 0.006 < |M(p π-) - m_Λ| < 0.009 GeV/c^2)."
NOTE_KS_TAG = "K_S0 candidate from K_S0 → π+ π-: best combination chosen by " \
  "smallest |M(π+ π-) - m_{K_S0}|; K_S0 signal region " \
  "|M(π+ π-) - m_{K_S0}| < 0.015 GeV/c^2."
NOTE_LBAR_RECOIL = "J/ψ → Λ Λbar selection: recoil mass against tagged Λ " \
  "required to satisfy |M_recoil(Λ) - m_Λbar| < 0.020 GeV/c^2."
NOTE_RXY = "Annihilation vertex R_xy ∈ [3.0, 3.5] cm to select the cooling-oil layer."
NOTE_POIL = "|P(p_oil)| < 0.04 GeV/c to suppress annihilations on " \
  "bound nucleons in 197Au/9Be/12C (Fermi momentum vs at-rest hydrogen)."

######################################################################
# Decay cards for each signal channel
######################################################################
decay_card_KsPip = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K_S0  pi+                             PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                              PHSP;
  Enddecay

  End
DECAYCARD

decay_card_KpPipPim = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi-                          PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Ks2PipPim = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K_S0  pi+  pi+  pi-                   PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                              PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Kp2Pip2Pim = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi+  pi-  pi-                PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Kp3Pip3Pim = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  pi+  pi+  pi+  pi-  pi-  pi-      PHSP;
  Enddecay

  End
DECAYCARD

decay_card_2KpKm = <<~DECAYCARD
  Decay J/psi
  1.0000  Lambda0  anti-Lambda0                 PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                               HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  K+  K+  K-                            PHSP;
  Enddecay

  End
DECAYCARD

######################################################################
# Exclusive MC samples
######################################################################
exMC_KsPip = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_KS0_pip"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_KsPip
  c.cross_section   = :default
end
exMC_KpPipPim = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kp_pip_pim"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_KpPipPim
  c.cross_section   = :default
end
exMC_Ks2PipPim = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_KS0_2pip_pim"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_Ks2PipPim
  c.cross_section   = :default
end
exMC_Kp2Pip2Pim = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kp_2pip_2pim"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_Kp2Pip2Pim
  c.cross_section   = :default
end
exMC_Kp3Pip3Pim = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_Kp_3pip_3pim"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_Kp3Pip3Pim
  c.cross_section   = :default
end
exMC_2KpKm = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Lbar_p_2Kp_Km"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_2KpKm
  c.cross_section   = :default
end

######################################################################
# Algorithm builder for each mode
######################################################################

# ------- Mode 1: Λbar p → K_S0 π+
alg1 = Algorithm.new("LbarpKSPip")
alg1.set_header(["LbarpKSPipAlg/LbarpKSPip.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel1 = Selection.new
sel1.select_track {
      cos_theta 0.93
      nChrp    ">=2"
      nChrn    ">=1"
      nTot     "==3"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :K_S0, :pip]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
      chi2_cut 200
    }
alg1.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:ks_tag,     NOTE_KS_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg1.with_decay_card(decay_card_KsPip).apply(sel1)
alg1.execute_on([jpsi_data, jpsi_incMC, exMC_KsPip])

# ------- Mode 2: Λbar p → K+ π+ π-
alg2 = Algorithm.new("LbarpKpPipPim")
alg2.set_header(["LbarpKpPipPimAlg/LbarpKpPipPim.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel2 = Selection.new
sel2.select_track {
      cos_theta 0.93
      nChrp    ">=2"
      nChrn    ">=2"
      nTot     "==4"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :kp, :pip, :pim]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
    }
alg2.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg2.with_decay_card(decay_card_KpPipPim).apply(sel2)
alg2.execute_on([jpsi_data, jpsi_incMC, exMC_KpPipPim])

# ------- Mode 3: Λbar p → K_S0 2π+ π-
alg3 = Algorithm.new("LbarpKS2PipPim")
alg3.set_header(["LbarpKS2PipPimAlg/LbarpKS2PipPim.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel3 = Selection.new
sel3.select_track {
      cos_theta 0.93
      nChrp    ">=3"
      nChrn    ">=2"
      nTot     "==5"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :K_S0, :pip, :pip, :pim]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
      chi2_cut 200
    }
alg3.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:ks_tag,     NOTE_KS_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg3.with_decay_card(decay_card_Ks2PipPim).apply(sel3)
alg3.execute_on([jpsi_data, jpsi_incMC, exMC_Ks2PipPim])

# ------- Mode 4: Λbar p → K+ 2π+ 2π-
alg4 = Algorithm.new("LbarpKp2Pip2Pim")
alg4.set_header(["LbarpKp2Pip2PimAlg/LbarpKp2Pip2Pim.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel4 = Selection.new
sel4.select_track {
      cos_theta 0.93
      nChrp    ">=3"
      nChrn    ">=3"
      nTot     "==6"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :kp, :pip, :pip, :pim, :pim]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
    }
alg4.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg4.with_decay_card(decay_card_Kp2Pip2Pim).apply(sel4)
alg4.execute_on([jpsi_data, jpsi_incMC, exMC_Kp2Pip2Pim])

# ------- Mode 5: Λbar p → K+ 3π+ 3π-
alg5 = Algorithm.new("LbarpKp3Pip3Pim")
alg5.set_header(["LbarpKp3Pip3PimAlg/LbarpKp3Pip3Pim.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel5 = Selection.new
sel5.select_track {
      cos_theta 0.93
      nChrp    ">=4"
      nChrn    ">=4"
      nTot     "==8"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :kp, :pip, :pip, :pip, :pim, :pim, :pim]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
    }
alg5.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg5.with_decay_card(decay_card_Kp3Pip3Pim).apply(sel5)
alg5.execute_on([jpsi_data, jpsi_incMC, exMC_Kp3Pip3Pim])

# ------- Mode 6: Λbar p → 2K+ K-
alg6 = Algorithm.new("Lbarp2KpKm")
alg6.set_header(["Lbarp2KpKmAlg/Lbarp2KpKm.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
sel6 = Selection.new
sel6.select_track {
      cos_theta 0.93
      nChrp    ">=3"
      nChrn    ">=2"
      nTot     "==5"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon,   against: [:pion]
      identify :pion,   against: [:kaon]
    }
    .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :kp, :kp, :km]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
      chi2_cut 200
    }
alg6.note(:lambda_tag, NOTE_LAM_TAG)
    .note(:lbar_recoil_window, NOTE_LBAR_RECOIL)
    .note(:annihilation_vertex_region, NOTE_RXY)
    .note(:poil_cut, NOTE_POIL)
alg6.with_decay_card(decay_card_2KpKm).apply(sel6)
alg6.execute_on([jpsi_data, jpsi_incMC, exMC_2KpKm])
