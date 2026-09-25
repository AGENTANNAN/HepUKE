### Dataset preparation ###
# psi(3686) data and inclusive MC
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for Channel I: psi(3686) -> gamma chi_cJ -> gamma anti-p- K*+ Lambda
# K*+ -> K+ pi0, pi0 -> gamma gamma, Lambda -> p+ pi-
# Includes chi_c0, chi_c1, chi_c2 with E1 radiative transition generators
decay_card_chi_cJ = <<~DECAYCARD
    Decay psi(2S)
    0.333 gamma chi_c0 P2GC0;
    0.333 gamma chi_c1 P2GC1;
    0.334 gamma chi_c2 P2GC2;
    Enddecay
    Decay chi_c0
    1.000 anti-p- K*+ Lambda PHSP;
    Enddecay
    Decay chi_c1
    1.000 anti-p- K*+ Lambda PHSP;
    Enddecay
    Decay chi_c2
    1.000 anti-p- K*+ Lambda PHSP;
    Enddecay
    Decay K*+
    1.000 K+ pi0 VSS;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay Lambda
    1.000 p+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card for Channel II: psi(3686) -> anti-p- K*+ Lambda (direct, no radiative photon)
decay_card_direct = <<~DECAYCARD
    Decay psi(2S)
    1.000 anti-p- K*+ Lambda PHSP;
    Enddecay
    Decay K*+
    1.000 K+ pi0 VSS;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay Lambda
    1.000 p+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_chi_cJ = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_chi_cJ_to_pbar_Kstar_Lambda"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_chi_cJ
  config.cross_section = :default
end

exMC_direct = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_to_pbar_Kstar_Lambda"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_direct
  config.cross_section = :default
end

### Channel I: psi(3686) -> gamma chi_cJ -> gamma anti-p- K*+ Lambda ###
# Three chi_cJ states share identical final states and selection criteria;
# a single Algorithm instance serves all three (chi_c0, chi_c1, chi_c2).

alg_chi_cJ = Algorithm.new("ChiCJToPbarKStarLambda")
alg_chi_cJ.set_header(["ChiCJToPbarKStarLambdaAlg/ChiCJToPbarKStarLambda.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .note(:pid_decay_track, "Only dE/dx information used for PID of pi- from Lambda decays; TOF unreachable for low-momentum Lambda daughters")
           .note(:background_veto, "Competing hypothesis veto in ROOT: CL(signal 5C) > CL(alt-5C no-radiative-gamma) AND CL(signal 5C) > CL(alt-4C no-pi0)")

sel_chi_cJ = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 5.0
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .remove([:kp <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  .select_isolated_photon {
    angle_to_prm_track 10.0
    angle_to_prp_track 10.0
    nGam ">=3"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 5C fit: gamma anti-p- K+ Lambda pi0 (pi0 mass constrained to nominal)
  .kinematic_fit([:gamma, :prm, :kp, :Lambda, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 70
  }
  # Competing 5C: anti-p- K+ Lambda pi0 (no radiative photon) -- stores chi2 for ROOT veto
  .kinematic_fit([:prm, :kp, :Lambda, :pi0]) {
    constrain_four_momentum
  }
  # Competing 4C: gamma anti-p- K+ Lambda (no pi0, all gammas as standalone) -- stores chi2 for ROOT veto
  .kinematic_fit([:gamma, :gamma, :gamma, :prm, :kp, :Lambda]) {
    constrain_four_momentum
  }

alg_chi_cJ.with_decay_card(decay_card_chi_cJ).apply(sel_chi_cJ)
alg_chi_cJ.execute_on([psip_data, psip_incMC, exMC_chi_cJ])

### Channel II: psi(3686) -> anti-p- K*+ Lambda (direct, no radiative photon) ###

alg_direct = Algorithm.new("PsiPToPbarKStarLambda")
alg_direct.set_header(["PsiPToPbarKStarLambdaAlg/PsiPToPbarKStarLambda.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .note(:pid_decay_track, "Only dE/dx information used for PID of pi- from Lambda decays; TOF unreachable for low-momentum Lambda daughters")
           .note(:background_veto, "Competing hypothesis veto in ROOT: CL(signal 5C) > CL(alt-5C with-radiative-gamma) AND CL(signal 5C) > CL(alt-4C no-pi0)")

sel_direct = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 5.0
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .remove([:kp <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  .select_isolated_photon {
    angle_to_prm_track 10.0
    angle_to_prp_track 10.0
    nGam ">=2"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 5C fit: anti-p- K+ Lambda pi0 (pi0 mass constrained to nominal)
  .kinematic_fit([:prm, :kp, :Lambda, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }
  # Competing 5C: gamma anti-p- K+ Lambda pi0 (with radiative photon) -- stores chi2 for ROOT veto
  .kinematic_fit([:gamma, :prm, :kp, :Lambda, :pi0]) {
    constrain_four_momentum
  }
  # Competing 4C: gamma anti-p- K+ Lambda (no pi0) -- stores chi2 for ROOT veto
  .kinematic_fit([:gamma, :prm, :kp, :Lambda]) {
    constrain_four_momentum
  }

alg_direct.with_decay_card(decay_card_direct).apply(sel_direct)
alg_direct.execute_on([psip_data, psip_incMC, exMC_direct])