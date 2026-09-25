# ============================================================================
# e+e- -> (D*Dbar*)^0 pi0  /  Z_c(4025)^0 -> (D*Dbar*)^0
# Partial reconstruction: one D (from D* -> D gamma) + bachelor pi0 are fully
# reconstructed; the recoiling Dbar* is inferred from the missing 4-momentum.
# Signal read from the RM(D pi0) spectrum.
# Data: 4.230 GeV (1092 pb-1) and 4.260 GeV (826 pb-1) + inclusive MC.
# Four D decay modes, each with its own 200k-event exclusive signal MC
# generated at BOTH energies.
# ============================================================================

### ---------------------------- Dataset preparation ----------------------- ###
data_4230  = DatasetManager.real_data.find("703_4230")   # 4.230 GeV real data
data_4260  = DatasetManager.real_data.find("703_4260")   # 4.260 GeV real data
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

### ---------------------------- Decay cards -------------------------------- ###
# Mode I: D0 -> K- pi+   (D*0 -> D0 gamma) ; both (D*Dbar*) sides present
decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Z_c(4025) pi0 PHSP;
    Enddecay

    Decay Z_c(4025)
    1.0000 D*0 anti-D*0 PHSP;
    Enddecay

    Decay D*0
    1.0000 D0 gamma PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 gamma PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: D0 -> K- pi+ pi0
decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Z_c(4025) pi0 PHSP;
    Enddecay

    Decay Z_c(4025)
    1.0000 D*0 anti-D*0 PHSP;
    Enddecay

    Decay D*0
    1.0000 D0 gamma PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 gamma PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: D0 -> K- pi+ pi+ pi-
decay_card_mode3 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Z_c(4025) pi0 PHSP;
    Enddecay

    Decay Z_c(4025)
    1.0000 D*0 anti-D*0 PHSP;
    Enddecay

    Decay D*0
    1.0000 D0 gamma PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 gamma PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi- pi+ PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV: D+ -> K- pi+ pi+  (charged mode, D*+ -> D+ gamma)
decay_card_mode4 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Z_c(4025) pi0 PHSP;
    Enddecay

    Decay Z_c(4025)
    1.0000 D*+ D*- PHSP;
    Enddecay

    Decay D*+
    1.0000 D+ gamma PHSP;
    Enddecay

    Decay D*-
    1.0000 D- gamma PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### --------- Exclusive signal MC: 200k events per mode at BOTH energies --- ###
# create_exclusive_mc_for runs the SAME signal MC over each energy point.
exMC_mode1 = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_D0toKPi"
  config.events        = 200000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_D0toKPiPi0"
  config.events        = 200000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_D0toKPiPiPi"
  config.events        = 200000
  config.decay_card    = decay_card_mode3
  config.cross_section = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_DplusToKPiPi"
  config.events        = 200000
  config.decay_card    = decay_card_mode4
  config.cross_section = :default
end

# Persist the generated MC configs (one file per energy point / mode)
[exMC_mode1, exMC_mode2, exMC_mode3, exMC_mode4].flatten.each do |m|
  m.save_to_config(format: :yaml, file_path: 'temp_for_test')
end

# Datasets common to every mode
common_data = [data_4230, data_4260, incMC_4230, incMC_4260]

### =========================== EVENT SELECTION ============================= ###

# ---------------------------------------------------------------------------
# MODE I : D0 -> K- pi+      (multiplicity 2+,2- ; PID: 1K+ 1K- 1pi+ 1pi-)
# ---------------------------------------------------------------------------
alg_name_mode1 = "Zc4025D0toKPi"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:efficiency_curve, "the CMS energy constant ECMS is set to 4.260 GeV; " \
               "the identical selection is also applied to the 4.230 GeV sample " \
               "(per-run ECMS handled by the framework)")
         .note(:background_veto, "events with M(D pi0) > 2.02 GeV/c**2 are rejected to " \
               "suppress D* -> D pi0 feed-down")
         .note(:pid_correction_method, "kaon/pion separation uses the probability method " \
               "(dE/dx and TOF) with probability cut 0.001")
         .note(:helix_correction, "the charged tracks forming the D are required to share a " \
               "common D vertex with chi2(D) < 100; the recoiling D* and Dbar* must share no " \
               "final-state particle and the D/Dbar pair is chosen by minimum summed chi2")
         .note(:efficiency_curve, "mass constraints: chi2(D) < 15 (< 20 when a pi0 is present) " \
               "and chi2(pi0) < 20; the tagged D and bachelor pi0 are selected by best mass " \
               "consistency to the nominal D0 mass (1.86484 GeV/c**2)")

sel_mode1 = Selection.new
  .select_track {
    cos_theta  0.93          # |cos(theta)| < 0.93
    Vz         10.0          # |Vz| < 10 cm
    Vr         1.0           # Vr < 1 cm
    nNet       "==0"         # net charge zero
    nChrp      "==2"         # 2+ multiplicity
    nChrn      "==2"         # 2- multiplicity
  }
  .select_photon {
    tdc_emc_start      0     # EMC timing window
    tdc_emc_end        14
    angle_to_track     10.0  # >= 10 deg from any charged track
    energyThreshold_b  0.025 # E > 25 MeV (barrel)
    energyThreshold_e  0.050 # E > 50 MeV (endcap)
    nGam               ">=2" # at least 2 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001                             # probability cut (dE/dx + TOF)
    identify :kaon, against: [:pion, :proton]  # K+ and K- (charge-conjugation shorthand)
    identify :pion, against: [:kaon, :proton]  # pi+ and pi-
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  # pi0 candidates built from photon pairs, mass window constrained to nominal pi0 mass
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)             # two-photon mass window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # constrain to m(pi0)
    chi2_cut 25                                                        # chi2 < 25
    npi0 ">=1"                                                         # at least one pi0
  }
  # Partial reconstruction: reconstruct the D (D0) and the bachelor pi0;
  # Dbar* inferred from the missing four-momentum -> recoil mass RM(D pi0)
  .partial_rec([5, 2]) {   # 5 = D0 (D*0 -> D0 gamma), 2 = bachelor pi0
    best_combination_by_mass :D0, 1.86484   # best consistency to nominal D0 mass
    require_recoil_mass 2.10, 2.22          # RM(D pi0) in [2.10, 2.22] GeV/c**2
  }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
root_files_mode1 = alg_mode1.execute_on(common_data + exMC_mode1)

# ---------------------------------------------------------------------------
# MODE II : D0 -> K- pi+ pi0  (multiplicity 2+,2- ; PID: 1K+ 1K- 1pi+ 1pi-)
# ---------------------------------------------------------------------------
alg_name_mode2 = "Zc4025D0toKPiPi0"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:background_veto, "events with M(D pi0) > 2.02 GeV/c**2 are rejected to " \
               "suppress D* -> D pi0 feed-down")
         .note(:pid_correction_method, "kaon/pion separation uses the probability method " \
               "(dE/dx and TOF) with probability cut 0.001")
         .note(:helix_correction, "track combinations forming the D are required to share a " \
               "common D vertex with chi2(D) < 100; D and Dbar share no final-state particle; " \
               "the pair is chosen by minimum summed chi2")
         .note(:efficiency_curve, "mass constraints: chi2(D) < 20 (pi0 present), chi2(pi0) < 20; " \
               "tagged D and bachelor pi0 chosen by best consistency to the nominal D0 mass " \
               "(1.86484 GeV/c**2)")
         .note(:efficiency_curve, "the bachelor pi0 is formed from the leftover photons after " \
               "the two D pi0's are assigned")

sel_mode2 = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nNet       "==0"
    nChrp      "==2"
    nChrn      "==2"
  }
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=4"   # >= 4 photons in the K- pi+ pi0 mode
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=3"          # at least three pi0 (two from the D's + bachelor)
  }
  .partial_rec([5, 2]) {   # 5 = D0 (with internal pi0), 2 = bachelor pi0
    best_combination_by_mass :D0, 1.86484
    require_recoil_mass 2.10, 2.22
  }

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
root_files_mode2 = alg_mode2.execute_on(common_data + exMC_mode2)

# ---------------------------------------------------------------------------
# MODE III : D0 -> K- pi+ pi+ pi-  (multiplicity 4+,4- ; PID: 1K+ 1K- 3pi+ 3pi-)
# ---------------------------------------------------------------------------
alg_name_mode3 = "Zc4025D0toKPiPiPi"
alg_mode3 = Algorithm.new(alg_name_mode3)
alg_mode3.set_header(["#{alg_name_mode3}Alg/#{alg_name_mode3}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:background_veto, "events with M(D pi0) > 2.02 GeV/c**2 are rejected to " \
               "suppress D* -> D pi0 feed-down")
         .note(:pid_correction_method, "kaon/pion separation uses the probability method " \
               "(dE/dx and TOF) with probability cut 0.001")
         .note(:helix_correction, "track combinations forming the D share a common D vertex " \
               "with chi2(D) < 100; D and Dbar share no final-state particle; the pair is " \
               "chosen by minimum summed chi2")
         .note(:efficiency_curve, "mass constraints: chi2(D) < 15, chi2(pi0) < 20; tagged D and " \
               "bachelor pi0 chosen by best consistency to the nominal D0 mass (1.86484 GeV/c**2)")

sel_mode3 = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nNet       "==0"
    nChrp      "==4"
    nChrn      "==4"
  }
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==3"
    npim "==3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .partial_rec([5, 2]) {   # 5 = D0, 2 = bachelor pi0
    best_combination_by_mass :D0, 1.86484
    require_recoil_mass 2.10, 2.22
  }

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)
root_files_mode3 = alg_mode3.execute_on(common_data + exMC_mode3)

# ---------------------------------------------------------------------------
# MODE IV : D+ -> K- pi+ pi+   (charged mode, D*+ -> D+ gamma)
#           multiplicity 3+,3- ; PID: 1K+ 1K- 2pi+ 2pi-
# ---------------------------------------------------------------------------
alg_name_mode4 = "Zc4025DplusToKPiPi"
alg_mode4 = Algorithm.new(alg_name_mode4)
alg_mode4.set_header(["#{alg_name_mode4}Alg/#{alg_name_mode4}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:background_veto, "events with M(D pi0) > 2.02 GeV/c**2 are rejected to " \
               "suppress D* -> D pi0 feed-down")
         .note(:pid_correction_method, "kaon/pion separation uses the probability method " \
               "(dE/dx and TOF) with probability cut 0.001")
         .note(:helix_correction, "track combinations forming the D share a common D vertex " \
               "with chi2(D) < 100; D and Dbar share no final-state particle; the pair is " \
               "chosen by minimum summed chi2")
         .note(:efficiency_curve, "mass constraints: chi2(D) < 15, chi2(pi0) < 20; tagged D and " \
               "bachelor pi0 chosen by best consistency to the nominal D+ mass (1.86962 GeV/c**2)")

sel_mode4 = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nNet       "==0"
    nChrp      "==3"
    nChrn      "==3"
  }
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .partial_rec([5, 2]) {   # 5 = D+, 2 = bachelor pi0
    best_combination_by_mass :Dplus, 1.86962  # nominal D+ mass
    require_recoil_mass 2.10, 2.22
  }

alg_mode4.with_decay_card(decay_card_mode4).apply(sel_mode4)
root_files_mode4 = alg_mode4.execute_on(common_data + exMC_mode4)