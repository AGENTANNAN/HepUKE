# =====================================================================
# e+e- -> D_s*+ D_s*-   (c.m. 4.130 - 4.230 GeV, 8 points)
# Relative branching fraction of  D_s*+ -> D_s+ pi0  vs  D_s*+ -> D_s+ gamma
# Data: BOSS 705 @ 4.130 / 4.160 GeV ; BOSS 703 @ 4.180 - 4.230 GeV
# Three D_s decay modes (I, II, III); D_s* -> D_s gamma or D_s pi0 (pi0 -> gamma gamma)
# =====================================================================

### ------------------------------ Datasets ------------------------------ ###
# Real data points (sample name = [BOSS version]_[CMS energy in MeV])
data_4130 = DatasetManager.real_data.find("705_4130")   # 4.130 GeV
data_4160 = DatasetManager.real_data.find("705_4160")   # 4.160 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

# Matching inclusive MC samples
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

real_data_points = [data_4130, data_4160, data_4180, data_4190,
                    data_4200, data_4210, data_4220, data_4230]
incMC_points     = [incMC_4130, incMC_4160, incMC_4180, incMC_4190,
                    incMC_4200, incMC_4210, incMC_4220, incMC_4230]

### ---------------------------- Decay cards ---------------------------- ###
# Two cards per D_s decay mode: one for the D_s* -> D_s gamma hypothesis and
# one for the D_s* -> D_s pi0 hypothesis (needed to form the two signal templates).

# ---- Mode I : D_s+ -> K+ K- pi+ , D_s- -> K+ K- pi- ----
decay_card_modeI_gamma = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ gamma PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- gamma PHSP;
  Enddecay

  Decay D_s+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_modeI_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ pi0 PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- pi0 PHSP;
  Enddecay

  Decay D_s+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Mode II : D_s+ -> K+ K- pi+ , D_s- -> K_S0 K- ----
decay_card_modeII_gamma = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ gamma PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- gamma PHSP;
  Enddecay

  Decay D_s+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K_S0 K- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_modeII_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ pi0 PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- pi0 PHSP;
  Enddecay

  Decay D_s+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K_S0 K- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Mode III : D_s+ -> K_S0 K+ , D_s- -> K_S0 K- ----
decay_card_modeIII_gamma = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ gamma PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- gamma PHSP;
  Enddecay

  Decay D_s+
  1.0000 K_S0 K+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K_S0 K- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_modeIII_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 D_s+ pi0 PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- pi0 PHSP;
  Enddecay

  Decay D_s+
  1.0000 K_S0 K+ PHSP;
  Enddecay

  Decay D_s-
  1.0000 K_S0 K- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### --------------------- Exclusive signal MC (500k) --------------------- ###
# The same signal MC is run at every c.m. energy point (energy scan) ->
# create_exclusive_mc_for over the eight real-data points.
exMC_modeI_gamma = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeI_gamma"
  config.events        = 500_000
  config.decay_card    = decay_card_modeI_gamma
  config.cross_section = :default
end

exMC_modeI_pi0 = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeI_pi0"
  config.events        = 500_000
  config.decay_card    = decay_card_modeI_pi0
  config.cross_section = :default
end

exMC_modeII_gamma = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeII_gamma"
  config.events        = 500_000
  config.decay_card    = decay_card_modeII_gamma
  config.cross_section = :default
end

exMC_modeII_pi0 = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeII_pi0"
  config.events        = 500_000
  config.decay_card    = decay_card_modeII_pi0
  config.cross_section = :default
end

exMC_modeIII_gamma = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeIII_gamma"
  config.events        = 500_000
  config.decay_card    = decay_card_modeIII_gamma
  config.cross_section = :default
end

exMC_modeIII_pi0 = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
  config.sample_name   = "DsStar_modeIII_pi0"
  config.events        = 500_000
  config.decay_card    = decay_card_modeIII_pi0
  config.cross_section = :default
end

### ========================= Event selection (BOSS) ========================= ###

# ------------------------------- Mode I ------------------------------------ #
# D_s+ -> K+ K- pi+ ,  D_s- -> K+ K- pi-  (2 pos / 2 neg / 4 total)
sel_modeI = Selection.new
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93
    Vz        10.0          # |Vz| < 10 cm
    Vr        1.0           # Vr < 1 cm
    nChrp     "==2"         # exactly 2 positive tracks
    nChrn     "==2"         # exactly 2 negative tracks
    nTot      "==4"         # exactly 4 charged tracks in total
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025 # photon energy >= 25 MeV (barrel)
    energyThreshold_e 0.050 # photon energy >= 50 MeV (endcap)
  }
  .pid(method: :probability) {
    prob_cut 0.001                      # PID probability level 0.001
    identify :kaon, against: [:pion, :electron]   # K+/- separated against pi and e
    identify :pion, against: [:kaon, :electron]   # pi+/- separated against K and e
  }
  # 2C kinematic fit on the D_s daughters: both invariant masses constrained
  # to the nominal D_s mass.
  .kinematic_fit([:kp, :km, :pip, :kp, :km, :pim]) {
    nominal
    invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:D_s)  # D_s+ -> K+ K- pi+
    invariant_mass_of(:kp, :km, :pim).constrain_to_nominal_mass_of(:D_s)  # D_s- -> K+ K- pi-
    chi2_cut 200
  }

# ------------------------------- Mode II ----------------------------------- #
# D_s+ -> K+ K- pi+ ,  D_s- -> K_S0 K-   (2 pos / 3 neg / 5 total)
sel_modeII = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"         # exactly 2 positive tracks
    nChrn     "==3"         # exactly 3 negative tracks
    nTot      "==5"         # exactly 5 charged tracks in total
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :electron]
    identify :pion, against: [:kaon, :electron]
  }
  # one K_S0 candidate: secondary-vertex fit of the pi+ pi- pair
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 2C kinematic fit: both D_s masses constrained
  .kinematic_fit([:kp, :km, :pip, :K_S0, :km]) {
    nominal
    invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:D_s)  # D_s+ -> K+ K- pi+
    invariant_mass_of(:K_S0, :km).constrain_to_nominal_mass_of(:D_s)      # D_s- -> K_S0 K-
    chi2_cut 200
  }

# ------------------------------- Mode III ---------------------------------- #
# D_s+ -> K_S0 K+ ,  D_s- -> K_S0 K-   (2 pos / 2 neg / 4 total)
sel_modeIII = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"         # exactly 2 positive tracks
    nChrn     "==2"         # exactly 2 negative tracks
    nTot      "==4"         # exactly 4 charged tracks in total
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :electron]
    identify :pion, against: [:kaon, :electron]
  }
  # two K_S0 candidates
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 2C kinematic fit: both D_s masses constrained
  .kinematic_fit([:K_S0, :kp, :K_S0, :km]) {
    nominal
    invariant_mass_of(:K_S0, :kp).constrain_to_nominal_mass_of(:D_s)   # D_s+ -> K_S0 K+
    invariant_mass_of(:K_S0, :km).constrain_to_nominal_mass_of(:D_s)   # D_s- -> K_S0 K-
    chi2_cut 200
  }

### --------------------------- Build algorithms ---------------------------- ###
# Mode I (gamma / pi0 hypotheses share the same reconstruction)
alg_modeI_gamma = Algorithm.new("DsStarModeIGamma")
alg_modeI_gamma.set_header(["DsStarModeIGammaAlg/DsStarModeIGamma.h"])
               .set_constant({"ECMS" => [:double, 4.180]})
               .note(:ecms_per_energy_point,
                     "data span eight c.m. energies 4.130-4.230 GeV; the ECMS constant is set to the representative 4.180 GeV value, BOSS re-runs per energy point")
               .note(:miss_mass_separation,
                     "D_s* -> D_s gamma vs D_s* -> D_s pi0 discrimination performed in ROOT from M_miss^2 (recoil against the two reconstructed D_s); not a BOSS cut")

alg_modeI_pi0 = Algorithm.new("DsStarModeIPi0")
alg_modeI_pi0.set_header(["DsStarModeIPi0Alg/DsStarModeIPi0.h"])
             .set_constant({"ECMS" => [:double, 4.180]})
             .note(:ecms_per_energy_point,
                   "data span eight c.m. energies 4.130-4.230 GeV; ECMS constant set to the representative 4.180 GeV")
             .note(:miss_mass_separation,
                   "D_s* -> D_s gamma vs D_s* -> D_s pi0 discrimination performed in ROOT from M_miss^2")

# Mode II
alg_modeII_gamma = Algorithm.new("DsStarModeIIGamma")
alg_modeII_gamma.set_header(["DsStarModeIIGammaAlg/DsStarModeIIGamma.h"])
                .set_constant({"ECMS" => [:double, 4.180]})
                .note(:ecms_per_energy_point,
                      "data span eight c.m. energies 4.130-4.230 GeV; ECMS constant set to the representative 4.180 GeV")
                .note(:ks0_quality,
                      "K_S0 candidates required to have invariant mass inside 0.485-0.510 GeV/c^2 and flight significance (L/sigma) above 2; these windows are applied on the secondary-vertex output and are not expressible in the DSL")
                .note(:miss_mass_separation,
                      "D_s* decay-channel discrimination performed in ROOT from M_miss^2")

alg_modeII_pi0 = Algorithm.new("DsStarModeIIPi0")
alg_modeII_pi0.set_header(["DsStarModeIIPi0Alg/DsStarModeIIPi0.h"])
              .set_constant({"ECMS" => [:double, 4.180]})
              .note(:ecms_per_energy_point,
                    "data span eight c.m. energies 4.130-4.230 GeV; ECMS constant set to the representative 4.180 GeV")
              .note(:ks0_quality,
                    "K_S0 candidates required to have invariant mass inside 0.485-0.510 GeV/c^2 and flight significance (L/sigma) above 2; applied on the secondary-vertex output")
              .note(:miss_mass_separation,
                    "D_s* decay-channel discrimination performed in ROOT from M_miss^2")

# Mode III
alg_modeIII_gamma = Algorithm.new("DsStarModeIIIGamma")
alg_modeIII_gamma.set_header(["DsStarModeIIIGammaAlg/DsStarModeIIIGamma.h"])
                 .set_constant({"ECMS" => [:double, 4.180]})
                 .note(:ecms_per_energy_point,
                       "data span eight c.m. energies 4.130-4.230 GeV; ECMS constant set to the representative 4.180 GeV")
                 .note(:ks0_quality,
                       "K_S0 candidates required to have invariant mass inside 0.485-0.510 GeV/c^2 and flight significance (L/sigma) above 2; applied on the secondary-vertex output")
                 .note(:ks0_pairing,
                       "with two K_S0 candidates per event the pi+pi- pairing ambiguity is resolved by choosing the combination whose two vertices are closest (best common-vertex quality)")
                 .note(:miss_mass_separation,
                       "D_s* decay-channel discrimination performed in ROOT from M_miss^2")

alg_modeIII_pi0 = Algorithm.new("DsStarModeIIIPi0")
alg_modeIII_pi0.set_header(["DsStarModeIIIPi0Alg/DsStarModeIIIPi0.h"])
               .set_constant({"ECMS" => [:double, 4.180]})
               .note(:ecms_per_energy_point,
                     "data span eight c.m. energies 4.130-4.230 GeV; ECMS constant set to the representative 4.180 GeV")
               .note(:ks0_quality,
                     "K_S0 candidates required to have invariant mass inside 0.485-0.510 GeV/c^2 and flight significance (L/sigma) above 2; applied on the secondary-vertex output")
               .note(:ks0_pairing,
                     "two K_S0 candidates: pi+pi- pairing ambiguity resolved by the combination with the closest vertices")
               .note(:miss_mass_separation,
                     "D_s* decay-channel discrimination performed in ROOT from M_miss^2")

# Attach decay cards and render the selection into the BOSS C++ algorithm
alg_modeI_gamma.with_decay_card(decay_card_modeI_gamma).apply(sel_modeI)
alg_modeI_pi0.with_decay_card(decay_card_modeI_pi0).apply(sel_modeI.dup)

alg_modeII_gamma.with_decay_card(decay_card_modeII_gamma).apply(sel_modeII)
alg_modeII_pi0.with_decay_card(decay_card_modeII_pi0).apply(sel_modeII.dup)

alg_modeIII_gamma.with_decay_card(decay_card_modeIII_gamma).apply(sel_modeIII)
alg_modeIII_pi0.with_decay_card(decay_card_modeIII_pi0).apply(sel_modeIII.dup)

# Cross-feed / simultaneous fit are ROOT-level operations
alg_modeI_gamma.note(:cross_feed_transfer,
                     "the D_s* -> D_s pi0 / D_s gamma ratio is extracted by a simultaneous fit over the eight energy points using MC-derived cross-feed transfer factors (ROOT level)")

### ------------------------------ Execution -------------------------------- ###
root_files_modeI_gamma   = alg_modeI_gamma.execute_on(real_data_points + incMC_points + exMC_modeI_gamma)
root_files_modeI_pi0     = alg_modeI_pi0.execute_on(real_data_points + incMC_points + exMC_modeI_pi0)

root_files_modeII_gamma  = alg_modeII_gamma.execute_on(real_data_points + incMC_points + exMC_modeII_gamma)
root_files_modeII_pi0    = alg_modeII_pi0.execute_on(real_data_points + incMC_points + exMC_modeII_pi0)

root_files_modeIII_gamma = alg_modeIII_gamma.execute_on(real_data_points + incMC_points + exMC_modeIII_gamma)
root_files_modeIII_pi0   = alg_modeIII_pi0.execute_on(real_data_points + incMC_points + exMC_modeIII_pi0)