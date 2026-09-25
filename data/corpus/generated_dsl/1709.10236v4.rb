# =====================================================================
# e+e- -> Lambda Lambdabar : Lambda effective form factor determination
#   Modes  : I  (threshold, Lambda->p pi-, Lambdabar->pbar pi+, 2 pions only)
#            II (threshold, Lambdabar->nbar pi0, pi0->gamma gamma)
#            III(2.400 / 2.800 / 3.080 GeV, Lambda->p pi-, Lambdabar->pbar pi+)
# =====================================================================

### ---------------------------- Datasets ---------------------------- ###
# Real data and inclusive MC (BOSS 713 R-scan samples)
data_thr   = DatasetManager.real_data.find("713_2232")    # sqrt(s) = 2.2324 GeV (just above threshold)
incMC_thr  = DatasetManager.inclusive_mc.find("713_2232")
data_2400  = DatasetManager.real_data.find("713_2396")    # 2.400 GeV scan point
incMC_2400 = DatasetManager.inclusive_mc.find("713_2396")
data_2800  = DatasetManager.real_data.find("713_2800")    # 2.800 GeV
incMC_2800 = DatasetManager.inclusive_mc.find("713_2800")
data_3080  = DatasetManager.real_data.find("713_3080")    # 3.080 GeV
incMC_3080 = DatasetManager.inclusive_mc.find("713_3080")

### --------------------------- Decay cards -------------------------- ###
# Threshold mode I : e+e- -> Lambda Lambdabar, Lambda->p pi-, Lambdabar->pbar pi+
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Threshold mode II : e+e- -> Lambda Lambdabar, Lambdabar->nbar pi0, pi0->gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-n- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Higher energies : ConExc (ISR Born cross section) card,
# e+e- -> Lambda Lambdabar, Lambda->p pi-, Lambdabar->pbar pi+
# (Particle vpho is injected per energy point by the DSL for multi-energy scans)
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc 74110;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

### -------------------------- Exclusive MC -------------------------- ###
# Threshold mode I : 200k events
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_LambdaLambdabar_thr_modeI"
  config.related_dataset = data_thr
  config.events          = 200_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

# Threshold mode II : 200k events
exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_LambdaLambdabar_thr_modeII"
  config.related_dataset = data_thr
  config.events          = 200_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ConExc signal MC, 200k events at each of 2.400, 2.800 and 3.080 GeV
exMC_conexc = DatasetManager.create_exclusive_mc_for([data_2400, data_2800, data_3080]) do |config|
  config.sample_name   = "exmc_LambdaLambdabar_conexc"
  config.events        = 200_000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

### ------------------- Mode I (threshold) BOSS ---------------------- ###
alg_modeI = Algorithm.new("LamLambarModeI")
alg_modeI.set_header(["LamLambarModeIAlg/LamLambarModeI.h"])
         .set_constant({"ECMS" => [:double, 2.2324]})

sel_modeI = Selection.new
  .select_track {                 # exactly two opposite-charge tracks, net charge 0
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .pid(method: :probability) {    # identify both tracks as pions (pi/K/p separation)
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip     "==1"
    npim     "==1"
  }
  .remove(:pip) { condition "three_momentum_of(:pip) < 0.08 || three_momentum_of(:pip) > 0.11" }
  .remove(:pim) { condition "three_momentum_of(:pim) < 0.08 || three_momentum_of(:pim) > 0.11" }
  .partial_rec([4, 6]) {
    # reconstruct only the two pions (pi- from Lambda and pi+ from Lambdabar);
    # Lambda and Lambdabar remain untagged, no kinematic fit is performed
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_thr, incMC_thr, exMC_modeI])

### ------------------- Mode II (threshold) BOSS --------------------- ###
alg_modeII = Algorithm.new("LamLambarModeII")
alg_modeII.set_header(["LamLambarModeIIAlg/LamLambarModeII.h"])
          .set_constant({"ECMS" => [:double, 2.2324]})

sel_modeII = Selection.new
  .select_track {                 # at most one good charged track
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      "<=1"
  }
  .select_photon {                # >=3 photons, E>25 MeV barrel / 50 MeV endcap, no track within 10 deg
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit of the photon pair to the pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0     ">=1"
  }
  .partial_miss([1]) {
    # reconstruct only nbar pi0 (i.e. the Lambdabar); the recoil Lambda stays untagged
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_thr, incMC_thr, exMC_modeII])

### --------- Mode III (2.400 / 2.800 / 3.080 GeV) common selection --- ###
sel_modeIII = Selection.new
  .select_track {                 # two positive and two negative tracks, net charge 0
    cos_theta 0.93
    Vz        30.0
    Vr        10.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .pid(method: :probability) {    # one p+, one pbar, one pi+, one pi- (p/pi/K separation)
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    identify :pion,   against: [:proton, :kaon]
    nprp     "==1"
    nprm     "==1"
    npip     "==1"
    npim     "==1"
  }
  .secondary_vertex_fit([:prp, :pim]) {       # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # Lambdabar -> pbar pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .for_each(:Lambda) {                        # Lambda mass window |M(p pi) - M_Lambda| < 0.01 GeV/c^2
    where { "abs(mass - 1.1157) > 0.01" }
    remove
  }
  .for_each(:Lambda_bar) {                    # Lambda_bar mass window
    where { "abs(mass - 1.1157) > 0.01" }
    remove
  }
  .kinematic_fit([:Lambda, :Lambda_bar]) {    # 4C fit to the initial e+e- four-momentum
    nominal
    constrain_four_momentum
    chi2_cut 200                              # loose BOSS cut; tight/angle/M(Lambda Lambdabar) cuts done in ROOT
  }

# 2.400 GeV
alg_2400 = Algorithm.new("LamLambarModeIII2400")
alg_2400.set_header(["LamLambarModeIII2400Alg/LamLambarModeIII2400.h"])
        .set_constant({"ECMS" => [:double, 2.400]})
alg_2400.with_decay_card(decay_card_conexc).apply(sel_modeIII.dup)
alg_2400.execute_on([data_2400, incMC_2400, exMC_conexc[0]])

# 2.800 GeV
alg_2800 = Algorithm.new("LamLambarModeIII2800")
alg_2800.set_header(["LamLambarModeIII2800Alg/LamLambarModeIII2800.h"])
        .set_constant({"ECMS" => [:double, 2.800]})
alg_2800.with_decay_card(decay_card_conexc).apply(sel_modeIII.dup)
alg_2800.execute_on([data_2800, incMC_2800, exMC_conexc[1]])

# 3.080 GeV
alg_3080 = Algorithm.new("LamLambarModeIII3080")
alg_3080.set_header(["LamLambarModeIII3080Alg/LamLambarModeIII3080.h"])
        .set_constant({"ECMS" => [:double, 3.080]})
alg_3080.with_decay_card(decay_card_conexc).apply(sel_modeIII.dup)
alg_3080.execute_on([data_3080, incMC_3080, exMC_conexc[2]])