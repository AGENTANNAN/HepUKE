### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC

# Decay card for the 4-pion channel: anti-n p -> 2pi+ 2pi-
# Anti-n source: J/psi -> p pi- anti-n
decay_card_4pi = <<~DECAYCARD
    Decay J/psi
    1.0 p+ pi- anti-n0 PHSP;
    Enddecay

    Decay anti-n0
    1.0 pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the 5-pion channel: anti-n p -> 2pi+ 2pi- pi0
decay_card_5pi = <<~DECAYCARD
    Decay J/psi
    1.0 p+ pi- anti-n0 PHSP;
    Enddecay

    Decay anti-n0
    1.0 pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the 6-pion channel: anti-n p -> 2pi+ 2pi- 2pi0
decay_card_6pi = <<~DECAYCARD
    Decay J/psi
    1.0 p+ pi- anti-n0 PHSP;
    Enddecay

    Decay anti-n0
    1.0 pi+ pi- pi+ pi- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for each of the three final states (500k events each)
exMC_4pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_antin_p_2pi2pim"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_4pi
  config.cross_section   = :default
end

exMC_5pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_antin_p_2pi2pimpi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_5pi
  config.cross_section   = :default
end

exMC_6pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_antin_p_2pi2pim2pi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_6pi
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------
# Mode I: anti-n p -> 2pi+ 2pi-  (no pi0, no photon requirement)
# ------------------------------------------------------------------
alg_name_4pi = "AntiNp4Pi"
alg_4pi = Algorithm.new(alg_name_4pi)
alg_4pi.set_header(["#{alg_name_4pi}Alg/#{alg_name_4pi}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:beam_pipe_target, "anti-n annihilates on the proton of the beam-pipe
         cooling oil; the target proton is not reconstructed and its Fermi motion /
         binding is not modelled in the exclusive MC, affecting the anti-n momentum
         resolution used for the five momentum-interval cross-section extraction")

sel_4pi = Selection.new
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        100.0   # |Vz| < 100 cm
    Vr        10.0    # Vr < 10 cm in the transverse plane
    nChrp     ">=3"   # at least 3 positive tracks
    nChrn     ">=3"   # at least 3 negative tracks
  }
  .select_photon {
    tdc_emc_start     0      # EMC TDC start
    tdc_emc_end       14     # EMC TDC end
    angle_to_track    10.0   # > 10 deg from any charged track
    energyThreshold_b 0.025  # > 25 MeV in the barrel
    energyThreshold_e 0.050  # > 50 MeV in the endcap
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / anti-p
    identify :pion,   against: [:kaon, :proton] # pi+ / pi-
    nprp ">=1"        # at least 1 proton
    npip ">=2"        # at least 2 pi+
    npim ">=2"        # at least 2 pi-
  }
  # 4C tag-side fit on p pi- with a missing anti-n
  .kinematic_fit([:prp, :pim]) {
    nominal
    miss_track_of :anti_n0     # anti-n escapes undetected
    constrain_four_momentum
    chi2_cut 200
  }

alg_4pi.with_decay_card(decay_card_4pi).apply(sel_4pi)

# ------------------------------------------------------------------
# Mode II: anti-n p -> 2pi+ 2pi- pi0  (>= 1 pi0, >= 2 photons)
# ------------------------------------------------------------------
alg_name_5pi = "AntiNp5Pi"
alg_5pi = Algorithm.new(alg_name_5pi)
alg_5pi.set_header(["#{alg_name_5pi}Alg/#{alg_name_5pi}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:beam_pipe_target, "anti-n annihilates on the proton of the beam-pipe
         cooling oil; the target proton is not reconstructed and its Fermi motion /
         binding is not modelled in the exclusive MC, affecting the anti-n momentum
         resolution used for the five momentum-interval cross-section extraction")

sel_5pi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     ">=3"
    nChrn     ">=3"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"    # at least 2 photons for the pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"
    npip ">=2"
    npim ">=2"
  }
  # Kalman fit to reconstruct pi0 -> gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    npi0 ">=1"                 # at least 1 pi0
  }
  # 4C tag-side fit on p pi- with a missing anti-n
  .kinematic_fit([:prp, :pim]) {
    nominal
    miss_track_of :anti_n0
    constrain_four_momentum
    chi2_cut 200
  }

alg_5pi.with_decay_card(decay_card_5pi).apply(sel_5pi)

# ------------------------------------------------------------------
# Mode III: anti-n p -> 2pi+ 2pi- 2pi0  (>= 2 pi0, >= 4 photons)
# ------------------------------------------------------------------
alg_name_6pi = "AntiNp6Pi"
alg_6pi = Algorithm.new(alg_name_6pi)
alg_6pi.set_header(["#{alg_name_6pi}Alg/#{alg_name_6pi}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:beam_pipe_target, "anti-n annihilates on the proton of the beam-pipe
         cooling oil; the target proton is not reconstructed and its Fermi motion /
         binding is not modelled in the exclusive MC, affecting the anti-n momentum
         resolution used for the five momentum-interval cross-section extraction")

sel_6pi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     ">=3"
    nChrn     ">=3"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"    # at least 4 photons for the two pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"
    npip ">=2"
    npim ">=2"
  }
  # Kalman fit to reconstruct pi0 -> gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    npi0 ">=2"                 # at least 2 pi0
  }
  # 4C tag-side fit on p pi- with a missing anti-n
  .kinematic_fit([:prp, :pim]) {
    nominal
    miss_track_of :anti_n0
    constrain_four_momentum
    chi2_cut 200
  }

alg_6pi.with_decay_card(decay_card_6pi).apply(sel_6pi)

### Execute on datasets ###
root_files_4pi = alg_4pi.execute_on([jpsi_data, jpsi_incMC, exMC_4pi])
root_files_5pi = alg_5pi.execute_on([jpsi_data, jpsi_incMC, exMC_5pi])
root_files_6pi = alg_6pi.execute_on([jpsi_data, jpsi_incMC, exMC_6pi])