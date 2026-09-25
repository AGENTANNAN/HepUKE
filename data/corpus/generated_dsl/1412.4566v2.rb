# ===========================================================================
# e+e- -> D*0 D0bar at sqrt(s) = 4.009 GeV
# The D*0 daughter (pi0 or gamma) is NOT reconstructed: the D0 D0bar system
# is fully reconstructed and the D0D0bar recoil mass separates the two
# D*0 decay channels.  Five tag modes (D0 decay, D0bar decay):
#   I   : K- pi+            , K+ pi-
#   II  : K- pi+            , K+ pi- pi0
#   III : K- pi+ pi0        , K+ pi-
#   IV  : K- pi+            , K+ pi- pi+ pi-
#   V   : K- pi+ pi+ pi-    , K+ pi-
# ===========================================================================

### Dataset description ###
data_4009  = DatasetManager.real_data.find("703_4009")      # 4.009 GeV, 482 pb^-1
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")   # matching inclusive MC

### Decay cards (EvtGen syntax) ###
# D*0 -> D0 pi0 (0.647) / D0 gamma (0.353) is common to all tag modes.

# Mode I : D0 -> K- pi+,  D0bar -> K+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    0.647 D0 pi0 PHSP;
    0.353 D0 gamma PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II : D0 -> K- pi+,  D0bar -> K+ pi- pi0
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    0.647 D0 pi0 PHSP;
    0.353 D0 gamma PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III : D0 -> K- pi+ pi0,  D0bar -> K+ pi-
decay_card_modeIII = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    0.647 D0 pi0 PHSP;
    0.353 D0 gamma PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV : D0 -> K- pi+,  D0bar -> K+ pi- pi+ pi-
decay_card_modeIV = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    0.647 D0 pi0 PHSP;
    0.353 D0 gamma PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode V : D0 -> K- pi+ pi+ pi-,  D0bar -> K+ pi-
decay_card_modeV = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    0.647 D0 pi0 PHSP;
    0.353 D0 gamma PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (200k events each) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dstar0D0bar_modeI"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dstar0D0bar_modeII"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dstar0D0bar_modeIII"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dstar0D0bar_modeIV"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_modeIV
  config.cross_section   = :default
end

exMC_modeV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dstar0D0bar_modeV"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_modeV
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# --------------------------------------------------------------------------
# Mode I : D0 -> K- pi+ , D0bar -> K+ pi-   (4 charged tracks, no pi0)
# --------------------------------------------------------------------------
alg_modeI_name = "Dstar0D0barModeI"
alg_modeI = Algorithm.new(alg_modeI_name)
alg_modeI.set_header(["#{alg_modeI_name}Alg/#{alg_modeI_name}.h"])
         .set_constant({"ECMS" => [:double, 4.009]})

sel_modeI = Selection.new
sel_modeI.select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nNet      "==0"     # zero net charge
    nChrp     "==2"     # 2 positive tracks
    nChrn     "==2"     # 2 negative tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001                            # 0.1% probability cut
    identify :kaon, against: [:pion]          # K+ and K- (pi/K separation)
    identify :pion, against: [:kaon]          # pi+ and pi- (pi/K separation)
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  # D0 and D0bar candidates constrained to the nominal D0 mass; best (smallest
  # total chi2) combination selected automatically
  .kinematic_fit([:km, :pip, :kp, :pim]) {
    nominal
    invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)   # D0  -> K- pi+
    invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)   # D0bar -> K+ pi-
    chi2_cut 30
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# --------------------------------------------------------------------------
# Mode II : D0 -> K- pi+ , D0bar -> K+ pi- pi0   (4 tracks + pi0)
# --------------------------------------------------------------------------
alg_modeII_name = "Dstar0D0barModeII"
alg_modeII = Algorithm.new(alg_modeII_name)
alg_modeII.set_header(["#{alg_modeII_name}Alg/#{alg_modeII_name}.h"])
          .set_constant({"ECMS" => [:double, 4.009]})

sel_modeII = Selection.new
sel_modeII.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     "==2"
    nChrn     "==2"
  }
  .select_photon {
    tdc_emc_start      0        # EMC time 0-700 ns
    tdc_emc_end        14
    angle_to_track     20.0     # at least 20 deg from any charged track
    energyThreshold_b  0.025    # > 25 MeV in the barrel
    energyThreshold_e  0.050    # > 50 MeV in the endcap
    nGam               ">=2"    # tag pi0 needs >= 2 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  # mass-constrained pi0 fit from the photon pair (|M(gg)-m(pi0)| < 15 MeV/c2)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # D0 and D0bar candidates constrained to the nominal D0 mass
  .kinematic_fit([:km, :pip, :kp, :pim, :pi0]) {
    nominal
    invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)         # D0    -> K- pi+
    invariant_mass_of(:kp, :pim, :pi0).constrain_to_nominal_mass_of(:D0)   # D0bar -> K+ pi- pi0
    chi2_cut 30
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# --------------------------------------------------------------------------
# Mode III : D0 -> K- pi+ pi0 , D0bar -> K+ pi-   (4 tracks + pi0)
# --------------------------------------------------------------------------
alg_modeIII_name = "Dstar0D0barModeIII"
alg_modeIII = Algorithm.new(alg_modeIII_name)
alg_modeIII.set_header(["#{alg_modeIII_name}Alg/#{alg_modeIII_name}.h"])
           .set_constant({"ECMS" => [:double, 4.009]})

sel_modeIII = Selection.new
sel_modeIII.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     "==2"
    nChrn     "==2"
  }
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     20.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:km, :pip, :pi0, :kp, :pim]) {
    nominal
    invariant_mass_of(:km, :pip, :pi0).constrain_to_nominal_mass_of(:D0)   # D0    -> K- pi+ pi0
    invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)         # D0bar -> K+ pi-
    chi2_cut 30
  }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)

# --------------------------------------------------------------------------
# Mode IV : D0 -> K- pi+ , D0bar -> K+ pi- pi+ pi-   (6 charged tracks)
# --------------------------------------------------------------------------
alg_modeIV_name = "Dstar0D0barModeIV"
alg_modeIV = Algorithm.new(alg_modeIV_name)
alg_modeIV.set_header(["#{alg_modeIV_name}Alg/#{alg_modeIV_name}.h"])
          .set_constant({"ECMS" => [:double, 4.009]})

sel_modeIV = Selection.new
sel_modeIV.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     "==3"     # 3 positive tracks
    nChrn     "==3"     # 3 negative tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:km, :pip, :kp, :pim, :pip, :pim]) {
    nominal
    invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)                  # D0    -> K- pi+
    invariant_mass_of(:kp, :pim, :pip, :pim).constrain_to_nominal_mass_of(:D0)      # D0bar -> K+ pi- pi+ pi-
    chi2_cut 30
  }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(sel_modeIV)

# --------------------------------------------------------------------------
# Mode V : D0 -> K- pi+ pi+ pi- , D0bar -> K+ pi-   (6 charged tracks)
# --------------------------------------------------------------------------
alg_modeV_name = "Dstar0D0barModeV"
alg_modeV = Algorithm.new(alg_modeV_name)
alg_modeV.set_header(["#{alg_modeV_name}Alg/#{alg_modeV_name}.h"])
         .set_constant({"ECMS" => [:double, 4.009]})

sel_modeV = Selection.new
sel_modeV.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     "==3"
    nChrn     "==3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:km, :pip, :pip, :pim, :kp, :pim]) {
    nominal
    invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)      # D0    -> K- pi+ pi+ pi-
    invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)                  # D0bar -> K+ pi-
    chi2_cut 30
  }

alg_modeV.with_decay_card(decay_card_modeV).apply(sel_modeV)

### Execute on data, inclusive MC and the corresponding signal MC ###
root_files_modeI   = alg_modeI.execute_on([data_4009, incMC_4009, exMC_modeI])
root_files_modeII  = alg_modeII.execute_on([data_4009, incMC_4009, exMC_modeII])
root_files_modeIII = alg_modeIII.execute_on([data_4009, incMC_4009, exMC_modeIII])
root_files_modeIV  = alg_modeIV.execute_on([data_4009, incMC_4009, exMC_modeIV])
root_files_modeV   = alg_modeV.execute_on([data_4009, incMC_4009, exMC_modeV])