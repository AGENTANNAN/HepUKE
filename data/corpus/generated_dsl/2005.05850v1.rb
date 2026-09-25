### Datasets ###
# Real data at the four scan points (BOSS 703): 4.467, 4.527, 4.575 and 4.600 GeV
data_4470 = DatasetManager.real_data.find("703_4470")
data_4530 = DatasetManager.real_data.find("703_4530")
data_4575 = DatasetManager.real_data.find("703_4575")
data_4600 = DatasetManager.real_data.find("703_4600")
data_points = [data_4470, data_4530, data_4575, data_4600]

# Corresponding inclusive MC samples
incMC_4470 = DatasetManager.inclusive_mc.find("703_4470")
incMC_4530 = DatasetManager.inclusive_mc.find("703_4530")
incMC_4575 = DatasetManager.inclusive_mc.find("703_4575")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_points = [incMC_4470, incMC_4530, incMC_4575, incMC_4600]

### Decay cards (EvtGen format) ###

# Mode I: e+e- -> Ds+ Ds1(2460)- ; Ds1(2460)- -> Ds*+ pi0, Ds*+ -> gamma Ds+, pi0 -> gamma gamma,
# Ds+ -> K+ K- pi+. Only the directly produced Ds+ is reconstructed; Ds1(2460)- is the recoil.
decay_card_DsDs1 = <<~DECAYCARD
  Alias Ds_recoil+ Ds+

  Decay psi(4260)
  1.0 Ds+ Ds1(2460)- PHSP;
  Enddecay

  Decay Ds1(2460)-
  1.0 Ds*+ pi0 PHSP;
  Enddecay

  Decay Ds*+
  1.0 gamma Ds_recoil+ PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K- pi+ PHSP;
  Enddecay

  Decay Ds_recoil+
  1.0 K+ K- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: e+e- -> Ds*+ Ds1(2460)- ; Ds*+ -> gamma Ds+, Ds+ -> K+ K- pi+,
# Ds1(2460)- -> Ds*+ pi0 -> gamma Ds+ pi0. The directly produced Ds*+ is reconstructed.
decay_card_DsstarDs1 = <<~DECAYCARD
  Alias Ds_star_recoil+ Ds*+
  Alias Ds_recoil+ Ds+

  Decay psi(4260)
  1.0 Ds*+ Ds1(2460)- PHSP;
  Enddecay

  Decay Ds*+
  1.0 gamma Ds+ PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K- pi+ PHSP;
  Enddecay

  Decay Ds1(2460)-
  1.0 Ds_star_recoil+ pi0 PHSP;
  Enddecay

  Decay Ds_star_recoil+
  1.0 gamma Ds_recoil+ PHSP;
  Enddecay

  Decay Ds_recoil+
  1.0 K+ K- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC ###
# Mode I: 200k events at each of the four energy points
exMC_mode1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_DsDs1"
  config.events        = 200000
  config.decay_card    = decay_card_DsDs1
  config.cross_section = :default
end

# Mode II: 100k events at 4.600 GeV only (small-scan points excluded)
exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DsstarDs1_4600"
  config.related_dataset = data_4600
  config.events          = 100000
  config.decay_card      = decay_card_DsstarDs1
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---------- Mode I: e+e- -> Ds+ Ds1(2460)- (no kinematic fit; recoil reconstruction) ----------
alg_name1 = "DsDs1"
alg_mode1 = Algorithm.new(alg_name1)
alg_mode1.set_header(["#{alg_name1}Alg/#{alg_name1}.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })
         .note(:ds_plus_resonance_windows,
               "Ds+ -> K+K-pi+ candidate selection uses intermediate-resonance windows " \
               "|M(K+K-)-m_phi|<15 MeV/c2 OR |M(K-pi+)-m_K*0|<84 MeV/c2 " \
               "(phi->K+K- and K*(892)0->K-pi+ sub-modes); not expressible in the current DSL")

sel_mode1 = Selection.new
  .select_track {                                    # charged-track selection
    cos_theta 0.93                                   # |cos(theta)| < 0.93
    Vz        10.0                                   # |Vz| < 10 cm
    Vr        1.0                                    # Vr < 1 cm
    nChrp     ">=2"                                  # at least 2 positive tracks
    nChrn     ">=1"                                  # at least 1 negative track
  }
  .pid(method: :probability) {                       # K/pi/p separation by probability method
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]        # K+ and K-
    identify :pion, against: [:kaon, :proton]        # pi+
    nkp  ">=1"                                       # at least one K+
    nkm  ">=1"                                       # at least one K-
    npip ">=1"                                       # at least one pi+
  }
  # Reconstruct the directly produced Ds+ from K+K-pi+ (best combination closest to the
  # nominal Ds+ mass) and infer Ds1(2460)- from the recoil four-momentum. No kinematic fit.
  .partial_rec([1]) {
    best_combination_by_mass :"Ds+", 1.96835
  }

# ---------- Mode II: e+e- -> Ds*+ Ds1(2460)- (2C mass-constrained fit) ----------
alg_name2 = "DsstarDs1"
alg_mode2 = Algorithm.new(alg_name2)
alg_mode2.set_header(["#{alg_name2}Alg/#{alg_name2}.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })

sel_mode2 = Selection.new
  .select_track {                                    # same tracking selection as Mode I
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=1"
  }
  .select_photon {                                   # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025                          # 25 MeV (barrel)
    energyThreshold_e 0.050                          # 50 MeV (endcap)
    angle_to_track    20.0                           # angle to nearest track > 20 deg
    nGam              ">=1"                          # at least one photon
  }
  .pid(method: :probability) {                       # same PID as Mode I
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  ">=1"
    nkm  ">=1"
    npip ">=1"
  }
  # 2C mass-constrained fit: constrain the K+K-pi+ system to the nominal Ds+ mass and the
  # gamma K+K-pi+ system to the nominal Ds*+ mass; chi2 < 10. Performed before the
  # recoil-mass study against the Ds*+ (which is done at ROOT level).
  .kinematic_fit([:gamma, :kp, :km, :pip]) {
    nominal
    invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:"Ds+")
    invariant_mass_of(:gamma, :kp, :km, :pip).constrain_to_nominal_mass_of(:"Ds*+")
    chi2_cut 10
  }

### Generate the algorithms ###
alg_mode1.with_decay_card(decay_card_DsDs1).apply(sel_mode1)
alg_mode2.with_decay_card(decay_card_DsstarDs1).apply(sel_mode2)

### Execute on datasets ###
root_files_mode1 = alg_mode1.execute_on(data_points + incMC_points + exMC_mode1)
root_files_mode2 = alg_mode2.execute_on([data_4600, incMC_4600, exMC_mode2])