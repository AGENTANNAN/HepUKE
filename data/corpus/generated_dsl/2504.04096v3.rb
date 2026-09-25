### Dataset description ###
# Representative energy point of the 59-point scan (4.009-4.950 GeV) implemented here
data_4260 = DatasetManager.real_data.find("703_4260")      # real data at 4.260 GeV (XYZ-I/II & R-scan sample)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")  # corresponding inclusive MC

# --- Decay card: e+e- -> pi+ pi- h_c, h_c -> gamma eta_c, eta_c -> K+ K- pi+ pi- ---
decay_card_kkpipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- h_c          PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c          PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- pi+ pi-       PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: e+e- -> pi+ pi- h_c, h_c -> gamma eta_c, eta_c -> pi+ pi- eta, eta -> gamma gamma ---
decay_card_pipieta = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- h_c          PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c          PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- eta          PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: e+e- -> pi+ pi- h_c, h_c -> gamma eta_c, eta_c -> 2(pi+ pi-) pi0, pi0 -> gamma gamma ---
decay_card_2pipipi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- h_c          PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c          PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Create the three 200k-event exclusive MC samples for the implemented eta_c modes
exMC_kkpipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4260_etac_kkpipi"
  config.related_dataset = data_4260
  config.events         = 200000
  config.decay_card     = decay_card_kkpipi
  config.cross_section  = :default
end

exMC_pipieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4260_etac_pipieta"
  config.related_dataset = data_4260
  config.events         = 200000
  config.decay_card     = decay_card_pipieta
  config.cross_section  = :default
end

exMC_2pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4260_etac_2pipipi0"
  config.related_dataset = data_4260
  config.events         = 200000
  config.decay_card     = decay_card_2pipipi0
  config.cross_section  = :default
end

### Event selection (BOSS) ###
# Mode I: eta_c -> K+ K- pi+ pi-   (final state pi+ pi- K+ K- gamma)

alg_name_kkpipi = "EtaCToKKPiPi"
alg_kkpipi = Algorithm.new(alg_name_kkpipi)
alg_kkpipi.set_header(["#{alg_name_kkpipi}Alg/#{alg_name_kkpipi}.h"])
          .set_constant({"ECMS" => [:double, 4.26]})
          .note(:scan_mode_deferred,
                "implemented for the representative 4.260 GeV point only; the full 59-point " \
                "4.009-4.950 GeV scan and the remaining 13 eta_c hadronic modes are deferred")

sel_kkpipi = Selection.new
  .select_track {            # charged track selection
    cos_theta 0.93           # |cos(theta)| < 0.93
    Vz        10.0           # |Vz| < 10 cm
    Vr        1.0            # Vr < 1 cm
    nChrp     ">=2"          # at least two positive tracks
    nChrn     ">=2"          # at least two negative tracks
    nNet      "==0"          # net charge zero
  }
  .select_photon {           # photon selection
    tdc_emc_start     0      # TDC window 0-14
    tdc_emc_end       14
    energyThreshold_b 0.025  # 25 MeV in the barrel
    energyThreshold_e 0.050  # 50 MeV in the endcap
    nGam              ">=1"  # at least one photon (from h_c -> gamma eta_c)
  }
  .pid(method: :probability) {   # probability PID
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ / K-
    identify :pion, against: [:kaon, :proton]   # pi+ / pi-
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  # Nominal 4C fit on pi+ pi- K+ K- gamma
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200            # loose cut; tight cut applied in ROOT
  }

alg_kkpipi.with_decay_card(decay_card_kkpipi).apply(sel_kkpipi)

# Mode II: eta_c -> pi+ pi- eta, eta -> gamma gamma   (final state pi+ pi- pi+ pi- gamma eta)

alg_name_pipieta = "EtaCToPiPiEta"
alg_pipieta = Algorithm.new(alg_name_pipieta)
alg_pipieta.set_header(["#{alg_name_pipieta}Alg/#{alg_name_pipieta}.h"])
           .set_constant({"ECMS" => [:double, 4.26]})
           .note(:scan_mode_deferred,
                 "implemented for the representative 4.260 GeV point only; the full 59-point " \
                 "4.009-4.950 GeV scan and the remaining 13 eta_c hadronic modes are deferred")

sel_pipieta = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"   # h_c -> gamma eta_c plus eta -> gamma gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  # Reconstruct eta from a photon pair (1-C mass constraint)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # Nominal 4C fit on pi+ pi- pi+ pi- gamma eta
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_pipieta.with_decay_card(decay_card_pipieta).apply(sel_pipieta)

# Mode III: eta_c -> 2(pi+ pi-) pi0, pi0 -> gamma gamma   (final state 3pi+ 3pi- gamma pi0)

alg_name_2pipipi0 = "EtaCTo2PiPiPi0"
alg_2pipipi0 = Algorithm.new(alg_name_2pipipi0)
alg_2pipipi0.set_header(["#{alg_name_2pipipi0}Alg/#{alg_name_2pipipi0}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})
            .note(:scan_mode_deferred,
                  "implemented for the representative 4.260 GeV point only; the full 59-point " \
                  "4.009-4.950 GeV scan and the remaining 13 eta_c hadronic modes are deferred")

sel_2pipipi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"          # at least three positive tracks
    nChrn     ">=3"          # at least three negative tracks
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"   # h_c -> gamma eta_c plus pi0 -> gamma gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==3"
    npim "==3"
  }
  # Reconstruct pi0 from a photon pair (1-C mass constraint)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 4C fit on 3pi+ 3pi- gamma pi0
  .kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim, :gamma, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_2pipipi0.with_decay_card(decay_card_2pipipi0).apply(sel_2pipipi0)

### Execute on datasets ###
root_files_kkpipi   = alg_kkpipi.execute_on([data_4260, incMC_4260, exMC_kkpipi])
root_files_pipieta  = alg_pipieta.execute_on([data_4260, incMC_4260, exMC_pipieta])
root_files_2pipipi0 = alg_2pipipi0.execute_on([data_4260, incMC_4260, exMC_2pipipi0])