# Analysis: Search for J/psi -> Ds- rho+ and J/psi -> D0bar K*0bar
# Ds- -> phi e- nu_e, phi -> K+K-, rho+ -> pi+ pi0 (pi0 -> gamma gamma)
# D0bar -> K+ e- nu_e, K*0bar -> K- pi+

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for mode 1: J/psi -> Ds- rho+
decay_card_Dsrho = <<~DECAYCARD
    Decay J/psi
    1.0000  D_s- rho+                          SVS;
    Enddecay

    Decay D_s-
    1.0000  phi e- anti-nu_e                   PHOTOS  ISGW2;
    Enddecay

    Decay phi
    1.0000  K+ K-                              VSS;
    Enddecay

    Decay rho+
    1.0000  pi+ pi0                            VSS;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for mode 2: J/psi -> D0bar K*0bar
decay_card_D0Kst = <<~DECAYCARD
    Decay J/psi
    1.0000  anti-D0 anti-K*0                   SVS;
    Enddecay

    Decay anti-D0
    1.0000  K+ e- anti-nu_e                    PHOTOS  ISGW2;
    Enddecay

    Decay anti-K*0
    1.0000  K- pi+                             VSS;
    Enddecay

    End
DECAYCARD

exMC_Dsrho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_Dsm_rhop"
  config.related_dataset = jpsi_data
  config.events          = 600000
  config.decay_card      = decay_card_Dsrho
  config.cross_section   = :default
end

exMC_D0Kst = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_D0bar_Kstar0bar"
  config.related_dataset = jpsi_data
  config.events          = 600000
  config.decay_card      = decay_card_D0Kst
  config.cross_section   = :default
end

### Mode 1: J/psi -> Ds- rho+ ###
alg_Dsrho = Algorithm.new("JpsiDsRho")
alg_Dsrho.set_header(["JpsiDsRhoAlg/JpsiDsRho.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_Dsrho = Selection.new
sel_Dsrho.select_track {
             cos_theta 0.93         # |cos(theta)| < 0.93
             Vz        20.0         # |Vz| < 20 cm
             Vr        2.0          # |Vxy| < 2 cm
             nGood     "==4"        # 4 charged tracks
             nNet      "==0"        # net charge zero
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14   # EMC timing 0 <= t <= 700 ns
             angle_to_track    20.0 # > 20 degrees to any charged track
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam              ">=2"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
             identify :kaon, against: [:pion, :proton]
             identify :pion, against: [:kaon, :proton]
             nkp ">=1"
             nkm ">=1"
             npip ">=1"
             nem  ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0 ">=1"
           }
          .kinematic_fit([:kp, :km, :pip, :em, :pi0]) {
             nominal
             miss_track_of(:nu_e)
             constrain_four_momentum
             invariant_mass_of(:kp, :km).within(1.01, 1.03)   # phi mass window
             invariant_mass_of(:pip, :pi0).within(0.62, 0.95) # rho+ mass window
             chi2_cut 200
           }

alg_Dsrho.with_decay_card(decay_card_Dsrho).apply(sel_Dsrho)
alg_Dsrho.note(:electron_ep_cut,   "Electron candidate requires E/cP > 0.8 and |cos(theta)| < 0.8")
         .note(:missing_momentum,  "Pmiss > 0.1 GeV/c to reduce backgrounds without missing neutrino")
         .note(:umiss_cut,         "|Umiss = Emiss - c*Pmiss| < 0.05 GeV to reduce multi-pi0/multi-gamma backgrounds")

alg_Dsrho.execute_on([jpsi_data, jpsi_incMC, exMC_Dsrho])

### Mode 2: J/psi -> D0bar K*0bar ###
alg_D0Kst = Algorithm.new("JpsiD0Kstar")
alg_D0Kst.set_header(["JpsiD0KstarAlg/JpsiD0Kstar.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_D0Kst = Selection.new
sel_D0Kst.select_track {
             cos_theta 0.93
             Vz        20.0
             Vr        2.0
             nGood     "==4"
             nNet      "==0"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    20.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
             identify :kaon, against: [:pion, :proton]
             identify :pion, against: [:kaon, :proton]
             nkp ">=1"
             nkm ">=1"
             npip ">=1"
             nem  ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 20    # any pi0 with chi2 < 20 is vetoed later
           }
          .kinematic_fit([:kp, :km, :pip, :em]) {
             nominal
             miss_track_of(:nu_e)
             constrain_four_momentum
             invariant_mass_of(:km, :pip).within(0.82, 0.98)  # K*0bar mass window
             chi2_cut 200
           }

alg_D0Kst.with_decay_card(decay_card_D0Kst).apply(sel_D0Kst)
alg_D0Kst.note(:pi0_veto,          "Veto events with a pi0 candidate having chi2_1C(gamma gamma) < 20 to suppress backgrounds with pi0s")
         .note(:electron_ep_cut,   "Electron candidate requires E/cP > 0.8 and |cos(theta)| < 0.8")
         .note(:missing_momentum,  "Pmiss > 0.1 GeV/c")
         .note(:umiss_cut,         "|Umiss = Emiss - c*Pmiss| < 0.02 GeV")

alg_D0Kst.execute_on([jpsi_data, jpsi_incMC, exMC_D0Kst])
