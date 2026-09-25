### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC

# Decay card for the signal channel: J/psi -> phi eta, phi -> K+K-, eta -> pi+pi-e+e-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 pi+ pi- e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the reference channel: J/psi -> phi eta, phi -> K+K-, eta -> gamma gamma
decay_card_ref = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Signal channel: J/psi -> phi(->K+K-) eta(->pi+pi-e+e-)
# ---------------------------------------------------------------------------
alg_signal_name = "PhiEtaToPipiEE"
alg_signal = Algorithm.new(alg_signal_name)
alg_signal.set_header(["#{alg_signal_name}Alg/#{alg_signal_name}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy = 3.097 GeV

sel_signal = Selection.new
sel_signal.select_track {                 # Charged track selection: six good tracks
    cos_theta 0.93                        # |cos(theta)| < 0.93
    Vz        10.0                        # |Vz| < 10 cm
    Vr        1.0                         # Vr < 1 cm
    nChrp     "==3"                       # K+, pi+, e+
    nChrn     "==3"                       # K-, pi-, e-
    nNet      "==0"                       # net charge zero
  }
  .pid(method: :probability) {            # Probability-based PID
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # tag the e+/e- (index_lp/index_lm)
    identify :kaon, against: [:pion]      # K+ and K- (charge-conjugation shorthand)
    identify :pion, against: [:kaon]      # pi+ and pi-
    nkp  "==1"                            # one K+
    nkm  "==1"                            # one K-
    npip "==1"                            # one pi+
    npim "==1"                            # one pi-
    nlp  "==1"                            # one lepton+ (electron)
    nlm  "==1"                            # one lepton- (electron)
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :lp, :lm]) {   # 4C fit over K+K-pi+pi-e+e-
    nominal                               # nominal fit (corrected four-momenta saved)
    constrain_four_momentum               # constrain total four-momentum to CMS energy
    chi2_cut 200                          # loose BOSS chi2 cut; tight cuts (20/40) applied in ROOT
  }

alg_signal.with_decay_card(decay_card_signal).apply(sel_signal)

# ---------------------------------------------------------------------------
# Reference channel: J/psi -> phi(->K+K-) eta(->gamma gamma)
# ---------------------------------------------------------------------------
alg_ref_name = "PhiEtaToGG"
alg_ref = Algorithm.new(alg_ref_name)
alg_ref.set_header(["#{alg_ref_name}Alg/#{alg_ref_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy = 3.097 GeV

sel_ref = Selection.new
sel_ref.select_track {                    # Two charged tracks with the same quality cuts
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"                       # K+
    nChrn     "==1"                       # K-
    nNet      "==0"                       # net charge zero
  }
  .select_photon {                        # Photon selection
    tdc_emc_start     0                   # EMC TDC start (unit: 700 ns)
    tdc_emc_end       14                  # EMC TDC end
    angle_to_track    20.0                # angle to nearest charged track > 20 degrees
    energyThreshold_b 0.025               # barrel energy threshold: 25 MeV
    energyThreshold_e 0.050               # endcap energy threshold: 50 MeV
    nGam              "==2"               # exactly two photons
  }
  .pid(method: :probability) {            # Kaon PID by probability
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K-
    nkp "==1"
    nkm "==1"
  }
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {   # 4C fit over K+K-gamma gamma
    nominal
    constrain_four_momentum
    chi2_cut 200                          # loose BOSS chi2 cut; tight cuts applied in ROOT
  }

# Upper photon-energy bound (1.4 GeV) has no dedicated DSL parameter — record it for ROOT/systematics
alg_ref.note(:efficiency_curve, "photon energy upper bound E < 1.4 GeV applied on the reference-channel
  eta -> gamma gamma selection (no dedicated BOSS photon-energy upper-cut parameter); enforced at ROOT level")

alg_ref.with_decay_card(decay_card_ref).apply(sel_ref)

### Execute on datasets ###
root_files_signal = alg_signal.execute_on([jpsi_data, jpsi_incMC])
root_files_ref    = alg_ref.execute_on([jpsi_data, jpsi_incMC])