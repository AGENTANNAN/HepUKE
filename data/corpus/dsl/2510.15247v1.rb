### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data, (10.087 +- 0.044)x10^9 events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC sample

### Decay cards ###
# Signal process: J/psi -> gamma eta_c, eta_c -> p pbar
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c                                     PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p-                                      PHSP;
    Enddecay

    End
DECAYCARD

# Dominant exclusive BOSS-side background: J/psi -> p pbar pi0 (based on a preliminary
# amplitude analysis result; the other main background, J/psi -> p pbar gamma^F with an
# FSR photon, is already modeled within the inclusive MC sample via photos)
decay_card_bkg_pi0pp = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- pi0                                  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_signal      = DatasetManager.create_exclusive_mc { |c| c.sample_name = "jpsi_gamma_etac_ppbar"; c.related_dataset = jpsi_data; c.events = 200000; c.decay_card = decay_card_signal;     c.cross_section = :default }
exMC_bkg_pi0pp   = DatasetManager.create_exclusive_mc { |c| c.sample_name = "jpsi_ppbar_pi0";         c.related_dataset = jpsi_data; c.events = 200000; c.decay_card = decay_card_bkg_pi0pp; c.cross_section = :default }

### Event selection (BOSS) ###
# J/psi -> gamma eta_c, eta_c -> p pbar
alg = Algorithm.new("JpsiGammaEtacPPbar")
alg.set_header(["JpsiGammaEtacPPbarAlg/JpsiGammaEtacPPbar.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {                       # two charged tracks with zero net charge
                  nChrp     "==1"                     # exactly one positively charged track (proton)
                  nChrn     "==1"                     # exactly one negatively charged track (anti-proton)
                  nNet      "==0"                     # zero net charge
                  cos_theta 0.93                      # |cos(theta)| < 0.93
                  Vz        10.0                      # |Vz| < 10 cm
                  Vr        1.0                       # Vr < 1 cm
                }
                .select_photon {                      # good photon candidates from isolated EMC clusters
                  energyThreshold_b 0.025              # E > 25 MeV in barrel (|cos theta| < 0.8)
                  energyThreshold_e 0.050              # E > 50 MeV in endcap (0.86 < |cos theta| < 0.92)
                  tdc_emc_start      0                 # EMC timing within [0, 700] ns after event start time
                  tdc_emc_end        14
                  nGam               ">=1"             # at least one photon candidate required
                }
                .pid(method: :probability) {           # PID using dE/dx and TOF
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]   # likelihood for proton greater than for pion and kaon
                }
                .select_isolated_photon {              # reject clusters near extrapolated (anti-)proton positions
                  angle_to_prp_track 20.0               # cone angle of 20 deg around proton
                  angle_to_prm_track 30.0               # cone angle of 30 deg around anti-proton
                  nGam               ">=1"
                }
                .kinematic_fit([:gamma, :prp, :prm]) {  # 4C kinematic fit under J/psi -> gamma p pbar hypothesis
                  nominal
                  constrain_four_momentum              # energy-momentum conservation between initial and final states
                  chi2_cut 200                         # loose BOSS-level cut; combination with minimum chi2_4C chosen automatically
                }

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg_pi0pp])
