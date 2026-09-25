# arXiv:1304.3205v1 — Search for the lepton flavor violation process J/psi -> e mu (BESIII)
# Final state: J/psi -> e+- mu-+  (two back-to-back charged tracks, no extra EMC activity)

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 225.3e6 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Inclusive J/psi MC

# Decay card for the LFV signal process
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ mu-    PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_Jpsi_emu"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                    cos_theta       0.8     # |cos(theta)| < 0.8
                    Vz              5.0     # closest approach < 5 cm along the beam direction
                    Vr              1.0     # closest approach < 1 cm transverse to the beam
                    nChrp           "==1"   # two well-measured tracks with zero net charge
                    nChrn           "==1"
                    nNet            "==0"
                }
               # Electron/muon separation from the MDC dE/dx, the EMC deposited energy
               # and the MUC penetration depth. The two tracks form an e/mu pair.
               .pid(method: :probability) {
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                   treat_as_electron_if_energy_above: 0.6
                    nlp    "==1"      # one lepton of either charge
                    nlm    "==1"
               }
               # Photon veto: the signal has no extra EMC activity. Showers in the gap
               # between barrel and endcap are not considered; showers closer than
               # 20 degrees to a charged track are attributed to that track.
               .select_photon {
                    tdc_emc_start     0
                    tdc_emc_end       14      # EMC cluster timing requirement against noise
                    energyThreshold_b 0.015   # E > 15 MeV in the barrel (|cos(theta)| < 0.80)
                    energyThreshold_e 0.015   # E > 15 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
                    angle_to_track    20.0    # > 20 degrees from the extrapolated tracks
                    nGam              "==0"   # veto events with any good photon candidate
               }
               # 4C kinematic fit of the two leptons to the J/psi four-momentum, which
               # also provides the E_vis/sqrt(s) and |sum p|/sqrt(s) quantities used to
               # define the signal region.
               .kinematic_fit([:lp, :lm]) {
                    nominal
                    constrain_four_momentum
                    chi2_cut 200
               }

my_algorithm
  .note(:background_veto, "cosmic-ray background rejected by requiring the TOF
    difference between the two charged tracks to be less than 1.0 ns")
  .note(:background_veto, "acollinearity and acoplanarity angles between the two
    charged tracks required to be less than 0.9 and 1.4 degrees respectively, to
    suppress e+e- -> e+e-(gamma) and e+e- -> mu+mu-(gamma) backgrounds")
  .note(:pid_correction_method, "electron identification requires no associated hits in
    the MUC, 0.95 < E/p < 1.50 (E is the EMC energy, p the MDC momentum) and
    |chi2_dE/dx(e)| < 1.8; muon identification requires |cos(theta)| < 0.75 (barrel MUC
    coverage), E/p < 0.5, 0.1 < E_EMC < 0.3 GeV, MUC penetration depth > 40 cm with
    chi2 < 100 when more than three MUC layers are hit, and chi2_dE/dx(e) < -1.8")
  .note(:background_veto, "the signal region is defined by 0.93 <= E_vis/sqrt(s) <= 1.10
    and |sum p|/sqrt(s) <= 0.10 (each about two standard deviations from MC); the yields
    are extracted from the scatter plot of these two variables in the ROOT analysis")
  .note(:efficiency_curve, "the selection was optimized in a blind analysis using a
    sensitivity figure of merit based on the Feldman-Cousins upper limit, with signal
    and six background MC samples (J/psi -> e+e-, mu+mu-, pi+pi-, K+K-,
    e+e- -> e+e-(gamma) and e+e- -> mu+mu-(gamma))")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
