# arXiv:1303.7360v1 — Study of eta' -> pi+ pi- l+ l- via J/psi -> gamma eta' (BESIII, 2.253e8 J/psi)
# Three independent final states, each with its own Algorithm object (Rule T1):
#   (I)   J/psi -> gamma eta', eta' -> pi+ pi- e+ e-
#   (II)  J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu-
#   (III) J/psi -> gamma eta', eta' -> gamma pi+ pi-   (normalization channel)

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 2.253e8 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Inclusive J/psi MC

# ---------------- Decay cards ----------------
decay_card_ee = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'    PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- e+ e-    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'    PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- mu+ mu-    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_gammapipi = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'    PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

exMC_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_etap_ee"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_ee
    config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_etap_mumu"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_mumu
    config.cross_section   = :default
end

exMC_gammapipi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_etap_gammapipi"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_gammapipi
    config.cross_section   = :default
end

### Event selection (BOSS) ###

# ==================== (I) J/psi -> gamma eta', eta' -> pi+ pi- e+ e- ====================
alg_name_I = "EtapToPiPiEE"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta       0.93    # |cos(theta)| < 0.93
              Vz              10.0    # within +/-10 cm of the IP along the beam direction
              Vr              1.0     # within 1 cm in the plane perpendicular to the beam
              nChrp           "==2"   # exactly four good charged tracks
              nChrn           "==2"
              nNet            "==0"   # zero net charge
          }
         # Combinatorial PID: the particle type of each track (pion or electron) is
         # assigned by minimising the sum of the TOF/dE/dx chi^2 over the tracks;
         # combined with the 4C chi^2 this forms chi^2_{gamma pi+ pi- l+ l-}.
         .pid(method: :chi2_sum) {
              chi_min_cut 4
              identify :pion, :electron
         }
         .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14      # EMC timing requirement vs. the collision
              energyThreshold_b 0.025   # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
              energyThreshold_e 0.050   # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
              angle_to_track    15.0    # photon at least 15 degrees from any good track
              nGam              ">=1"   # at least one good photon
         }
         # 4C kinematic fit under the gamma pi+ pi- e+ e- hypothesis; the combination
         # (including the choice of the best photon) with the minimum chi^2 is kept.
         .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {
              nominal
              constrain_four_momentum
              chi2_cut 75
         }
         # Competing hypothesis: 4C fit under gamma 2(pi+ pi-), used to reject the
         # J/psi -> gamma 2(pi+ pi-) background. No chi2_cut, no nominal: only the
         # chi^2 value is stored and compared in the ROOT analysis.
         .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {
              use_track_index_from_nominal_kmfit
              constrain_four_momentum
         }

alg_modeI
  .note(:background_veto, "J/psi -> gamma 2(pi+ pi-) background rejected by requiring
    chi2_{gamma 2(pi+pi-)} > chi2_{gamma pi+ pi- e+ e-}; the competing 4C fit above stores
    the second chi^2 value for the ROOT-level comparison. The second hypothesis treats
    the lepton pair as a pion pair, which the BOSS participant list cannot re-label here")
  .note(:background_veto, "non-eta' background (J/psi -> pi+pi-pi0 and J/psi -> gamma pi+pi-pi0)
    estimated from the eta' sideband 0.88 < M(pi+pi-e+e-) < 0.90 GeV/c^2 or
    1.02 < M(pi+pi-e+e-) < 1.04 GeV/c^2; the eta' mass window
    |M(pi+pi-e+e-) - m(eta')| < 0.02 GeV/c^2 is applied in the ROOT analysis after the
    kinematic fit")
  .note(:efficiency_curve, "detection efficiency from the dedicated VMD-based eta' -> pi+pi-l+l-
    generator with infinite-width corrections and pseudoscalar mixing; form-factor
    (VMD factor) variation is a systematic uncertainty")
  .with_decay_card(decay_card_ee)
  .apply(sel_modeI)

# ==================== (II) J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu- ====================
alg_name_II = "EtapToPiPiMuMu"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
              cos_theta       0.93
              Vz              10.0
              Vr              1.0
              nChrp           "==2"
              nChrn           "==2"
              nNet            "==0"
          }
          .pid(method: :chi2_sum) {
              chi_min_cut 4
              identify :pion, :muon
          }
          .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              angle_to_track    15.0
              nGam              ">=1"
          }
          .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) {
              nominal
              constrain_four_momentum
              chi2_cut 75
          }
          .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) {
              use_track_index_from_nominal_kmfit
              constrain_four_momentum
          }

alg_modeII
  .note(:background_veto, "J/psi -> gamma 2(pi+ pi-) background rejected by requiring
    chi2_{gamma 2(pi+pi-)} > chi2_{gamma pi+ pi- mu+ mu-}; the stored competing chi^2
    (second fit block) is compared in the ROOT analysis")
  .note(:background_veto, "main backgrounds J/psi -> pi0 pi+ pi- pi+ pi- and
    J/psi -> gamma pi+ pi- pi+ pi-, estimated with inclusive J/psi MC")
  .note(:efficiency_curve, "upper limit on the number of signal events obtained from
    unbinned maximum-likelihood fits to M(pi+pi-mu+mu-) with an MC-determined eta'
    line shape and a second-order Chebychev background")
  .with_decay_card(decay_card_mumu)
  .apply(sel_modeII)

# ==================== (III) J/psi -> gamma eta', eta' -> gamma pi+ pi- ====================
alg_name_III = "EtapToGammaPiPi"
alg_modeIII = Algorithm.new(alg_name_III)
alg_modeIII.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeIII = Selection.new
sel_modeIII.select_track {
               cos_theta       0.93    # |cos(theta)| < 0.93
               Vz              10.0
               Vr              1.0
               nChrp           "==2"   # exactly four good charged tracks
               nChrn           "==2"
               nNet            "==0"
           }
           # No PID is applied in this mode.
           .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    15.0
               nGam              ">=2"   # radiative gamma and the gamma from eta'
           }
           # 4C kinematic fit under the J/psi -> pi+ pi- gamma gamma hypothesis;
           # for events with more than two photon candidates the combination with the
           # minimum chi^2 is retained. The pi0 veto M(gamma gamma) > 0.16 GeV/c^2
           # removes 94% of the background with an efficiency loss of 0.73%.
           .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
               nominal
               constrain_four_momentum
               invariant_mass_of(:gamma, :gamma).out_of(0.0, 0.16)
               chi2_cut 75
           }

alg_modeIII
  .note(:background_veto, "background events with a pi0 in the final state removed by
    the requirement M(gamma gamma) > 0.16 GeV/c^2 on the two selected photons")
  .note(:efficiency_curve, "the combination of gamma pi+ pi- invariant mass closest to
    the nominal eta' mass is chosen to reconstruct the eta'; the radiative photon from
    J/psi decays carries a unique energy of 1.4 GeV and is distinguished from the
    photon from eta' decays")
  .note(:background_veto, "eta' -> gamma pi+ pi- simulated with the rho0 -> pi+ pi-
    resonant contribution plus the non-resonant 'box anomaly' contribution, tuned to
    the measured di-pion mass spectrum")
  .with_decay_card(decay_card_gammapipi)
  .apply(sel_modeIII)

### Job submission ###
root_files_I   = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_ee])
root_files_II  = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])
root_files_III = alg_modeIII.execute_on([jpsi_data, jpsi_incMC, exMC_gammapipi])
