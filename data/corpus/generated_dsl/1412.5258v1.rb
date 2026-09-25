# ============================================================================
# J/psi -> eta Y(2175) (Y(2175) -> phi f0(980) -> phi pi+ pi-),
# J/psi -> phi f1(1285) (-> phi eta pi+ pi-),
# J/psi -> phi eta(1405) (-> phi eta pi+ pi-)
# Common final state: eta phi pi+ pi-  with eta -> gamma gamma, phi -> K+ K-
# i.e. K+ K- pi+ pi- gamma gamma
# ============================================================================

### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")       # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC (background)

# ---------------------------------------------------------------------------
# Decay cards (EvtGen format) for the three exclusive signal modes
# ---------------------------------------------------------------------------
# Mode 1 (nominal Y(2175)): J/psi -> eta Y(2175), Y(2175) -> phi f0(980), f0(980) -> pi+ pi-
decay_card_etaY2175 = <<~DECAYCARD
    Decay J/psi
    1.0000 eta Y(2175) PHSP;
    Enddecay

    Decay Y(2175)
    1.0000 phi f0(980) PHSP;
    Enddecay

    Decay f0(980)
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Mode 2: J/psi -> phi f1(1285), f1(1285) -> eta pi+ pi-
decay_card_phif1 = <<~DECAYCARD
    Decay J/psi
    1.0000 phi f1(1285) PHSP;
    Enddecay

    Decay f1(1285)
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Mode 3: J/psi -> phi eta(1405), eta(1405) -> eta pi+ pi-
decay_card_phieta1405 = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta(1405) PHSP;
    Enddecay

    Decay eta(1405)
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples (200k events each) for the three signal modes
# ---------------------------------------------------------------------------
exMC_etaY2175 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_etaY2175_phif0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_etaY2175
  config.cross_section   = :default
end

exMC_phif1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phif1_1285"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_phif1
  config.cross_section   = :default
end

exMC_phieta1405 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phieta1405"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_phieta1405
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All three modes share the same final state K+ K- pi+ pi- gamma gamma
# and hence the same selection; a single Algorithm is used.
alg_name  = "etaY2175phif0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = J/psi mass (GeV)

# Generator-level details that cannot be expressed in the HepScript decay-card
# syntax (P-wave eta-Y(2175), S-wave phi-f0(980) and pi+pi-, Flatte f0(980)
# lineshape with BESII parameters) are captured as a note for the MC production.
algorithm.note(:mc_generation_model,
  "nominal Y(2175) signal generated with a P-wave eta-Y(2175) decay, S-wave " \
  "phi-f0(980) and pi+pi- decays, and a Flatte f0(980) lineshape using BESII " \
  "parameters; PHSP is used as the decay-card placeholder since the dedicated " \
  "EvtGen model tokens are beyond the DSL surface")

event_selection = Selection.new
event_selection.select_track {
    cos_theta   0.93     # |cos(theta)| < 0.93 for charged tracks
    Vz          20.0     # |Vz| < 20 cm
    Vr          2.0      # Vr < 2 cm in the transverse plane
    nChrp       "==2"    # exactly two positively charged tracks
    nChrn       "==2"    # exactly two negatively charged tracks
    nNet        "==0"    # net charge zero
  }
  .select_photon {
    tdc_emc_start     0      # EMC timing window 0 ... 700 ns
    tdc_emc_end       14
    angle_to_track    10.0   # >10 degrees from the nearest charged track
    energyThreshold_b 0.025  # barrel (|cos(theta)| < 0.80) energy > 25 MeV
    energyThreshold_e 0.050  # endcap (0.86 < |cos(theta)| < 0.92) energy > 50 MeV
    nGam              ">=2"  # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001   # confidence (probability) cut
    # each track is assigned to the pi/K/p hypothesis with the highest confidence
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp   "==1"      # one K+
    nkm   "==1"      # one K-
    npip  "==1"      # one pi+
    npim  "==1"      # one pi-
  }
  # 4C kinematic fit to the K+ K- pi+ pi- gamma gamma hypothesis.
  # The fit iterates over all photon-pair combinations and keeps the smallest chi2.
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    # f0(980) mass region used in the fit: 0.90 < M(pi+ pi-) < 1.05 GeV/c^2
    invariant_mass_of(:pip, :pim).between(0.90, 1.05)
    chi2_cut 200
  }
# eta (0.528 < M(gamma gamma) < 0.566) and phi (1.006 < M(K+K-) < 1.032) windows
# are applied AFTER the kinematic fit -> handled in the ROOT analysis stage.

algorithm.with_decay_card(decay_card_etaY2175).apply(event_selection)

# Execute on real data, inclusive MC, and the three signal MC samples
root_files = algorithm.execute_on([jpsi_data, jpsi_incMC,
                                   exMC_etaY2175, exMC_phif1, exMC_phieta1405])