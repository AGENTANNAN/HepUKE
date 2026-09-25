### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC at 3.097 GeV

# Decay card (mode 1): J/psi -> omega K+ K- eta (phase space),
# omega -> pi+ pi- pi0 with Dalitz distribution, pi0 -> gamma gamma, eta -> gamma gamma
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000 omega K+ K- eta PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card (mode 2): J/psi -> omega X(1870) with X(1870) -> K+ K- eta,
# omega -> pi+ pi- pi0 with Dalitz distribution, pi0 -> gamma gamma, eta -> gamma gamma
decay_card_X1870 = <<~DECAYCARD
    Decay J/psi
    1.0000 omega X(1870) PHSP;
    Enddecay

    Decay X(1870)
    1.0000 K+ K- eta PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each of the two signal modes
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_omegaKKeta_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end

exMC_X1870 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_omegaX1870_KKeta"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_X1870
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Both signal modes share the identical final state K+ K- pi+ pi- gamma gamma gamma gamma,
# so a single algorithm instance serves both (shared-final-state case).
alg_name = "OmegaKKeta"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})  # J/psi center-of-mass energy (GeV)
            .note(:signal_lineshape,
                  "in the exclusive MC the omega is generated with the OMEGA_DALITZ model " \
                  "and the X(1870) is generated as a Breit-Wigner resonance")

event_selection = Selection.new
  .select_track {                    # Charged track selection
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        10.0                   # |Vz| < 10 cm
    Vr        1.0                    # Vr < 1 cm
    nChrp     "==2"                  # exactly two positive tracks
    nChrn     "==2"                  # exactly two negative tracks
    nNet      "==0"                  # net charge zero
  }
  .select_photon {                   # Photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0           # min angle to nearest charged track (deg)
    energyThreshold_b 0.025          # 25 MeV in the barrel
    energyThreshold_e 0.050          # 50 MeV in the endcap
    nGam              ">=4"          # at least four photons (2 from pi0, 2 from eta)
  }
  .pid(method: :probability) {       # PID by the probability method
    prob_cut 0.001                   # probability > 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K- vs pi and p
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
    nkp  "==1"                       # one K+
    nkm  "==1"                       # one K-
    npip "==1"                       # one pi+
    npim "==1"                       # one pi-
  }
  # 5C kinematic fit: K+ K- pi+ pi- + 4 gamma, 4-momentum conservation (4C)
  # plus one gamma-gamma pair constrained to the pi0 mass (1C)
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    nominal                          # nominal fit -> its corrected four-momenta are used
    constrain_four_momentum          # 4-momentum conservation
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # gamma-gamma -> pi0 mass
    chi2_cut 80                      # chi^2 < 80
  }
  # NOTE: the omega / eta mass windows and the eta' veto, and the search for X(1870)
  # in M(K+ K- eta), are applied on kinematic-fit-corrected variables and therefore
  # belong to the ROOT-level analysis, not to this BOSS selection.

my_Algorithm.with_decay_card(decay_card_phsp).apply(event_selection)

root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_phsp, exMC_X1870])