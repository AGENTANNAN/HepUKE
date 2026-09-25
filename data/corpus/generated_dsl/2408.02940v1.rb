# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # 2.712-billion-event psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")      # 3.650 GeV continuum real data
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")   # 3.650 GeV continuum inclusive MC

# Decay card for the full signal chain (EvtGen format):
#   psi(2S) -> gamma eta_c(2S), eta_c(2S) -> K+ K- eta, eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 K+ K- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive phase-space MC sample for the full signal chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi2S_gamma_etac2S_KKeta_exclusive_mc"
  config.related_dataset = psip_data     # associated real dataset for matching conditions
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaEtac2SKKeta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})         # sqrt(s) = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                                             # Charged track selection
    cos_theta 0.93                                            # |cos(theta)| < 0.93
    Vz        10.0                                            # |Vz| < 10 cm
    Vr        1.0                                             # Vr < 1 cm
    nTot      "==2"                                           # exactly two charged tracks
    nChrp     ">=1"                                           # at least one positive track
    nChrn     ">=1"                                           # at least one negative track
    nNet      "==0"                                           # net charge zero
  }
  .select_photon {                                            # Photon selection
    tdc_emc_start     0                                       # EMC TDC window start
    tdc_emc_end       14                                      # EMC TDC window end
    angle_to_track    10.0                                    # > 10 deg from any charged track
    energyThreshold_b 0.025                                   # > 25 MeV (barrel region)
    energyThreshold_e 0.040                                   # > 40 MeV (endcap region; low threshold keeps the M1 radiative photon)
    nGam ">=3 && nGam<=6"                                     # between 3 and 6 photons
  }
  .pid(method: :probability) {                                # Particle identification
    prob_cut 0.001                                            # pid probability > 0.001
    identify :kaon, against: [:pion, :proton]                 # K+ and K- against pion / proton
    nkp "==1"
    nkm "==1"
  }
  # 4C kinematic fit of the K+ K- gamma gamma gamma (3-photon) hypothesis -> chi2(3gamma)
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    chi2_cut 20
  }
  # Competing 4C fit of the K+ K- gamma gamma gamma gamma (4-photon) hypothesis.
  # No chi2_cut / no nominal -> stores chi2(4gamma) for the ROOT-level veto chi2(3gamma) < chi2(4gamma)
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Nominal 5C fit: 4C + eta mass constraint on the two eta-daughter photons
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).within(0.51, 0.57)                  # eta signal window 0.51 < M(gamma gamma) < 0.57 GeV
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # 5C eta mass constraint
    chi2_cut 200
  }

# BOSS-side procedures that cannot be expressed in the selection DSL
my_algorithm
  .note(:background_veto, "competing 4-photon 4C hypothesis chi2 is stored; ROOT-level veto chi2(3gamma) < chi2(4gamma) suppresses psi(3686) -> gamma gamma K+ K- eta")
  .note(:photon_assignment, "the lowest-energy photon among the three is taken as the radiative photon from psi(2S) -> gamma eta_c(2S)")
  .note(:pi0_veto, "2D pi0 veto applied; eta' photon-energy veto rejecting photons with energy outside [0.156, 0.196] GeV; M(3gamma) > 0.6 GeV; M(K+K-) < 3.0 GeV and M(K+K-) outside [1.007, 1.033] GeV")
  .note(:track_momentum, "charged track momentum <= 2 GeV/c")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC (both energies) and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_signal])