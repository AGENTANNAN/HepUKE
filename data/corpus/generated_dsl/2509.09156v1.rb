# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")      # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card for the signal chain
# psi(3686) -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi-, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000  gamma eta(1405)          PHSP;
    Enddecay

    Decay eta(1405)
    1.000  f0(980) pi0              PHSP;
    Enddecay

    Decay f0(980)
    1.000  pi+ pi-                  PHSP;
    Enddecay

    Decay pi0
    1.000  gamma gamma              PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal chain: 5,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_gamma_eta1405_f0pi0"
    config.related_dataset = psip_data
    config.events          = 5_000_000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsiPToGammaEta1405"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # charged track selection
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        10.0                  # |Vz| < 10 cm along the beam direction
        Vr        1.0                   # Vxy < 1 cm in the transverse plane
        nChrp     ">=1"                 # at least two charged tracks overall ...
        nChrn     ">=1"                 # ... (one positive and one negative)
    }
    .select_photon {                    # photon selection
        tdc_emc_start     0             # TDC start time
        tdc_emc_end       14            # TDC end time
        energyThreshold_b 0.025         # E > 25 MeV in the EMC barrel
        energyThreshold_e 0.050         # E > 50 MeV in the EMC endcap
        angle_to_track    10.0          # >= 10 deg between photon and nearest charged track
        nGam              ">=3"         # at least three photons
    }
    .pid(method: :probability) {        # PID with the probability method
        prob_cut 0.001                  # PID probability > 0.001
        identify :pion, against: [:electron, :kaon, :proton]
        npip ">=1"                      # at least one pi+
        npim ">=1"                      # at least one pi-
    }
    # Nominal 4C kinematic fit to pi+ pi- gamma gamma gamma.
    # The gamma gamma pair closest to m_pi0 is the pi0; the remaining photon is the
    # radiative photon from psi(3686) -> gamma eta(1405).
    .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {
        nominal                         # nominal fit: its corrected four-momenta are used
        constrain_four_momentum         # 4C energy-momentum constraint to the CMS
        chi2_cut 20                     # chi^2 < 20
        invariant_mass_of(:gamma, :gamma).within(0.110, 0.160)     # |M(gg) - m_pi0| < 0.025 GeV
        invariant_mass_of(:gamma, :gamma).out_of(0.52, 0.57)       # veto M(gg) in [0.52, 0.57] GeV
        invariant_mass_of(:pip, :pim, :gamma).out_of(3.0469, 3.1469) # |M(pi+pi-gamma) - m_J/psi| > 0.05 GeV
    }

my_Algorithm
    # The omega veto |M(pi0 gamma) - m_omega| > 0.04 GeV/c^2 has no DSL form: the pi0 is not
    # a participant of the nominal 4C fit (whose final state is pi+ pi- gamma gamma gamma),
    # so invariant_mass_of(:pi0, :gamma) is not expressible inside the fit block.
    .note(:background_veto, "veto |M(pi0 gamma) - m_omega| > 0.04 GeV/c^2, where pi0 is the "
      "gamma-gamma pair at the pi0 mass; not expressible in the DSL because pi0 is not a "
      "participant of the nominal 4C kinematic fit (final state pi+ pi- gamma gamma gamma), "
      "so the veto is applied at ROOT level on the fitted four-momenta")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])