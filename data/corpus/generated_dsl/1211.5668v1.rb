### Dataset description ###
# J/psi (3.097 GeV) real data and inclusive MC, BESIII sample name 708_3097
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma X(1810), X(1810) -> omega phi,
# omega -> pi+ pi- pi0 (Dalitz), phi -> K+ K- (VSS), pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma X(1810) PHSP;
    Enddecay

    Decay X(1810)
    1.0000 omega phi PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_X1810_omega_phi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiGammaX1810OmegaPhi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        10.0                  # |Vz| < 10 cm
        Vr        1.0                   # Vr < 1 cm
        nChrp     "==2"                 # exactly two positive tracks
        nChrn     "==2"                 # exactly two negative tracks
        nNet      "==0"                 # net charge zero
    }
    .select_photon {                    # Photon selection
        tdc_emc_start     0             # EMC TDC window start
        tdc_emc_end       14            # EMC TDC window end
        angle_to_track    10.0          # > 10 deg from nearest charged track
        energyThreshold_b 0.025         # 25 MeV barrel threshold
        energyThreshold_e 0.050         # 50 MeV endcap threshold
        nGam              ">=3"         # at least three photons
    }
    .pid(method: :probability) {        # PID via probability method
        prob_cut 0.001                  # probability > 0.001
        identify :kaon, against: [:pion, :proton]   # K+/K- vs pions and protons
        identify :pion, against: [:kaon, :proton]   # pi+/pi- vs kaons and protons
        nkp  "==1"                      # exactly one K+
        nkm  "==1"                      # exactly one K-
        npip "==1"                      # exactly one pi+
        npim "==1"                      # exactly one pi-
    }
    # Nominal 4C kinematic fit to gamma gamma gamma K+ K- pi+ pi-
    .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
        nominal                          # nominal fit
        constrain_four_momentum          # 4C energy-momentum constraint
        chi2_cut 200                     # loose cut; tighter chi2 < 40 applied later in ROOT
    }
    # Competing hypothesis: 2-photon hypothesis (stores chi2 for ROOT-level veto)
    .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
        constrain_four_momentum
    }
    # Competing hypothesis: 4-photon hypothesis (stores chi2 for ROOT-level veto)
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
        constrain_four_momentum
    }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])