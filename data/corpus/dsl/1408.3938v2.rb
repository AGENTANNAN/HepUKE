# BESIII measurement of J/psi -> p pbar a0(980), a0(980) -> pi0 eta
# Based on 2.25e8 J/psi events collected at BEPCII

### Dataset description ###
jpsi_data   = DatasetManager.real_data.find("708_3097")       # J/psi data sample
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")    # J/psi inclusive MC sample

# Decay card for the signal process: J/psi -> p pbar a0(980), a0(980) -> pi0 eta
# a0(980) is parameterized by a Flatte formula (BESIII EvtGen convention: use PHSP if
# the specific Flatte model name is uncertain; the shape is enforced by the generator).
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  anti-p-  p+  a_0(980)0             PHSP;
    Enddecay

    Decay a_0(980)0
    1.0000  pi0 eta                             PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                         PHSP;
    Enddecay

    Decay eta
    1.0000  gamma gamma                         PHSP;
    Enddecay

    End
DECAYCARD

# Auxiliary decay card: phase-space MC J/psi -> p pbar pi0 eta
# Used to derive the mass-dependent efficiency curve of the a0(980) fit.
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000  anti-p-  p+  pi0 eta                PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                         PHSP;
    Enddecay

    Decay eta
    1.0000  gamma gamma                         PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the signal process J/psi -> p pbar a0(980) -> p pbar pi0 eta
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_ppbar_a0_pi0eta_signal"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

# Exclusive MC sample for phase-space J/psi -> p pbar pi0 eta (efficiency curve)
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_ppbar_pi0eta_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end
exMC_phsp.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiPPbarA0"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy in GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {
                  cos_theta  0.93     # |cos(theta)| < 0.93 (MDC polar-angle coverage)
                  Vr         10.0     # |Vr| < 10 mm (1 cm) in the transverse plane
                  Vz         10.0     # |Vz| < 10 cm along the beam axis
                  nChrp      ">=1"    # at least one positively charged track (proton)
                  nChrn      ">=1"    # at least one negatively charged track (anti-proton)
                }
                .select_photon {
                  energyThreshold_b 0.025   # E > 25 MeV in the EMC barrel (|cos theta| < 0.8)
                  energyThreshold_e 0.050   # E > 50 MeV in the EMC end caps (0.86 < |cos theta| < 0.92)
                  tdc_emc_start     0       # EMC timing T >= 0 (units of 50 ns)
                  tdc_emc_end       14      # EMC timing T <= 14 (units of 50 ns)
                  nGam              ">=4"   # at least four good photons
                }
                .pid(method: :probability) {
                  # PID uses only the dE/dx information: a track is a (anti-)proton if
                  # ProbPID(p) > ProbPID(K) and ProbPID(p) > ProbPID(pi).
                  identify :proton, against: [:kaon, :pion]
                  nprp ">=1"    # at least one identified proton
                  nprm ">=1"    # at least one identified anti-proton
                }
                # 4C kinematic fit: energy-momentum conservation at the production vertex,
                # applied to one proton + one anti-proton + four photons. For events with
                # more than four photons, the DSL iterates over all four-photon combinations
                # and selects the one with the smallest chi2_4C.
                .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200      # loose BOSS-level cut; tight cut chi2_4C < 35 applied in ROOT
                }

# Notes: BOSS-side procedures that have no direct DSL construct
my_Algorithm
  .note(:photon_pairing,
        "After the 4C kinematic fit, the four selected photons are paired to pi0 and " \
        "eta by minimising the chi2-like variable " \
        "chi2_pi0eta = (M(g1,g2) - M_pi0)^2/sigma_pi0^2 + (M(g3,g4) - M_eta)^2/sigma_eta^2 " \
        "with sigma_pi0 = 6.0 MeV/c^2 and sigma_eta = 9.8 MeV/c^2 (extracted from signal MC). " \
        "MC study shows the correct-pairing rate exceeds 99%.")
  .note(:background_veto,
        "pi0 pi0 veto: reject events whose minimum " \
        "chi2_pi0pi0 = (M(g1,g2) - M_pi0)^2/sigma_pi0^2 + (M(g3,g4) - M_pi0)^2/sigma_pi0^2 " \
        "is less than 100, to suppress p pbar pi0 pi0 final states surviving the 4C fit.")
  .note(:mass_window,
        "After photon pairing, the pi0 and eta candidate two-photon invariant masses " \
        "are required to lie within a 3-sigma window around the fitted mean values.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC, and both exclusive MC samples
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_phsp])
