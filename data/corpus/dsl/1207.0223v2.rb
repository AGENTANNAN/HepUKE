# BESIII paper arXiv:1207.0223v2
# Partial wave analysis of psi(3686) -> p pbar pi0, based on 106 x 10^6 psi(3686)
# events (160 pb^-1) collected with BESIII.
# Final state: p pbar pi0 (pi0 -> gamma gamma) -> two charged tracks + two photons.
# Nominal fit: 5C kinematic fit (4C plus the pi0 mass constraint on the two photons).

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data, 106 M events (160 pb^-1)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC (10^8 events) for background studies
off_data   = DatasetManager.real_data.find("709_3650")     # continuum data at sqrt(s) = 3.650 GeV (42 pb^-1)
off_incMC  = DatasetManager.inclusive_mc.find("709_3650")

### Decay cards (EvtGen format) ###
# Signal: psi(3686) -> p pbar pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000  p+  anti-p-  pi0         PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Continuum background: e+e- -> gamma* -> p pbar pi0 (studied with the 3.650 GeV data)
decay_card_continuum = <<~DECAYCARD
    Decay psi(4260)
    1.0000  p+  anti-p-  pi0         PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_ppbar_pi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'exMC_psip_ppbar_pi0')

exMC_continuum = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "continuum_gammastar_to_ppbar_pi0"
  config.related_dataset = off_data
  config.events          = 100_000
  config.decay_card      = decay_card_continuum
  config.cross_section   = :default
end
exMC_continuum.save_to_config(format: :yaml, file_path: 'exMC_continuum_ppbar_pi0')

### Event selection (BOSS) ###
alg_name = "PsipToPpbarPi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.8      # |cos(theta)| < 0.8 for each charged track
                  Vz        20.0     # closest approach within +/- 20 cm of the IP along the beam direction
                  Vr        2.0      # closest approach within 2 cm of the beam axis
                  nChrp     "==1"    # exactly two charged tracks of opposite charge
                  nChrn     "==1"
                  nNet      "==0"
                }
               .select_photon {
                  tdc_emc_start     0     # EMC timing requirement
                  tdc_emc_end      14
                  angle_to_track   10.0   # angle between photon candidate and the proton > 10 degrees
                  energyThreshold_b 0.025 # E_min = 25 MeV in the barrel EMC
                  energyThreshold_e 0.050 # E_min = 50 MeV in the endcap EMC
                  nGam             ">=2"  # at least two photons from the pi0 -> gamma gamma decay
                }
               .select_isolated_photon {
                  angle_to_prp_track 10.0  # angle between photon candidate and proton > 10 degrees
                  angle_to_prm_track 30.0  # stricter 30 degree cut against the anti-proton,
                                           # excluding photons from anti-proton annihilation
                  nGam               ">=2"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # largest PID probability assigned per track:
                                                             # one track proton, the other anti-proton
                  nprp "==1"
                  nprm "==1"
                }
               # 4C kinematic fit: sum of four-momenta of all particles constrained to the
               # energy and three momentum components of the initial e+e- system.
               .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
                  constrain_four_momentum
               }
               # 5C kinematic fit: one additional constraint of the pi0 mass on the two photons,
               # providing more accurate momentum information for the final states. When more than
               # two photons are found, all possible p pbar gamma gamma combinations are considered
               # and the one yielding the smallest chi2_5C is retained (automatic in the DSL).
               .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200   # loose BOSS cut; the optimal chi2_5C requirement is applied in ROOT
               }

alg.note(:jpsi_veto,
         "A cut of |M(p pbar) - M(J/psi)| > 40 MeV/c^2 is applied to exclude events with a p pbar pair " \
         "arising from J/psi decay: the observed J/psi width is far larger than its natural width due " \
         "to detector resolution, which would cause problems in the partial wave analysis. This cut " \
         "uses the kinematic-fit-corrected four-momenta and is therefore applied at the ROOT level.")
   .note(:pi0_sideband,
         "Background from psi(3686) decays is estimated with pi0 sideband events defined by " \
         "30 MeV/c^2 < |M(gamma gamma) - 135 MeV/c^2| < 45 MeV/c^2. Only 26 events are found in the " \
         "sideband region; the mass window is evaluated with fit-corrected four-momenta and belongs " \
         "to the ROOT-level analysis.")
   .note(:background_estimation,
         "Two background sources are studied: (1) psi(3686) decays, estimated both from 10^8 " \
         "MC-simulated psi(3686) events (40 events survive the selection, mainly due to " \
         "misidentified or lost photons) and from the pi0 sideband (26 events); (2) the continuum " \
         "process e+e- -> gamma* -> p pbar pi0, studied with the 42 pb^-1 of data at sqrt(s) = 3650 MeV, " \
         "which normalizes to 447 background events. The continuum contribution accounts for about " \
         "95% of the total background.")
   .note(:event_yield,
         "4988 events survive the event selection criteria; the signal efficiency from MC is 25.8% " \
         "and the branching fraction obtained is B(psi(3686) -> p pbar pi0) = (1.65 +/- 0.03 +/- 0.15) x 10^-4.")
   .note(:partial_wave_analysis,
         "A partial wave analysis is performed with an unbinned maximum likelihood fit. The " \
         "amplitudes A_i for all possible partial waves are constructed with the relativistic " \
         "covariant tensor amplitude formalism, and the total transition probability per event is " \
         "omega = |sum_i c_i A_i|^2. The free parameters c_i are determined by maximizing " \
         "ln(L) = sum_i ln(omega(xi_i) epsilon(xi_i) / integral dxi omega(xi) epsilon(xi)), where " \
         "epsilon(xi) is the detection efficiency. Background from the pi0 sideband and continuum " \
         "processes is removed by subtracting the corresponding log-likelihood values; no " \
         "interference between continuum processes and psi(3686) decays is considered. Each N* is " \
         "parameterized with a Breit-Wigner function with free mass and width. Nineteen intermediate " \
         "resonances (including N(940), N(1440), N(1520), N(1535), N(1650), N(1720), N(1885), " \
         "N(2065), N(2090), N(2100) and phase space) are considered; eight states are found to be " \
         "significant, among them the two new resonances N(2300) (1/2+) and N(2570) (5/2-). This fit " \
         "is not expressible in the DSL.")
   .note(:systematic_uncertainties,
         "Common systematic uncertainties applicable to all branching fraction measurements: number " \
         "of psi(3686) events 4%, MDC tracking 4% (two charged tracks), particle identification 2% " \
         "(both proton and anti-proton), photon detection efficiency 2%, kinematic fit 7%; total from " \
         "these common sources 9.4%. A second category concerns the fitting procedure: additional " \
         "possible resonances, different Breit-Wigner parameterizations for the partial wave " \
         "amplitudes, background estimation, the J/psi exclusion cut, and the differences in the " \
         "Input-Output check; these apply to the mass, width and branching fraction measurements of " \
         "the intermediate states.")

alg.with_decay_card(decay_card_signal).apply(event_selection)

root_files = alg.execute_on([psip_data, psip_incMC, off_data, off_incMC, exMC_signal, exMC_continuum])
