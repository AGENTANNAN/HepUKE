### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at 3.097 GeV (BOSS 7.0.8)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive J/psi MC sample

# Decay card for the signal process: J/psi -> gamma eta eta, eta -> gamma gamma (5 photons, PHSP)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the background process: J/psi -> gamma eta pi0 pi0, eta -> gamma gamma, pi0 -> gamma gamma (5 photons, PHSP)
decay_card_background = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta pi0 pi0 PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_eta_eta"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

# Exclusive MC for the background process (500k events)
exMC_background = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_eta_pi0_pi0"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_background
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "GammaEtaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV

# Signal and background share the identical five-photon final state and selection,
# so a single Algorithm instance is used for both.
event_selection = Selection.new
event_selection.select_track {          # Charged-track selection: no charged tracks in the final state
                  cos_theta 0.92        # |cos(theta)| < 0.92
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm
                  nChrp     "==0"       # zero positive tracks
                  nChrn     "==0"       # zero negative tracks
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0       # EMC time window start (700 ns units)
                  tdc_emc_end       14      # EMC time window end
                  angle_to_track    10.0    # Min angle to any charged track (degrees)
                  energyThreshold_b 0.025   # > 25 MeV in the EMC barrel
                  energyThreshold_e 0.050   # > 50 MeV in the EMC endcap
                  nGam              ">=5"   # at least five photons
                }
               # No PID is applied: the final state is all neutral.
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {  # Nominal 4C fit to the five-photon hypothesis
                  nominal                 # mark as the nominal fit
                  constrain_four_momentum # constrain total four-momentum to the CMS energy
                  chi2_cut 200            # loose cut here; the paper's chi2 < 50 is applied at ROOT level
                }

# Generate the complete algorithm for the signal decay card (shared by both processes)
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute the algorithm on real data, inclusive MC, signal exclusive MC and background exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_background])