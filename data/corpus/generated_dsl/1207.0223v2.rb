# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data   = DatasetManager.real_data.find("709_3686")      # psi(3686) real data (~106M events, 160 pb^-1)
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")   # psi(3686) inclusive MC
cont_data   = DatasetManager.real_data.find("709_3650")      # 3.650 GeV continuum data (42 pb^-1)
cont_incMC  = DatasetManager.inclusive_mc.find("709_3650")   # 3.650 GeV continuum inclusive MC

# Decay card for the signal psi(3686) -> p pbar pi0, pi0 -> gamma gamma (phase-space)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the continuum e+e- -> gamma* -> p pbar pi0 (phase-space).
# No intermediate resonance is produced -> use psi(4260) as the KKMC top mother.
decay_card_continuum = <<~DECAYCARD
    Decay psi(4260)
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k phase-space signal events on the psi(3686)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ppbarpi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC: 100k phase-space continuum events at 3.650 GeV
exMC_continuum = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3650_ppbarpi0_continuum"
  config.related_dataset = cont_data
  config.events          = 100000
  config.decay_card      = decay_card_continuum
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "ppbarpi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})        # 3.686 GeV centre-of-mass energy
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
      cos_theta 0.8                            # |cos(theta)| < 0.8
      Vz        20.0                           # |Vz| < 20 cm
      Vr        2.0                            # Vr < 2 cm
      nChrp     "==1"                          # exactly one positive track (p)
      nChrn     "==1"                          # exactly one negative track (pbar)
      nNet      "==0"                          # net charge zero
    }
  .select_photon {                             # Photon selection
      tdc_emc_start     0                      # EMC timing window start
      tdc_emc_end       14                     # EMC timing window end
      angle_to_track    10.0                   # > 10 deg from nearest charged track (proton)
      energyThreshold_b 0.025                  # > 25 MeV in the barrel
      energyThreshold_e 0.050                  # > 50 MeV in the endcap
      nGam              ">=2"                  # at least two photons
    }
  .pid(method: :probability) {                 # Particle identification (probability method)
      prob_cut 0.001                           # PID probability > 0.001
      identify :proton, against: [:kaon, :pion] # separate protons from kaons and pions (both charges)
      nprp "==1"                               # exactly one proton
      nprm "==1"                               # exactly one antiproton
    }
  .select_isolated_photon {                    # Isolated photon selection (reject p/pbar shower remnants)
      angle_to_prp_track 10.0                  # >= 10 deg from the proton track
      angle_to_prm_track 30.0                  # >= 30 deg from the antiproton track
      nGam ">=2"                               # at least two photons remain
    }
  # 4C kinematic fit to the p pbar gamma gamma system (stores chi2_4C for later comparison)
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
      constrain_four_momentum                  # 4-momentum conservation
    }
  # Nominal 5C fit: 4C + pi0 mass constraint on the two photons;
  # the fit automatically keeps the smallest-chi2 combination among photon pairs.
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 5C (pi0 mass constraint)
      chi2_cut 200                             # loose chi2 cut (tight cut applied in ROOT)
    }

# Generate the algorithm for the signal decay card and render the selection
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on all datasets (signal & continuum, data & MC, plus the two exclusive MC samples)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_signal, exMC_continuum])