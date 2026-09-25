# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Real data at the three energy points: 3.773 GeV (psi(3770)), 3.650 GeV, and the 3.780 GeV scan point
data_3773 = DatasetManager.real_data.find("712_3773")
data_3650 = DatasetManager.real_data.find("709_3650")
data_3780 = DatasetManager.real_data.find("712_3780")

# Inclusive MC at 3.773 and 3.650 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# Signal decay card: e+e- -> psi(3770) -> p pbar pi0, pi0 -> gamma gamma (phase space)
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.000 p+ anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Background decay card: psi(3770) -> p pbar pi0 gamma, pi0 -> gamma gamma (phase space)
decay_card_bkg = <<~DECAYCARD
  Decay psi(3770)
  1.000 p+ anti-p- pi0 gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

energy_points = [data_3773, data_3650, data_3780]

# 200k-event exclusive signal MC generated on each energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name = "exmc_ppbarpi0_signal"
  config.events      = 200_000
  config.decay_card  = decay_card_signal
  config.cross_section = :default
end

# 20k-event exclusive background MC generated on each energy point
exMC_bkg = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name = "exmc_ppbarpi0gamma_bkg"
  config.events      = 20_000
  config.decay_card  = decay_card_bkg
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarPi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # nominal psi(3770) energy
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:multi_energy_ecms, "Identical selection is applied at 3.773, 3.650 and 3.780 GeV; ECMS is fixed to the nominal psi(3770) value 3.773 GeV here and must be re-configured per energy point when running the 3.650 and 3.780 GeV kinematic fits.")

event_selection = Selection.new
  .select_track {                # charged track selection
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     "==1"              # exactly one positive track
    nChrn     "==1"              # exactly one negative track
    nNet      "==0"              # net charge zero
  }
  .select_photon {               # photon selection
    tdc_emc_start     0          # EMC timing window start
    tdc_emc_end       14         # EMC timing window end
    energyThreshold_b 0.025      # barrel energy > 25 MeV (|cos(theta)| < 0.8)
    energyThreshold_e 0.050      # endcap energy > 50 MeV (0.86 < |cos(theta)| < 0.92)
    angle_to_track    10.0       # more than 10 degrees from the nearest charged track
    nGam              ">=2"      # at least two photon candidates
  }
  .pid(method: :probability) {   # PID via dE/dx and TOF (probability method)
    prob_cut 0.001               # confidence > 0.001
    identify :proton, against: [:pion, :kaon]   # identify proton and anti-proton against pi and K
    nprp ">=1"                   # at least one proton
    nprm ">=1"                   # at least one anti-proton
  }
  .remove(:prp) { condition "pt_of(:prp) < 0.3" }   # drop protons with pT < 0.3 GeV/c
  .remove(:prm) { condition "pt_of(:prm) < 0.3" }   # drop anti-protons with pT < 0.3 GeV/c
  .select_isolated_photon {      # isolated photon selection against (anti-)proton showers
    angle_to_prp_track 10.0      # angle to proton track > 10 degrees
    angle_to_prm_track 30.0      # angle to anti-proton track > 30 degrees
    nGam               ">=2"     # at least two isolated photons
  }
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {    # 5C kinematic fit to p pbar gamma gamma
    nominal                      # nominal fit (corrected four-momenta are used downstream)
    constrain_four_momentum      # 4C energy-momentum conservation
    invariant_mass_of(:gamma, :gamma).within(0.115, 0.155)              # +-3 sigma pi0 mass window on the photon pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # 1C pi0 mass constraint (all gamma gamma combinations iterated, smallest chi2 kept)
    chi2_cut 200                 # loose BOSS chi2 cut; the tight optimal cut is applied in the ROOT analysis
  }

# Generate the algorithm for the signal process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, signal exclusive MC and background exclusive MC
root_files = my_algorithm.execute_on([data_3773, data_3650, data_3780, incMC_3773, incMC_3650] + exMC_signal + exMC_bkg)