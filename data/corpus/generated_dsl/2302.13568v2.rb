# ============================================================
# J/psi (3.097 GeV) -> Sigma+ anti-Sigma-
#   Mode (a): Sigma+ -> p gamma , anti-Sigma- -> anti-p pi0 (pi0 -> gamma gamma)
#   Mode (b): Sigma+ -> p pi0   , anti-Sigma- -> anti-p gamma
# Both charge-conjugate modes share the same final state (p anti-p gamma gamma gamma)
# and the same selection chain -> one shared Algorithm.
# ============================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 3.097 GeV J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC sample

# Decay card for the signal mode Sigma+ -> p gamma, anti-Sigma- -> anti-p pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ gamma PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# One-million-event exclusive MC for the signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_sigmap_sigmam"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToSigmaSigmaBar"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

# Common selection chain for both charge-conjugate modes
event_selection = Selection.new
  .select_track {                 # Charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     ">=1"               # at least one positive track
    nChrn     ">=1"               # at least one negative track
  }
  .select_photon {                # Photon selection
    tdc_emc_start     0           # EMC TDC window start
    tdc_emc_end       14          # EMC TDC window end
    energyThreshold_b 0.025       # barrel energy threshold 25 MeV
    energyThreshold_e 0.050       # endcap energy threshold 50 MeV
    angle_to_track    10.0        # angle to nearest charged track > 10 degrees
    nGam              ">=3"       # at least three photons
  }
  .pid(method: :probability) {    # Proton identification (probability method)
    prob_cut 0.001                # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]  # p+ and anti-p- separated from K and pi
    nprp ">=1"                    # at least one proton candidate
    nprm ">=1"                    # at least one anti-proton candidate
  }
  # Nominal 5C fit: 4-momentum conservation + one gamma-gamma mass constraint to the pi0
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
  }
  # Competing 6C hypothesis (two pi0 mass constraints on the same photons); no chi2_cut,
  # no nominal -> the competing chi2 is stored for a later ROOT-level veto
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

# BOSS-side procedures that have no dedicated DSL construct
algorithm
  .note(:background_veto, "Delta(1232)+ -> p pi0 suppressed by requiring the proton decay-length significance L/sigma_L > 1.5")
  .note(:candidate_selection, "Single-tag anti-Sigma- candidate required to satisfy |M(anti-p pi0) - M(Sigma-)| < 13.5 MeV/c^2, best candidate chosen by minimum mass difference; if several photon combinations remain, the one with smallest 5C chi2 is kept")
  .note(:signal_region, "Double-tag signal region defined as 0.215 < p_p (Sigma+ rest frame) < 0.235 GeV/c; the competing Sigma+ -> p pi0 (two-pi0 6C) hypothesis is vetoed in ROOT when its chi2 is smaller than the nominal 5C chi2")

# Generate the algorithm for the signal process in the decay card
algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])