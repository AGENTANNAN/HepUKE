### Dataset preparation ###
# Seven real-data ψ(3686) scan points (BOSS 709)
data_points = [
  DatasetManager.real_data.find("709_3682"),
  DatasetManager.real_data.find("709_3683"),
  DatasetManager.real_data.find("709_3684"),
  DatasetManager.real_data.find("709_3685"),
  DatasetManager.real_data.find("709_3687"),
  DatasetManager.real_data.find("709_3691"),
  DatasetManager.real_data.find("709_3710")
]

# Corresponding inclusive MC samples
incMC_points = [
  DatasetManager.inclusive_mc.find("709_3682"),
  DatasetManager.inclusive_mc.find("709_3683"),
  DatasetManager.inclusive_mc.find("709_3684"),
  DatasetManager.inclusive_mc.find("709_3685"),
  DatasetManager.inclusive_mc.find("709_3687"),
  DatasetManager.inclusive_mc.find("709_3691"),
  DatasetManager.inclusive_mc.find("709_3710")
]

# Decay card for the signal process e+e- -> Sigma+ anti-Sigma-
# (psi(4260) as top mother, KKMC convention), Sigma+ -> p pi0,
# anti-Sigma- -> anti-p pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Sigma+ anti-Sigma- PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: same process generated at each scan energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_SigmaPairPol"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "SigmaPairPol"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93 (no Vz / Vr cut)
    nChrp     ">=1"     # at least one positive track
    nChrn     ">=1"     # at least one negative track
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start    0      # EMC timing lower edge
    tdc_emc_end      14     # EMC timing upper edge
    energyThreshold_b 0.025 # barrel energy threshold (25 MeV)
    energyThreshold_e 0.050 # endcap energy threshold (50 MeV)
    nGam  ">=4"             # at least four photons
  }
  .assign({:chrgp => :prp, :chrgn => :prm}) # no PID: all charged tracks are proton candidates
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.5" } # keep only p > 0.5 GeV/c
  .remove(:prm) { condition "three_momentum_of(:prm) < 0.5" } # keep only p > 0.5 GeV/c
  # Form two pi0 from photon pairs, each gamma-gamma mass constrained to m(pi0)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==2"    # exactly two pi0
  }
  # 6C kinematic fit to p pbar pi0 pi0 (4-momentum + two pi0 mass constraints)
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_signal)