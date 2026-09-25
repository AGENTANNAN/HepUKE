# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# The e+e- -> Sigma+ Sigma- Born-cross-section measurement uses the 41 BESIII XYZ
# energy points (3.510-4.951 GeV, 24.1 fb^-1) together with the matching inclusive MC.
scan_sample_names = %w[
  703_3510 703_3810 703_3900 703_4009 705_4130
  705_4160 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4237 703_4245 703_4246
  703_4260 703_4270 703_4280 705_4290 703_4310
  705_4315 705_4340 703_4360 705_4380 703_4390
  703_4420 705_4440 703_4470 703_4530 703_4575
  703_4600 706_4610 706_4620 706_4640 706_4660
  706_4680 706_4700 707_4740 707_4750 707_4780
  707_4946
]
data_points  = scan_sample_names.map { |name| DatasetManager.real_data.find(name) }
incMC_points = scan_sample_names.map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay card: e+e- -> Sigma+ Sigma- (KKMC + PHSP), with the full chain
# Sigma+ -> p+ pi0, anti-Sigma- -> anti-p- pi0, pi0 -> gamma gamma.
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

# 100k-event exclusive signal MC generated for every energy point
exMC_signals = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_sigmap_sigmam"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "SigmaPlusSigmaMinus"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.230]})   # central scan energy (3.510-4.951 GeV)
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     ">=1"                     # at least one positive track
    nChrn     ">=1"                     # at least one negative track
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start 0                     # EMC timing 0-700 ns (14 x 50 ns)
    tdc_emc_end 14
    angle_to_track 10.0                 # > 10 deg from any charged track
    energyThreshold_b 0.025             # E > 25 MeV (barrel)
    energyThreshold_e 0.050             # E > 50 MeV (endcap)
    nGam ">=4"                          # at least four photons
  }
  .pid(method: :probability) {          # PID: protons/anti-protons vs K and pi
    prob_cut 0.001                      # probability method, prob > 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"                          # at least one proton
    nprm ">=1"                          # at least one anti-proton
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: each gamma-gamma pair -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"                          # at least two pi0 candidates
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {  # 6C fit: 4C energy-momentum + two pi0 mass constraints
    nominal                             # nominal fit: corrected four-momenta are kept
    constrain_four_momentum
    chi2_cut 200
  }

# BOSS-side procedures that the DSL cannot express formally
algorithm
  .note(:pid_correction_method,
        "proton/anti-proton identification is realised in practice as a momentum " \
        "threshold p > 0.5 GeV/c combined with E/p < 0.8 to suppress the Bhabha " \
        "(e+e- -> e+e-) background; the momentum/E-p requirement cannot be " \
        "expressed by the DSL PID block directly")
  .note(:efficiency_curve,
        "the 4C constraint is evaluated with the CMS four-momentum of the individual " \
        "energy point; only a single ECMS constant is available, set to the central " \
        "scan energy (4.230 GeV), while the per-point beam energy is supplied by the framework")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)

root_files = algorithm.execute_on(data_points + incMC_points + exMC_signals)