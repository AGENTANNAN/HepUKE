# Core DSL classes and dependencies are loaded automatically at execution

### Dataset preparation ###
# Six XYZ energy points (BOSS 7.0.3): 4.178, 4.189, 4.199, 4.209, 4.219 and 4.226 GeV
data_points = [
  DatasetManager.real_data.find("703_4180"),   # 4178 MeV
  DatasetManager.real_data.find("703_4190"),   # 4189 MeV
  DatasetManager.real_data.find("703_4200"),   # 4199 MeV
  DatasetManager.real_data.find("703_4210"),   # 4209 MeV
  DatasetManager.real_data.find("703_4220"),   # 4219 MeV
  DatasetManager.real_data.find("703_4230")    # 4226 MeV
]

# Corresponding inclusive MC sample at each energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the normalization mode e+e- -> Ds*+ Ds- (+ c.c.),
# Ds*+ -> gamma Ds+, Ds+ -> K+ K- pi+.
# Only the Ds+ side is reconstructed; the Ds- is treated inclusively
# (no decay block is given for it in the card).
decay_card_normalization = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds*+ Ds- PHSP;
    Enddecay

    Decay Ds*+
    1.0000 gamma Ds+ PHSP;
    Enddecay

    Decay Ds+
    1.0000 K+ K- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the normalization mode, one sample per energy point
exMC_normalization = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsstar_ds_norm"   # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_normalization
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "DsStarDsRecoil"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.226]})  # nominal point; scan covers 4.178 - 4.226 GeV

event_selection = Selection.new
  .select_track {                 # charged track selection
      cos_theta   0.93            # |cos(theta)| < 0.93
      Vz          10.0            # |Vz| < 10 cm
      Vr          1.0             # Vr < 1 cm in the transverse plane
      nChrp       ">=1"           # at least one positively charged track
  }
  .select_photon {                # photon selection (the soft gamma from Ds*+)
      tdc_emc_start     0         # EMC timing 0 - 14 (~0 - 700 ns)
      tdc_emc_end       14
      energyThreshold_b 0.025     # barrel energy > 25 MeV
      energyThreshold_e 0.050     # endcap energy > 50 MeV
      angle_to_track    10.0      # > 10 degrees from the nearest charged track
      nGam              ">=1"     # at least one photon
  }
  .pid(method: :probability) {    # probability PID with a 0.001 cut
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # separate pi from K and p
      identify :kaon, against: [:pion, :proton]   # separate K from pi and p
      nkp  ">=1"                  # at least one K+
      nkm  ">=1"                  # at least one K-
      npip ">=1"                  # at least one pi+
  }
  # Partial reconstruction: recID 2 = Ds- is left undetected, everything else
  # (Ds*+ 1, gamma 3, Ds+ 4, K+ 5, K- 6, pi+ 7) is tagged.  The recoil is taken
  # against the root-level tagged particle Ds*+ = Ds+ + gamma, so the window is
  # the required M_rec(Ds+ gamma) band.
  .partial_miss([2]) do
      require_recoil_mass 1.95, 1.99          # M_rec(Ds+ gamma) in (1.95, 1.99) GeV
  end

algorithm
  .note(:conexc_line_shape,
        "exclusive-MC line shape: the e+e- -> Ds*+ Ds- cross section over the six energy points " \
        "follows the ConExc parametrisation (ISR up to second order with vacuum-polarisation " \
        "correction); in the decay card the Ds- is treated inclusively, i.e. no decay block is " \
        "specified for it.")
  .note(:delta_e_selection,
        "of the surviving Ds*+ candidates the combination with the smallest |DeltaE| is kept, " \
        "DeltaE = (E_Ds+ + E_gamma + E_rec) - E0; the normalization-mode window " \
        "DeltaE in (-0.030, +0.020) GeV is applied to that best combination.  This " \
        "minimum-|DeltaE| combination choice has no dedicated DSL primitive.")
  .note(:kinematic_fit_cpp,
        "a 1C kinematic fit constraining M_rec(Ds+ gamma) to the nominal Ds- mass is implemented " \
        "directly in the C++ algorithm; partial reconstruction replaces the standard DSL " \
        "kinematic fit, so no kinematic_fit block is emitted here.")
  .note(:mass_bands,
        "in the same C++ step the mass bands M_rec(Ds+) in (2.100, 2.130) GeV and " \
        "M(Ds+ gamma) in (2.095, 2.130) GeV are required.")
  .note(:signal_mode_specifications,
        "the seven signal modes Ds+ -> K+ eta', eta' pi+, K+ eta, eta pi+, K+ K_S0, K_S0 pi+, " \
        "K+ pi0 are not part of this normalization-mode specification; each requires its own " \
        "Algorithm object and Selection chain.")
  .note(:other_channel_selections,
        "the per-mode DeltaE windows and the pi0 / eta / eta' / K_S0 selections used by the " \
        "other channels are listed in the analysis note; only the normalization-mode DeltaE " \
        "window is applied here.")
  .with_decay_card(decay_card_normalization)
  .apply(event_selection)

# Run the algorithm on the six real-data points, their inclusive MC and the
# per-point normalization-mode exclusive MC
root_files = algorithm.execute_on(data_points + incMC_points + exMC_normalization)