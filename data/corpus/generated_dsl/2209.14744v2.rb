# ============================================================
# Dataset preparation
# ============================================================
# Real-data points in the 4.23-4.70 GeV range (BOSS convention: <boss>_<Ecms in MeV>)
data_points = %w[
  703_4230 703_4237 703_4245 703_4246 703_4260 703_4270 703_4280
  703_4310 703_4360 703_4390 703_4420 703_4470 703_4530 703_4575
  703_4600 705_4290 705_4315 705_4340 705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
].map { |name| DatasetManager.real_data.find(name) }

# Corresponding inclusive MC samples
incmc_points = %w[
  703_4230 703_4237 703_4246 703_4260 703_4270 703_4280
  703_4360 703_4420 703_4600 706_4610 706_4620 706_4640
  706_4660 706_4680 706_4700
].map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay card: e+e- -> pi0 pi0 psi2(3823)
#                psi2(3823) -> gamma chi_c1
#                chi_c1     -> gamma J/psi
#                J/psi      -> l+ l-   (l = e or mu)
#                pi0        -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi0 pi0 psi(3823) PHSP;
  Enddecay

  Decay psi(3823)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  0.5000 e+   e-   PHOTOS VLL;
  0.5000 mu+  mu-  PHOTOS VLL;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC per energy point (same card / xs across the scan)
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "psipp_pi0pi0_psi3823"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ============================================================
# Event selection (BOSS)
# ============================================================
alg_name = "Psi2ToGammaChiC1"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})

event_selection = Selection.new
  # Exactly one positive and one negative charged track
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  # At least five good photons
  .select_photon {
    nGam              ">=5"
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  }
  # High-momentum lepton identification: p > 1.0 GeV -> lepton
  # (electron if EMC energy > 1.1 GeV, otherwise muon); require one l+ and one l-
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.1
    nlp "==1"
    nlm "==1"
  }
  # 4C kinematic fit on l+l- + five photons. One photon from psi2(3823) -> gamma chi_c1
  # is allowed to be missing (handled as a zero-mass missing particle in the fit,
  # i.e. the partial-reconstruction part of the analysis).
  .kinematic_fit([:lp, :lm, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    miss_track_of :gamma
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 15
  }

# Attach the signal decay card and render the BOSS selection
algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on data, inclusive MC and the per-point exclusive signal MC
algorithm.execute_on(data_points + incmc_points + exMC_signal)