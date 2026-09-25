# =============================================================================
# e+e- -> gamma X(3872), X(3872) -> rho0 J/psi, rho0 -> pi+pi-,
# J/psi -> e+e- or mu+mu-
# Energy points: 4.009, 4.229 (4.230), 4.260, 4.360 GeV
# BOSS part: dataset preparation + event selection up to the final 4C fit
# =============================================================================

### Datasets ###
# Four real-data sets (BOSS 7.0.3) at the scan points, with matching inclusive MC.
data_4009 = DatasetManager.real_data.find("703_4009")   # ~4.008 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # ~4.226 GeV (4.229 GeV point)
data_4260 = DatasetManager.real_data.find("703_4260")   # ~4.258 GeV
data_4360 = DatasetManager.real_data.find("703_4360")   # ~4.358 GeV
energy_points = [data_4009, data_4230, data_4260, data_4360]

incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMCs = [incMC_4009, incMC_4230, incMC_4260, incMC_4360]

### Decay cards (EvtGen format) ###
# Signal: e+e- -> gamma X(3872), X(3872) -> rho0 J/psi, rho0 -> pi+pi-, J/psi -> e+e-
decay_card_signal_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 rho0 J/psi PHSP;
  Enddecay

  Decay rho0
  1.000 pi+ pi- VSS;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal: same chain but J/psi -> mu+mu-
decay_card_signal_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 rho0 J/psi PHSP;
  Enddecay

  Decay rho0
  1.000 pi+ pi- VSS;
  Enddecay

  Decay J/psi
  1.000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Dominant peaking background: e+e- -> pi+pi- J/psi (ISR photon from initial state),
# J/psi -> e+e-  (same visible final state as the e-channel signal)
decay_card_bkg_pipiJpsi = <<~DECAYCARD
  Decay psi(4260)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

### Exclusive MC: 100k events per energy point, per J/psi lepton mode, plus the ISR background ###
exMCs_signal_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gammaX3872_rhoJpsi_ee"    # one MC per energy point (auto-suffixed)
  config.events        = 100_000
  config.decay_card    = decay_card_signal_ee
  config.cross_section = :default
end

exMCs_signal_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gammaX3872_rhoJpsi_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_signal_mumu
  config.cross_section = :default
end

exMCs_bkg = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipiJpsi_isr_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_pipiJpsi
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Same final state (gamma pi+ pi- l+ l-) and identical selection for both lepton modes,
# so a single Algorithm covers the e-channel and the mu-channel (header built from the
# e-channel decay card).
alg_name = "GammaX3872RhoJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})   # beam energy default (per-point value set at job level)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                     # Charged track selection
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        10.0                    # |Vz| < 10 cm
    Vr        1.0                     # Vr < 1 cm (transverse plane)
    nChrp     "==2"                   # exactly 2 positive tracks
    nChrn     "==2"                   # exactly 2 negative tracks
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # Photon candidate selection
    tdc_emc_start     0               # EMC timing window 0-14
    tdc_emc_end       14
    angle_to_track    10.0            # > 10 deg from any charged track
    energyThreshold_b 0.025           # barrel energy > 25 MeV
    energyThreshold_e 0.050           # endcap energy > 50 MeV
    nGam              ">=1"           # at least one photon
  }
  .pid(method: :probability) {        # PID: probability method
    prob_cut 0.001                    # probability > 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6   # p>1.0 -> lepton; EMC E>0.6 -> e, else mu
    identify :pion, against: [:kaon, :proton]   # separate pions from kaons/protons
    npip "==1"                        # 1 pi+
    npim "==1"                        # 1 pi-
    nlp  "==1"                        # 1 l+
    nlm  "==1"                        # 1 l-
  }
  # 4C kinematic fit to gamma pi+ pi- l+ l- against the nominal beam four-momentum.
  # The fit iterates over photon/candidate combinations and keeps the smallest chi2.
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
    nominal                           # nominal fit; corrected four-momenta are saved
    constrain_four_momentum           # 4C energy-momentum constraint to CMS
    chi2_cut 60                       # chi2 < 60
    invariant_mass_of(:gamma, :pip, :pim).within(0.6, 20.0)   # loose gamma pi+ pi- (X(3872)) window
    invariant_mass_of(:lp, :lm).within(3.08, 3.12)            # J/psi mass window
  }

# BOSS-side procedures without a dedicated DSL construct
my_algorithm
  .note(:background_veto, "events with |cos(angle between the two pions)| >= 0.98 are rejected
    to suppress radiative Bhabha (e+e- -> gamma e+e-) and dimuon (e+e- -> gamma mu+mu-) backgrounds")
  .note(:radiative_photon_selection, "the highest-energy photon candidate is taken as the radiative
    photon; the 4C fit iterates over the remaining photon candidates and retains the combination
    with the smallest chi2 (default fit behaviour)")

# Render the selection into the BOSS algorithm package
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)

# Execute on the four energy points, their inclusive MC, and all exclusive MC samples
root_files = my_algorithm.execute_on(
  energy_points + incMCs +
  exMCs_signal_ee + exMCs_signal_mumu + exMCs_bkg
)