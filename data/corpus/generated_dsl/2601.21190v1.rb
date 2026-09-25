# =============================================================================
# BOSS event-selection DSL
# Search for e+e- -> eta eta psi(2S)  and  e+e- -> eta psi0(4360) -> eta eta psi(2S)
# at sqrt(s) = 4.840, 4.918 and 4.951 GeV
#
# Partial reconstruction of the final state: only ONE eta -> gamma gamma is
# reconstructed (via a 1C Kalman fit); the second eta is left missing and is
# handled by the 4C energy-momentum constraint of the nominal kinematic fit.
# psi(2S) -> pi+ pi- J/psi,  J/psi -> e+ e-  or  mu+ mu-.
# =============================================================================

### ---------------------- Dataset preparation ---------------------- ###
# Real data at the three energy points (0.9 fb^-1 total, BOSS 7.0.7)
data_4840 = DatasetManager.real_data.find("707_4840")   # 4.840 GeV
data_4914 = DatasetManager.real_data.find("707_4914")   # 4.918 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # 4.951 GeV
# Matching inclusive MC samples at the same energies
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

### Decay cards (EvtGen) ###
# Mode 1 : e+e- -> eta eta psi(2S), J/psi -> e+ e-
decay_card_direct_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0  eta eta psi(2S)  PHSP;
  Enddecay

  Decay psi(2S)
  1.0  pi+ pi- J/psi  PHSP;
  Enddecay

  Decay J/psi
  1.0  e+ e-  PHOTOS VLL;
  Enddecay

  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2 : e+e- -> eta eta psi(2S), J/psi -> mu+ mu-
decay_card_direct_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0  eta eta psi(2S)  PHSP;
  Enddecay

  Decay psi(2S)
  1.0  pi+ pi- J/psi  PHSP;
  Enddecay

  Decay J/psi
  1.0  mu+ mu-  PHOTOS VLL;
  Enddecay

  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode 3 : e+e- -> eta psi0(4360) -> eta eta psi(2S), J/psi -> e+ e-
decay_card_res_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0  eta psi(4360)  PHSP;
  Enddecay

  Decay psi(4360)
  1.0  eta psi(2S)  PHSP;
  Enddecay

  Decay psi(2S)
  1.0  pi+ pi- J/psi  PHSP;
  Enddecay

  Decay J/psi
  1.0  e+ e-  PHOTOS VLL;
  Enddecay

  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode 4 : e+e- -> eta psi0(4360) -> eta eta psi(2S), J/psi -> mu+ mu-
decay_card_res_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0  eta psi(4360)  PHSP;
  Enddecay

  Decay psi(4360)
  1.0  eta psi(2S)  PHSP;
  Enddecay

  Decay psi(2S)
  1.0  pi+ pi- J/psi  PHSP;
  Enddecay

  Decay J/psi
  1.0  mu+ mu-  PHOTOS VLL;
  Enddecay

  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC (100k events each) ###
# Modes 1 & 2 : the direct eta eta psi(2S) process, generated at ALL THREE points
exMC_direct_ee = DatasetManager.create_exclusive_mc_for([data_4840, data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_etaetapsi2s_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_direct_ee
  config.cross_section = :default
end

exMC_direct_mumu = DatasetManager.create_exclusive_mc_for([data_4840, data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_etaetapsi2s_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_direct_mumu
  config.cross_section = :default
end

# Modes 3 & 4 : the resonant eta psi0(4360) process, generated at the TWO higher
# energy points only
exMC_res_ee = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_etapsi4360_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_res_ee
  config.cross_section = :default
end

exMC_res_mumu = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_etapsi4360_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_res_mumu
  config.cross_section = :default
end

### ---------------------- Event selection (BOSS) ---------------------- ###
# All four signal modes share the same final state (eta + pi+ pi- + l+ l- with a
# missing second eta) and the same selection, so a single Algorithm is used, as in
# the shared-final-state case.
alg_name = "EtaEtaPsi2S"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.918]})   # nominal CMS energy (per-run beam energy used by the fit)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  # Charged tracks : four tracks, two positive, two negative, net charge zero
  .select_track {
      cos_theta   0.93    # |cos(theta)| < 0.93
      Vz          10.0    # |Vz| < 10 cm
      Vr          1.0     # Vr < 1 cm
      nChrp       "==2"   # exactly two positive tracks
      nChrn       "==2"   # exactly two negative tracks
      nNet        "==0"   # net charge zero
  }
  # Photons : at least two, 25 MeV (barrel) / 50 MeV (endcap) energy threshold
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"
  }
  # PID (probability method) : pions separated from kaons; high-momentum tracks
  # treated as leptons (electron / muon)
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon]
      npip ">=1"
      npim ">=1"
  }
  # 1C Kalman fit : constrain the gamma gamma pair to the eta mass and require at
  # least one eta with M(gamma gamma) in [499.5, 576.9] MeV
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      invariant_mass_of(:gamma, :gamma).within(0.4995, 0.5769)   # M(gamma gamma) window
      chi2_cut 200
      neta ">=1"
  }
  # Nominal 4C kinematic fit of the reconstructed eta pi+ pi- l+ l- system, with the
  # second eta missing (energy-momentum conservation to the CMS)
  .kinematic_fit([:eta, :pip, :pim, :lp, :lm]) {
      nominal
      constrain_four_momentum
      miss_track_of(:eta)   # the second eta is not reconstructed (partial reconstruction)
      chi2_cut 200
  }

# BOSS-side procedures that cannot be expressed with the current DSL
my_algorithm
  .note(:muon_muc_cut, "at least one muon is required to have MUC hit depth > 30 cm
    to suppress pi -> mu misidentification; the MUC hit-depth variable is not
    available in the DSL track/PID blocks")
  .note(:lepton_pid_threshold, "high-momentum leptons follow electron if EMC energy
    > 0.6 GeV and muon if EMC energy < 0.4 GeV; identify_high_momentum_leptons splits
    electron/muon at a single energy value (0.6 GeV), so the muon-side threshold is
    applied in the ROOT analysis")

my_algorithm.with_decay_card(decay_card_direct_ee).apply(event_selection)

# Execute on the three real-data sets, their matching inclusive MC, and all four
# exclusive-MC groups
all_datasets = [data_4840, data_4914, data_4946,
                incMC_4840, incMC_4914, incMC_4946] +
               exMC_direct_ee + exMC_direct_mumu + exMC_res_ee + exMC_res_mumu

root_files = my_algorithm.execute_on(all_datasets)