# =============================================================================
# e+e- -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> l+ l-
#   Mode I  : J/psi -> e+ e-
#   Mode II : J/psi -> mu+ mu-
# Real data + inclusive MC : full scan 3.773 - 4.600 GeV
# =============================================================================

### ---------------------- Dataset description ---------------------- ###
# Full energy scan used in this analysis (BOSS_version_energy_tag)
scan_samples = %w[712_3773 703_4009 703_4180 703_4190 703_4200 703_4210
                  703_4220 703_4230 703_4246 703_4260 703_4270 703_4280
                  703_4360 703_4420 703_4600]

scan_data  = scan_samples.map { |tag| DatasetManager.real_data.find(tag) }   # real data of the scan
scan_incMC = scan_samples.map { |tag| DatasetManager.inclusive_mc.find(tag) } # matching inclusive MC

# Decay card - Mode I : J/psi -> e+ e-
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Decay card - Mode II : J/psi -> mu+ mu-
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC for each J/psi decay mode, at every scan point
exMCs_modeI = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_gammachic1_jpsi_ee"    # auto-suffixed per energy point
  config.events        = 100000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMCs_modeII = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_gammachic1_jpsi_mumu"  # auto-suffixed per energy point
  config.events        = 100000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

# Persist the MC configurations for later use
exMCs_modeI.each  { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }
exMCs_modeII.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### ---------------------- Event selection (BOSS) ---------------------- ###
# Selection common to both modes: charged tracks, soft-track removal, photons
event_selection_common = Selection.new
  .select_track {                      # charged-track selection
    cos_theta 0.93                     # |cos(theta)| < 0.93
    Vz        10.0                     # |Vz| < 10 cm
    Vr        1.0                      # Vr < 1 cm
    nChrp     ">=1"                    # at least one positively charged track
    nChrn     ">=1"                    # at least one negatively charged track
    nNet      "==0"                    # net charge zero
  }
  .remove(:chrgp) { condition "three_momentum_of(:chrgp) < 1.0" }  # drop tracks with p < 1.0 GeV
  .remove(:chrgn) { condition "three_momentum_of(:chrgn) < 1.0" }  # drop tracks with p < 1.0 GeV
  .select_photon {                     # photon selection
    tdc_emc_start     0                # TDC window 0-14 (700 ns units)
    tdc_emc_end       14
    angle_to_track    10.0             # >= 10 degrees from the nearest charged track
    energyThreshold_b 0.025            # 25 MeV in the barrel
    energyThreshold_e 0.050            # 50 MeV in the endcap
    nGam              ">=2"            # at least two photons
  }

# --------------------------------------------------------------------- #
# Mode I : J/psi -> e+ e-
# --------------------------------------------------------------------- #
alg_name_modeI = "GammaChiC1JpsiEE"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection_modeI = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.001                                     # PID probability cut
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0   # EMC energy > 1.0 GeV -> electron
    nlp ">=1"                                          # at least one e+
    nlm ">=1"                                          # at least one e-
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {         # 4C fit on gamma gamma l+ l-
    nominal
    invariant_mass_of(:lp, :lm).within(3.08, 3.12)     # J/psi mass window
    invariant_mass_of(:gamma, :gamma).out_of(0.5179, 0.5779)  # eta veto: |M(gg) - m_eta| > 0.03 GeV
    constrain_four_momentum
    chi2_cut 40                                        # chi2 < 40
  }

alg_modeI
  .note(:background_veto, "radiative Bhabha events vetoed: require cos(theta_{e gamma}) < 0.86 and |cos(theta_gamma)| < 0.8 for both photons; these angular criteria have no BOSS-DSL expression and are applied on the stored four-momenta")
  .note(:beam_energy_per_run, "the data set spans the 3.773-4.600 GeV scan; ECMS used by the 4C fit must be updated to the per-run beam energy of each energy point (only a nominal 4.260 GeV is declared here)")
  .with_decay_card(decay_card_modeI)
  .apply(event_selection_modeI)

# --------------------------------------------------------------------- #
# Mode II : J/psi -> mu+ mu-
# --------------------------------------------------------------------- #
alg_name_modeII = "GammaChiC1JpsiMuMu"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})

event_selection_modeII = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.001                                     # PID probability cut
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.4   # EMC energy > 0.4 GeV -> electron
    nlp ">=1"                                          # at least one mu+
    nlm ">=1"                                          # at least one mu-
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {         # 4C fit on gamma gamma l+ l-
    nominal
    invariant_mass_of(:lp, :lm).within(3.08, 3.12)     # J/psi mass window
    invariant_mass_of(:gamma, :gamma).out_of(0.1200, 0.1500)  # pi0 veto: |M(gg) - m_pi0| > 0.015 GeV
    constrain_four_momentum
    chi2_cut 40                                        # chi2 < 40
  }

alg_modeII
  .note(:background_veto, "pi+pi-pi0 background suppressed: require |cos(theta_gamma)| < 0.8 for both photons (the |M(gamma gamma) - m_pi0| > 0.015 GeV/c^2 part is applied as the invariant-mass veto in the kinematic-fit block); the angular criterion has no BOSS-DSL expression and is applied on the stored four-momenta")
  .note(:beam_energy_per_run, "the data set spans the 3.773-4.600 GeV scan; ECMS used by the 4C fit must be updated to the per-run beam energy of each energy point (only a nominal 4.260 GeV is declared here)")
  .with_decay_card(decay_card_modeII)
  .apply(event_selection_modeII)

### ---------------------- Execution ---------------------- ###
# The chi_c1,2 signal is finally extracted from the M(gamma J/psi) distribution (ROOT level)
root_files_modeI  = alg_modeI.execute_on(scan_data + scan_incMC + exMCs_modeI)
root_files_modeII = alg_modeII.execute_on(scan_data + scan_incMC + exMCs_modeII)