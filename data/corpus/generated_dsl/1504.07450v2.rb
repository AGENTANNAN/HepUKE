# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3770) real data (round03 + round04, total ~2916.94 pb^-1) and inclusive MC
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay cards for the four signal modes: psi(3770) -> gamma chi_cJ, chi_cJ -> gamma J/psi, J/psi -> l+ l-
decay_card_chic1_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma chi_c1 PHSP;
    Enddecay
    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_chic2_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma chi_c2 PHSP;
    Enddecay
    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_chic1_mumu = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma chi_c1 PHSP;
    Enddecay
    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_chic2_mumu = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma chi_c2 PHSP;
    Enddecay
    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# Exclusive MC: 500k events for each of the four signal modes
exMC_chic1_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_chic1_ee"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_chic1_ee
  config.cross_section   = :default
end

exMC_chic2_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_chic2_ee"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_chic2_ee
  config.cross_section   = :default
end

exMC_chic1_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_chic1_mumu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_chic1_mumu
  config.cross_section   = :default
end

exMC_chic2_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_chic2_mumu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_chic2_mumu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# One independent Algorithm per signal mode (chi_c1/chi_c2 x e+e-/mu+mu-);
# all four share the identical selection chain (lepton flavour is handled
# automatically by identify_high_momentum_leptons).
alg_name_chic1_ee = "Chic1GammaEE"
alg_chic1_ee = Algorithm.new(alg_name_chic1_ee)
alg_chic1_ee.set_header(["#{alg_name_chic1_ee}Alg/#{alg_name_chic1_ee}.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .set_alias({"std::vector<double>" => "Vdouble"})

alg_name_chic2_ee = "Chic2GammaEE"
alg_chic2_ee = Algorithm.new(alg_name_chic2_ee)
alg_chic2_ee.set_header(["#{alg_name_chic2_ee}Alg/#{alg_name_chic2_ee}.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .set_alias({"std::vector<double>" => "Vdouble"})

alg_name_chic1_mumu = "Chic1GammaMuMu"
alg_chic1_mumu = Algorithm.new(alg_name_chic1_mumu)
alg_chic1_mumu.set_header(["#{alg_name_chic1_mumu}Alg/#{alg_name_chic1_mumu}.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .set_alias({"std::vector<double>" => "Vdouble"})

alg_name_chic2_mumu = "Chic2GammaMuMu"
alg_chic2_mumu = Algorithm.new(alg_name_chic2_mumu)
alg_chic2_mumu.set_header(["#{alg_name_chic2_mumu}Alg/#{alg_name_chic2_mumu}.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection chain shared by all four modes
event_selection_common = Selection.new
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        15.0                    # |Vz| < 15 cm
      Vr        1.0                     # Vr < 1 cm (transverse plane)
      nChrp     "==1"                   # exactly one positive track
      nChrn     "==1"                   # exactly one negative track
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start time
      tdc_emc_end       14              # EMC TDC end time
      angle_to_track    10.0            # > 10 deg from nearest charged track
      energyThreshold_b 0.050           # > 50 MeV in the barrel region
      energyThreshold_e 0.050           # > 50 MeV in the endcap region
      nGam              ">=2"           # at least two photons
  }
  .pid(method: :probability) {          # Particle identification
      prob_cut 0.001                    # PID probability > 0.001
      # high-momentum (p > 1.2 GeV/c) tracks treated as leptons; electron if EMC energy > 0.6 GeV, else muon
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,
                                     treat_as_electron_if_energy_above: 0.6
      nlp "==1"                         # exactly one l+
      nlm "==1"                         # exactly one l-
  }
  # 5C kinematic fit: gamma gamma l+ l-; 4-momentum conservation (4C) + M(l+l-) = M(J/psi) (1C)
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
      nominal                           # nominal fit — fitted four-momenta are stored
      constrain_four_momentum           # constrain total four-momentum to the c.m. values
      invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # J/psi mass constraint
      # veto gamma-gamma invariant masses in the pi0 and eta windows
      invariant_mass_of(:gamma, :gamma).out_of(0.124, 0.146)  # pi0 window
      invariant_mass_of(:gamma, :gamma).out_of(0.537, 0.558)  # eta window
      chi2_cut 200                      # loose chi^2 < 200 (tight cut applied in ROOT)
  }

# Each mode gets its own copy of the shared selection chain
alg_chic1_ee.with_decay_card(decay_card_chic1_ee).apply(event_selection_common.dup)
alg_chic2_ee.with_decay_card(decay_card_chic2_ee).apply(event_selection_common.dup)
alg_chic1_mumu.with_decay_card(decay_card_chic1_mumu).apply(event_selection_common.dup)
alg_chic2_mumu.with_decay_card(decay_card_chic2_mumu).apply(event_selection_common.dup)

# Execute the algorithms on real data, inclusive MC, and the matching signal MC
root_files_chic1_ee = alg_chic1_ee.execute_on([data_3773, incMC_3773, exMC_chic1_ee])
root_files_chic2_ee = alg_chic2_ee.execute_on([data_3773, incMC_3773, exMC_chic2_ee])
root_files_chic1_mumu = alg_chic1_mumu.execute_on([data_3773, incMC_3773, exMC_chic1_mumu])
root_files_chic2_mumu = alg_chic2_mumu.execute_on([data_3773, incMC_3773, exMC_chic2_mumu])