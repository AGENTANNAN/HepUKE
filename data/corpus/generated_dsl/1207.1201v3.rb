# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # 225.2 M J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# ---------------- Decay cards (EvtGen format) ----------------
# Mode 1 (signal): J/psi -> Lambda anti-Sigma0, anti-Sigma0 -> gamma anti-Lambda0
decay_card_Lambda_Sigma0 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Sigma0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 2 (signal, c.c.): J/psi -> anti-Lambda0 Sigma0, Sigma0 -> gamma Lambda0
decay_card_Lambdabar_Sigma0 = <<~DECAYCARD
    Decay J/psi
    1.0000 anti-Lambda0 Sigma0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 3 (signal): J/psi -> gamma eta_c, eta_c -> Lambda anti-Lambda0
decay_card_gamma_etac = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 4 (signal): J/psi -> Lambda0 anti-Lambda(1520)0, anti-Lambda(1520)0 -> gamma anti-Lambda0
decay_card_Lambda_Lambdabar1520 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda(1520)0 PHSP;
    Enddecay

    Decay anti-Lambda(1520)0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 5 (signal): J/psi -> Lambda(1520)0 anti-Lambda0, Lambda(1520)0 -> gamma Lambda0
decay_card_Lambda1520_Lambdabar = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda(1520)0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda(1520)0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 6 (background): J/psi -> Lambda0 anti-Lambda0
decay_card_Lambda_Lambdabar = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 7 (background): J/psi -> Sigma0 anti-Sigma0
decay_card_Sigma0_Sigmabar0 = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode 8 (background): J/psi -> Lambda0 anti-Lambda0 pi0
decay_card_Lambda_Lambdabar_pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# ---------------- Exclusive MC for each of the eight modes (100k events) ----------------
exMC_Lambda_Sigma0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambda_Sigmabar0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambda_Sigma0
  config.cross_section   = :default
end

exMC_Lambdabar_Sigma0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambdabar_Sigma0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambdabar_Sigma0
  config.cross_section   = :default
end

exMC_gamma_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etac"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_gamma_etac
  config.cross_section   = :default
end

exMC_Lambda_Lambdabar1520 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambda_Lambdabar1520"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambda_Lambdabar1520
  config.cross_section   = :default
end

exMC_Lambda1520_Lambdabar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambda1520_Lambdabar"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambda1520_Lambdabar
  config.cross_section   = :default
end

exMC_Lambda_Lambdabar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambda_Lambdabar"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambda_Lambdabar
  config.cross_section   = :default
end

exMC_Sigma0_Sigmabar0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Sigma0_Sigmabar0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Sigma0_Sigmabar0
  config.cross_section   = :default
end

exMC_Lambda_Lambdabar_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Lambda_Lambdabar_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_Lambda_Lambdabar_pi0
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All eight modes share the same final state gamma Lambda anti-Lambda -> gamma p pi- anti-p pi+,
# hence a single Algorithm instance with one common event selection is used.
alg_name = "gLLbar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {            # charged track selection
      cos_theta 0.93                      # |cos(theta)| < 0.93
      Vz        20.0                      # |Vz| < 20 cm
      Vr        10.0                      # Vr < 10 cm
      nChrp     "==2"                     # exactly two positive tracks
      nChrn     "==2"                     # exactly two negative tracks
      nNet      "==0"                     # net charge zero
  }
  .select_photon {                        # photon selection
      tdc_emc_start     0                 # EMC timing start
      tdc_emc_end       14                # EMC timing end
      angle_to_track    5.0               # at least 5 degrees from any charged track
      energyThreshold_b 0.025             # barrel energy > 25 MeV
      energyThreshold_e 0.050             # endcap energy > 50 MeV
      nGam              ">=1"             # require at least one photon
  }
  .pid(method: :probability) {            # particle identification for the Lambda daughters
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]   # identify p+ and p-bar (charge-conjugation shorthand)
      nprp ">=1"                          # at least one proton
      nprm ">=1"                          # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])       # remove identified (anti-)protons from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})       # treat all remaining tracks as pions
  .select_isolated_photon {               # isolated photon selection against (anti-)proton tracks
      angle_to_prp_track 5.0              # at least 5 degrees from the nearest proton
      angle_to_prm_track 30.0             # at least 30 degrees from the nearest anti-proton
      nGam ">=1"                          # require at least one isolated photon
  }
  .secondary_vertex_fit([:prp, :pim]) {   # reconstruct Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference   # pick the combination with mass closest to Lambda
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {   # reconstruct anti-Lambda -> anti-p pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {   # nominal 4C fit to the gamma Lambda anti-Lambda hypothesis
      nominal                       # this is the nominal fit; corrected four-momenta are stored
      constrain_four_momentum       # 4C energy-momentum conservation
      chi2_cut 200                  # loose BOSS cut; the paper's chi2 < 45 is applied later at ROOT level
      # The photon combination with the smallest chi2 is selected automatically.
  }

# One algorithm covers all eight modes (identical final state and selection);
# the decay card of the main signal channel defines the header kinematic variables.
my_Algorithm.with_decay_card(decay_card_Lambda_Sigma0).apply(event_selection)

# Execute on real data, inclusive MC, and all eight exclusive MC samples
root_files = my_Algorithm.execute_on([
  jpsi_data,
  jpsi_incMC,
  exMC_Lambda_Sigma0,
  exMC_Lambdabar_Sigma0,
  exMC_gamma_etac,
  exMC_Lambda_Lambdabar1520,
  exMC_Lambda1520_Lambdabar,
  exMC_Lambda_Lambdabar,
  exMC_Sigma0_Sigmabar0,
  exMC_Lambda_Lambdabar_pi0
])