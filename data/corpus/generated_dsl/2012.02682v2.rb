# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# 14 BOSS 703 data sets collected between 4.178 and 4.600 GeV (energy scan)
data_4180 = DatasetManager.real_data.find("703_4180")   # 4178 MeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4188.8 MeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4198.9 MeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4209.2 MeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4218.7 MeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4226.3 MeV
data_4237 = DatasetManager.real_data.find("703_4237")   # 4235.7 MeV
data_4246 = DatasetManager.real_data.find("703_4246")   # 4243.8 MeV
data_4260 = DatasetManager.real_data.find("703_4260")   # 4258.0 MeV
data_4270 = DatasetManager.real_data.find("703_4270")   # 4266.8 MeV
data_4280 = DatasetManager.real_data.find("703_4280")   # 4277.7 MeV
data_4360 = DatasetManager.real_data.find("703_4360")   # 4358.3 MeV
data_4420 = DatasetManager.real_data.find("703_4420")   # 4415.6 MeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.5 MeV

# Corresponding inclusive MC samples at each energy point
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4237 = DatasetManager.inclusive_mc.find("703_4237")
incMC_4246 = DatasetManager.inclusive_mc.find("703_4246")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4270 = DatasetManager.inclusive_mc.find("703_4270")
incMC_4280 = DatasetManager.inclusive_mc.find("703_4280")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# All real data points (passed together to the energy-scan MC builder)
data_points = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230,
               data_4237, data_4246, data_4260, data_4270, data_4280, data_4360,
               data_4420, data_4600]
incMC_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230,
                incMC_4237, incMC_4246, incMC_4260, incMC_4270, incMC_4280, incMC_4360,
                incMC_4420, incMC_4600]

# Decay cards for the three signal modes: e+e- -> chi_cJ pi+ pi-, chi_cJ -> gamma J/psi, J/psi -> e+ e-
decay_card_chic0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c0 pi+ pi- PHSP;
    Enddecay

    Decay chi_c0
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c1 pi+ pi- PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c2 pi+ pi- PHSP;
    Enddecay

    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each chi_cJ mode at every energy point (one MC per data point).
# The mu+mu- final state is NOT generated separately - it is treated as a competing fit hypothesis.
exMC_chic0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_chic0_pipijpsi_ee"
  config.events        = 500000
  config.decay_card    = decay_card_chic0
  config.cross_section = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_chic1_pipijpsi_ee"
  config.events        = 500000
  config.decay_card    = decay_card_chic1
  config.cross_section = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_chic2_pipijpsi_ee"
  config.events        = 500000
  config.decay_card    = decay_card_chic2
  config.cross_section = :default
end

### Event selection (BOSS) ###
# The three chi_cJ share the identical final state (pi+ pi- e+ e- gamma) and selection,
# so a single Algorithm instance covers all three modes (separated only at the ROOT level).
alg_name = "ChiCJPiPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})   # representative CMS energy of the scan
            .set_alias({"std::vector<double>" => "Vdouble"})
            # The E/p-based e/mu separation and the additional background vetoes have no
            # dedicated DSL primitive and are recorded for the systematic-uncertainty step.
            .note(:pid_e_over_p, "e/mu separation in PID uses E_EMC/P_MDC > 0.7 for electrons and
              < 0.3 for muons, together with the momentum separation (pion p < 1.0 GeV/c,
              lepton p > 1.0 GeV/c); the DSL identify_high_momentum_leptons primitive applies an
              EMC-energy based lepton/e discrimination instead, so the analytic E/p criteria are
              applied in the generated BOSS code")
            .note(:background_veto, "after the 5C kinematic fits the following background vetoes are
              applied: cos(alpha_pip_pim) < 0.98; vetoes against eta/eta' J/psi and omega
              intermediate states; and recoil-mass windows against psi(2S) / X(3872)")

# Full event selection chain
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
      cos_theta  0.93                   # |cos(theta)| < 0.93
      Vz         10.0                   # |Vz| < 10 cm
      Vr         1.0                    # Vr < 1 cm
      nChrp      "==2"                  # exactly two positively charged tracks
      nChrn      "==2"                  # exactly two negatively charged tracks
      nNet       "==0"                  # net charge zero
    }
  .select_photon {                      # Photon selection (the gamma from chi_cJ -> gamma J/psi)
      tdc_emc_start      0              # EMC timing start
      tdc_emc_end        14             # EMC timing end (14 x 50 ns = 700 ns)
      angle_to_track     20.0           # > 20 degrees from the nearest charged track
      energyThreshold_b  0.025          # barrel E > 25 MeV (|cos(theta)| < 0.80)
      energyThreshold_e  0.050          # endcap E > 50 MeV (0.86 < |cos(theta)| < 0.92)
      nGam               ">=1"          # at least one good photon
    }
  .pid(method: :probability) {          # PID, probability method
      prob_cut 0.001                    # probability > 0.001
      # momentum-based separation: pion p < 1.0 GeV/c, lepton p > 1.0 GeV/c
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]   # pi+ and pi- against K / p
      npip "==1"                        # exactly one pi+
      npim "==1"                        # exactly one pi-
      nlp  "==1"                        # exactly one lepton (l+)
      nlm  "==1"                        # exactly one lepton (l-)
    }
  # Nominal 5C fit: 4-momentum conservation + J/psi mass constraint on pi+ pi- l+ l- gamma
  .kinematic_fit([:pip, :pim, :lp, :lm, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
      chi2_cut 50                       # chi2_5C < 50
    }
  # Competing hypothesis: refit the same event as pi+ pi- mu+ mu- gamma (muon mass hypothesis),
  # store its chi2_5C for the ROOT-level best-candidate selection (no chi2_cut, not nominal)
  .assign({:lp => :mup, :lm => :mum})   # treat the identified leptons as muons for the second fit
  .kinematic_fit([:pip, :pim, :mup, :mum, :gamma]) {
      constrain_four_momentum
      invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
    }

# Generate the algorithm for the shared final state, then execute on all data / MC samples
my_algorithm.with_decay_card(decay_card_chic0).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points +
                                     exMC_chic0 + exMC_chic1 + exMC_chic2)