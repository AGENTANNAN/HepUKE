# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# XYZ + R-scan real data at 40 energy points, sqrt(s) = 3.7730 - 4.7008 GeV
scan_data = [
  DatasetManager.real_data.find("712_3773"),   # 3773 MeV (psi(3770))
  DatasetManager.real_data.find("703_3810"),   # 3807.7 MeV
  DatasetManager.real_data.find("703_3872"),   # 3867.4 MeV
  DatasetManager.real_data.find("703_3900"),   # 3896.2 MeV
  DatasetManager.real_data.find("703_4009"),   # 4007.6 MeV
  DatasetManager.real_data.find("703_4090"),   # 4085.5 MeV
  DatasetManager.real_data.find("705_4130"),   # 4128.5 MeV
  DatasetManager.real_data.find("705_4160"),   # 4157.4 MeV
  DatasetManager.real_data.find("703_4180"),   # 4178 MeV
  DatasetManager.real_data.find("703_4190"),   # 4188.8 MeV
  DatasetManager.real_data.find("703_4200"),   # 4198.9 MeV
  DatasetManager.real_data.find("703_4210"),   # 4209.2 MeV
  DatasetManager.real_data.find("703_4220"),   # 4218.7 MeV
  DatasetManager.real_data.find("703_4230"),   # 4226.3 MeV
  DatasetManager.real_data.find("703_4237"),   # 4235.7 MeV
  DatasetManager.real_data.find("703_4245"),   # 4241.7 MeV
  DatasetManager.real_data.find("703_4246"),   # 4243.8 MeV
  DatasetManager.real_data.find("703_4260"),   # 4258.0 MeV
  DatasetManager.real_data.find("703_4270"),   # 4266.8 MeV
  DatasetManager.real_data.find("703_4280"),   # 4277.7 MeV
  DatasetManager.real_data.find("705_4290"),   # 4287.9 MeV
  DatasetManager.real_data.find("703_4310"),   # 4307.9 MeV
  DatasetManager.real_data.find("705_4315"),   # 4312.1 MeV
  DatasetManager.real_data.find("705_4340"),   # 4337.4 MeV
  DatasetManager.real_data.find("703_4360"),   # 4358.3 MeV
  DatasetManager.real_data.find("705_4380"),   # 4377.4 MeV
  DatasetManager.real_data.find("703_4390"),   # 4387.4 MeV
  DatasetManager.real_data.find("705_4400"),   # 4396.5 MeV
  DatasetManager.real_data.find("703_4420"),   # 4415.6 MeV
  DatasetManager.real_data.find("705_4440"),   # 4436.2 MeV
  DatasetManager.real_data.find("703_4470"),   # 4467.1 MeV
  DatasetManager.real_data.find("703_4530"),   # 4527.1 MeV
  DatasetManager.real_data.find("703_4575"),   # 4574.5 MeV
  DatasetManager.real_data.find("703_4600"),   # 4599.5 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.9 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.0 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.9 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.2 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.9 MeV
  DatasetManager.real_data.find("706_4700")    # 4698.8 MeV
]

# Matching inclusive MC samples (available at the corresponding energy points)
scan_incMC = [
  DatasetManager.inclusive_mc.find("712_3773"),
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

# Decay card (mode I): e+e- -> pi+ pi- J/psi, J/psi -> e+ e-
# KKMC-generated continuum production, psi(4260) as BESIII top-mother convention.
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card (mode II): e+e- -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  mu+  mu-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Signal MC generated at every energy point of the scan (one ExclusiveMC per point)
exMC_scan_ee = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_scan_pipijpsi_ee"   # becomes ..._<boss>_<energy> per point
  config.events        = 100000
  config.decay_card    = decay_card_signal_ee
  config.cross_section = :default
end

exMC_scan_mumu = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_scan_pipijpsi_mumu"
  config.events        = 100000
  config.decay_card    = decay_card_signal_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsi"                              # pi+ pi- J/psi
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.237]})   # nominal central scan energy (MeV->GeV); per-point value from the dataset
            .set_alias({"std::vector<double>" => "Vdouble"})

# Both J/psi decay modes share the identical final-state topology (pi+ pi- l+ l-)
# and the identical selection, so a single Algorithm / Selection chain is used.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93     # |cos(theta)| < 0.93
                  Vz        10.0     # |Vz| < 10 cm
                  Vr        1.0      # Vr < 1 cm
                  nChrp     "==2"    # exactly two positive charged tracks
                  nChrn     "==2"    # exactly two negative charged tracks
                  nNet      "==0"    # net charge zero  -> 4 charged tracks in total
                }
               .pid(method: :probability) {
                  # High-momentum tracks (p > 1.0 GeV) treated as leptons;
                  # an electron if the EMC energy > 0.6 GeV, otherwise a muon.
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  # Pions identified against kaons and protons (both charges at once)
                  identify :pion, against: [:kaon, :proton]
                  npip "==1"       # one pi+
                  npim "==1"       # one pi-
                  nlp  "==1"       # one l+
                  nlm  "==1"       # one l-
                }
               # 4C kinematic fit to pi+ pi- l+ l- , nominal fit (corrected four-momenta saved)
               .kinematic_fit([:pip, :pim, :lp, :lm]) do
                  nominal
                  constrain_four_momentum            # 4C energy-momentum conservation
                  chi2_cut 60                        # chi^2 < 60
                }

# Iterative KKMC-based ISR correction (no formal DSL construct; preserved for the
# systematic-uncertainties / cross-section extraction stage)
my_algorithm
  .note(:isr_correction, "iterative KKMC-based ISR correction applied to extract the Born
    cross section: the ISR radiator from KKMC is folded with the measured pi+pi-J/psi line
    shape and iterated until the (1+delta) correction factor is self-consistent")

# Generate the complete algorithm for the signal process in the decay card
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)

# Further requirements (cos(pi+,pi-)<0.98, cos(pi+-,e-+)<0.98, pion dE/dx discriminators,
# BDT against two-photon background, EMC electron-muon separation, inverse-variance mode
# combination) act on kinematic-fit-corrected variables and are applied at the ROOT stage.

# Execute on real data, inclusive MC and signal MC for both J/psi decay modes
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMC_scan_ee + exMC_scan_mumu)