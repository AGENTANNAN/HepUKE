# Core DSL classes and dependencies are loaded automatically at execution
### Dataset preparation ###
# Real data at the major c.m. energy points of the 3.645–3.891 GeV scan.
data_3650 = DatasetManager.real_data.find("709_3650")   # 3.650 GeV
data_3682 = DatasetManager.real_data.find("709_3682")   # 3.682 GeV
data_3686 = DatasetManager.real_data.find("709_3686")   # 3.686 GeV (psi(2S))
data_3773 = DatasetManager.real_data.find("712_3773")   # 3.773 GeV (psi(3770))
data_3780 = DatasetManager.real_data.find("712_3780")   # 3.800 GeV

# Corresponding inclusive MC samples for each major energy point.
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")
incMC_3682 = DatasetManager.inclusive_mc.find("709_3682")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3780 = DatasetManager.inclusive_mc.find("712_3780")

# Exclusive MC decay cards at 3.686 GeV for the two J/psi leptonic modes.
# Signal generation: KKMC provides e+e- -> psi(3686)/psi(3770) with ISR (the top
# mother psi(2S) is the KKMC-level resonance), then EvtGen decays the psi(2S)
# producing the J/psi, which decays leptonically.
decay_card_ee = <<~DECAYCARD
    Decay psi(2S)
    1.0000 J/psi pi+ pi-          PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-                  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.0000 J/psi pi+ pi-          PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-                PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each of the two J/psi leptonic modes at 3.686 GeV.
exMC_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_inclJpsi_ee"
    config.related_dataset = data_3686
    config.events          = 500000
    config.decay_card      = decay_card_ee
    config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_inclJpsi_mumu"
    config.related_dataset = data_3686
    config.events          = 500000
    config.decay_card      = decay_card_mumu
    config.cross_section   = :default
end

# Save the signal-MC configurations for later use.
exMC_ee.save_to_config(format: :yaml, file_path: 'temp_for_test')
exMC_mumu.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "inclJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # c.m. energy in GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection for the two leptonic channels (l+ l- final state, l = e or mu).
event_selection = Selection.new
event_selection.select_track {          # Charged-track quality cuts
                  cos_theta   0.93      # |cos(theta)| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # |Vr| < 1 cm
                  nChrp       ">=1"     # At least one positively charged track
                  nChrn       ">=1"     # At least one negatively charged track
                }
               .select_photon {          # Photon selection
                  energyThreshold_b 0.025  # >= 25 MeV (barrel)
                  energyThreshold_e 0.050  # >= 50 MeV (endcap)
                  tdc_emc_start     0      # TDC window start
                  tdc_emc_end       14     # TDC window end
                  angle_to_track    10.0   # > 10 deg from nearest charged track
                  nGam              ">=1"  # At least one good photon required
                }
               .pid(method: :probability) {   # PID: high-momentum leptons (e / mu)
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.7  # p>1.0 GeV -> lepton; EMC eraw>0.7 -> electron, else muon
                  nlp   ">=1"   # At least one l+
                  nlm   ">=1"   # At least one l-
                }
# No conventional kinematic fit is performed for this analysis: the BOSS selection
# stops here and the J/psi yield is extracted from a fit to the l+l- invariant-mass
# spectrum in the ROOT stage.

# The full 69-point 3.645-3.891 GeV scan is assembled by combining the scan samples
# (not a single named dataset); the major points above are used as the anchors.
my_algorithm.note(:scan_sample_assembly,
  "the 69-point 3.645-3.891 GeV scan sample is assembled by merging the scan data
   files by run number; only the five major energy points (3.650, 3.682, 3.686,
   3.773, 3.800 GeV) are addressed as named datasets here")
# Additional quality cuts are applied in the ROOT stage (round-trip to the ROOT spec).
my_algorithm.note(:root_stage_quality_cuts,
  "no kinematic fit; J/psi yield from a fit to the l+l- invariant-mass spectrum.
   Extra quality cuts imposed in ROOT: e E/p>0.7, mu 0.05<E/p<0.35, lepton
   |cos(theta)|<0.81, opening angle <179 deg, lepton momentum in
   [1 GeV, 0.47*ECMS], pi/K separation from dE/dx+TOF (CL(pi)>CL(K)), and the
   event topology (exactly two tracks plus >=1 photon, or 3-4 charged tracks).
   Full-chain selection efficiency is 58.8-60.8%.")

# Generate the algorithm. One algorithm suffices for both channels because the
# combined e/mu lepton lists share one selection; the decay card provides the
# kinematic variables for the l+l- system.
my_algorithm.with_decay_card(decay_card_ee).apply(event_selection)

# Execute on all real data points, their inclusive MC, and both exclusive signal MCs.
root_files = my_algorithm.execute_on([
    data_3650, data_3682, data_3686, data_3773, data_3780,
    incMC_3650, incMC_3682, incMC_3686, incMC_3773, incMC_3780,
    exMC_ee, exMC_mumu
])