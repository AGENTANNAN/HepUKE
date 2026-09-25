# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset preparation ###
# Five identifiable R-scan points (2950, 2981, 3000, 3020, 3080 MeV) of the 26-point φη scan
data_points = [
  DatasetManager.real_data.find("713_2950"),
  DatasetManager.real_data.find("713_2981"),
  DatasetManager.real_data.find("713_3000"),
  DatasetManager.real_data.find("713_3020"),
  DatasetManager.real_data.find("713_3080")
]

# Matching inclusive MC samples for the five points
incMC_points = [
  DatasetManager.inclusive_mc.find("713_2950"),
  DatasetManager.inclusive_mc.find("713_2981"),
  DatasetManager.inclusive_mc.find("713_3000"),
  DatasetManager.inclusive_mc.find("713_3020"),
  DatasetManager.inclusive_mc.find("713_3080")
]

# ConExc decay card for e+e- -> phi eta (ConExc mode 23).
# The literal token `ConExc` switches the DSL to the no-KKMC template and injects
# `Particle vpho <ECMS>` per energy point, so no `Particle vpho` line is written here.
decay_card_phi_eta = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 23;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 50k-event exclusive signal MC per scan point (same card/cross section, per-point energies)
exMCs_phi_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phieta_conexc"
  config.events        = 50_000
  config.decay_card    = decay_card_phi_eta
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PhiEta"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.080]})  # per-point ECMS follows each dataset in the scan
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        100.0    # |Vz| < 100 cm
      Vr        10.0     # Vr < 10
      nChrp     ">=1"    # at least one positively charged track
      nChrn     ">=1"    # at least one negatively charged track
      nNet      "==0"    # net charge zero
    }
    .select_photon {
      tdc_emc_start     0      # EMC time window 0-700 ns
      tdc_emc_end       14
      angle_to_track    10.0   # angle to nearest charged track > 10 degrees
      energyThreshold_b 0.025  # barrel energy threshold 25 MeV
      energyThreshold_e 0.050  # endcap energy threshold 50 MeV
      nGam              ">=2"  # at least two photons (eta -> gamma gamma)
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]  # K+ and K- vs pions/protons
      nkp ">=1"                                  # at least one K+
      nkm ">=1"                                  # at least one K-
    }
    .kinematic_fit([:kp, :km, :gamma, :gamma]) {
      nominal                  # nominal 4C fit
      constrain_four_momentum  # four-momentum constraint
      chi2_cut 200             # loose BOSS cut; the tighter chi2 < 85 is applied in ROOT
    }
# The eta and phi mass windows (|M(gg)-M_eta| < 30 MeV, phi window) and the paper's
# tighter chi2 < 85 are post-fit cuts applied in the ROOT analysis -> out of BOSS scope.

alg.with_decay_card(decay_card_phi_eta).apply(event_selection)

root_files = alg.execute_on(data_points + incMC_points + exMCs_phi_eta)