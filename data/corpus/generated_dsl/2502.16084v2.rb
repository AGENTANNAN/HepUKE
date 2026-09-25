# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Five 713 R-scan energy points for the single-inclusive π± / K± continuum measurement
data_2000 = DatasetManager.real_data.find("713_2000")   # √s = 2.000 GeV
data_2200 = DatasetManager.real_data.find("713_2200")   # √s = 2.200 GeV
data_2396 = DatasetManager.real_data.find("713_2396")   # √s = 2.396 GeV
data_2644 = DatasetManager.real_data.find("713_2644")   # √s = 2.644 GeV
data_2900 = DatasetManager.real_data.find("713_2900")   # √s = 2.900 GeV
rscan_points = [data_2000, data_2200, data_2396, data_2644, data_2900]

# Matching inclusive MC samples for the same five points
incMC_2000 = DatasetManager.inclusive_mc.find("713_2000")
incMC_2200 = DatasetManager.inclusive_mc.find("713_2200")
incMC_2396 = DatasetManager.inclusive_mc.find("713_2396")
incMC_2644 = DatasetManager.inclusive_mc.find("713_2644")
incMC_2900 = DatasetManager.inclusive_mc.find("713_2900")
rscan_incMC = [incMC_2000, incMC_2200, incMC_2396, incMC_2644, incMC_2900]

# ConExc continuum decay card — R-value mode 74110 (light hadrons).
# ConExc is auto-detected from the literal token "ConExc": the no-KKMC simulation
# template is used and "Particle vpho <ECMS> 0.0" is injected per energy point,
# so no explicit Particle vpho line is written here.
continuum_decay_card = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 74110;
    Enddecay
    End
DECAYCARD

# One 1M-event exclusive continuum MC per R-scan point (same card, different √s)
exMCs_continuum = DatasetManager.create_exclusive_mc_for(rscan_points) do |config|
  config.sample_name   = "continuum_piK_rvalue"   # per-point suffix appended automatically
  config.events        = 1_000_000
  config.decay_card    = continuum_decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "ContinuumPiK"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.900]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:ecms_scan, "ECMS differs at each of the five R-scan points (2.000, 2.200, 2.396, 2.644, 2.900 GeV); the set_constant value is a placeholder for the multi-energy scan")
            .note(:background_veto, "Hadronic event selection (Bhabha and γγ rejection, prong counting, QED background suppression) is performed in ROOT, not in BOSS")

# Inclusive track / photon selection; π± and K± are extracted inclusively from the
# selected hadronic events, so no decay channel is reconstructed and no kinematic fit is applied.
event_selection = Selection.new
event_selection.select_track {            # Charged track quality cuts
                  cos_theta 0.93          # |cosθ| < 0.93
                  Vz        10.0          # |Vz| < 10 cm
                  Vr        1.0           # Vr < 1 cm
                }
               .select_photon {           # Photon selection
                  tdc_emc_start     0     # EMC timing window [0, 14]
                  tdc_emc_end       14
                  energyThreshold_b 0.025 # 25 MeV in the barrel
                  energyThreshold_e 0.050 # 50 MeV in the endcap
                  angle_to_track    10.0  # ≥ 10° from any charged track
                }
               .pid(method: :probability) {  # combined dE/dx + TOF probability PID
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]  # π± separated from K and p
                  identify :kaon, against: [:pion, :proton]  # K± separated from π and p
                }

my_algorithm.with_decay_card(continuum_decay_card).apply(event_selection)
root_files = my_algorithm.execute_on(rscan_points + rscan_incMC + exMCs_continuum)