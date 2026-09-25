# =====================================================================
# BOSS / dataset-preparation DSL
# Born cross section for e+e- -> f1(1285) pi+ pi- via ISR,
# f1(1285) -> pi+ pi- eta, eta -> gamma gamma,
# at 45 CM energies (3.808 - 4.951 GeV); 7 representative scan points.
# =====================================================================

### Dataset description ###
# Real data at the representative scan points
data_3810 = DatasetManager.real_data.find("703_3810")   # 3.810 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.180 GeV
data_4260 = DatasetManager.real_data.find("703_4260")   # 4.260 GeV
data_4420 = DatasetManager.real_data.find("703_4420")   # 4.420 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.680 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # 4.946 GeV

# Corresponding inclusive MC
incMC_3810 = DatasetManager.inclusive_mc.find("703_3810")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

scan_points = [data_3810, data_4180, data_4260, data_4420, data_4600, data_4680, data_4946]
scan_incMC  = [incMC_3810, incMC_4180, incMC_4260, incMC_4420, incMC_4600, incMC_4680, incMC_4946]

# ConExc decay card for the ISR continuum signal e+e- -> f1(1285) pi+ pi-.
# The literal `ConExc` token selects the no-KKMC ISR template; the DSL injects
# `Particle vpho <ECMS> 0.0` per energy point (no explicit vpho line required).
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1.000 ConExc f1(1285) pi+ pi-;
    Enddecay

    Decay f1(1285)
    1.000 pi+ pi- eta;
    Enddecay

    Decay eta
    1.000 gamma gamma;
    Enddecay

    End
DECAYCARD

# One 100k-event ConExc exclusive MC per representative scan point
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_conexc_f1pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name     = "F1PiPiISR"
my_algorithm = Algorithm.new(alg_name)
my_algorithm
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({"ECMS" => [:double, 4.260]})   # representative; see per_point_ecms note
  .set_alias({"std::vector<double>" => "Vdouble"})
  .note(:per_point_ecms, "ISR continuum scan over 45 CM energies (3.808-4.951 GeV); the ECMS constant is only representative and the generator / kinematic fit must use the measured per-run CMS energy of each dataset.")
  .note(:efficiency_curve, "ISR and vacuum-polarisation correction factors are read from the ConExc generator log and applied per energy point.")

event_selection = Selection.new
event_selection
  .select_track {                 # exactly four charged tracks, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                # good photons; eta -> gamma gamma needs >= 2
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  .pid(method: :probability) {    # probability PID separating pions from kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]   # pi+ and pi- (charge-conjugation shorthand)
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                   # eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).within(0.5229, 0.5729)  # |M(gg)-M(eta)| < 25 MeV
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim, :eta]) {            # 4C fit (best chi2 combination auto-selected)
    nominal
    constrain_four_momentum
    chi2_cut 200   # loose at BOSS level; tightened to < 100 in the final (ROOT) analysis
  }

my_algorithm.with_decay_card(decay_card_conexc).apply(event_selection)
root_files = my_algorithm.execute_on(scan_points + scan_incMC + exMCs_signal)