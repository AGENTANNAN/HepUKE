# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# e+e- -> Sigma0 anti-Sigma0 studied over the BESIII scan energies 3.650 - 4.946 GeV.
scan_samples = %w[
  709_3650
  703_3810
  703_3872
  703_3900
  703_4009
  703_4090
  705_4130
  705_4160
  703_4180
  703_4190
  703_4200
  703_4210
  703_4220
  703_4230
  703_4237
  703_4245
  703_4246
  703_4260
  703_4270
  703_4280
  705_4290
  703_4310
  705_4315
  705_4340
  703_4360
  705_4380
  703_4390
  705_4400
  703_4420
  705_4440
  703_4470
  703_4530
  703_4575
  703_4600
  706_4610
  706_4620
  706_4640
  706_4660
  706_4680
  706_4700
  707_4740
  707_4750
  707_4780
  707_4840
  707_4914
  707_4946
]

# Real data at each scan point
scan_data  = scan_samples.map { |s| DatasetManager.real_data.find(s) }
# Matching inclusive MC for background / continuum study
scan_incMC = scan_samples.map { |s| DatasetManager.inclusive_mc.find(s) }

# ConExc decay card (continuum R-scan / Born-cross-section measurement).
# The DSL auto-detects the literal "ConExc" token, switches to the no-KKMC
# simulation template and injects "Particle vpho <ECMS> 0.0" per energy point,
# so Particle vpho must NOT be written here. Mode 3 = Sigma0 anti-Sigma0.
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1 ConExc 3;
  Enddecay

  Decay Sigma0
  1.0 Lambda0 gamma PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0 anti-Lambda0 gamma PHSP;
  Enddecay

  Decay Lambda0
  1.0 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC at every scan point (one MC per energy point,
# sharing decay card / cross section / event count)
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_sigma0_sigma0bar_conexc"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Sigma0Sigma0bar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.650]})  # representative value; per-point ECMS injected with the ConExc card

# BOSS-side procedure with no DSL counterpart: ISR / vacuum-polarization factors
# are read from the ConExc generator log for each energy point.
my_algorithm.note(:isr_vp_correction,
  "ISR and vacuum-polarization correction factors are extracted from the ConExc
   generator log per scan point and applied when forming the Born cross section")

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz        10.0        # |Vz| < 10 cm
    Vr        1.0         # Vr < 1 cm
    nChrp     ">=2"       # at least two positively charged tracks
    nChrn     ">=2"       # at least two negatively charged tracks
  }
  .select_photon {
    tdc_emc_start     0      # EMC TDC window 0 - 14
    tdc_emc_end       14
    angle_to_track    10.0   # angle to nearest charged track > 10 degrees
    energyThreshold_b 0.025  # barrel energy > 25 MeV
    energyThreshold_e 0.050  # endcap energy > 50 MeV
    nGam              ">=2"  # at least two photons
  }
  # 1st PID pass: locate (anti-)protons, then remove them from the charged lists
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p (charge-conjugation shorthand)
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # 2nd PID pass: pi+/pi- from the remaining tracks
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # pi+ and pi-
    npip ">=1"
    npim ">=1"
  }
  # Lambda -> p pi- via secondary vertex fit (mass-difference minimisation)
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # anti-Lambda -> anti-p pi+ via secondary vertex fit
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Final 4C kinematic fit of Sigma0 anti-Sigma0 -> Lambda anti-Lambda gamma gamma
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMC_signal)