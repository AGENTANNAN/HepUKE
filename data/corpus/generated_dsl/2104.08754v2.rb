# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# e+e- -> phi Lambda anti-Lambda measured at 28 energy points from 3.51 to 4.60 GeV
# (BOSS 703/705 real data, total 19.5 fb-1).
data_points = %w[
  703_3810 703_3872 703_3900 703_4009 703_4090
  703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4245 703_4246 703_4260
  703_4270 703_4280 703_4310 703_4360 703_4390
  703_4420 703_4470 703_4530 703_4575 703_4600
  705_4130 705_4160 705_4290
].map { |name| DatasetManager.real_data.find(name) }

data_4180  = DatasetManager.real_data.find("703_4180")     # 4.178 GeV reference point
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")  # inclusive MC at 4.178 GeV

# Decay card for the signal process e+e- -> phi Lambda anti-Lambda (EvtGen format).
# RecID mapping: 0=psi(4260) 1=phi 2=Lambda0 3=anti-Lambda0 4=K+ 5=K- 6=p+ 7=pi- 8=anti-p- 9=pi+
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0 phi Lambda0 anti-Lambda0 PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay Lambda0
  1.0 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# Exclusive MC for e+e- -> phi Lambda anti-Lambda at the 4.178 GeV point
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_phiLambdaLambdabar"
  config.related_dataset = data_4180
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PhiLambdaLambdabar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that cannot be expressed in the formal DSL
my_algorithm
  .note(:secondary_vertex_chi2, "Lambda -> p+ pi- secondary vertex fit required to satisfy chi2 < 100")
  .note(:flight_distance_significance, "require flight-distance significance > 2 sigma between the interaction point and the Lambda decay vertex")
  .note(:partial_rec_chi2, "partial reconstruction of the Lambda K+ K- system with the missing mass constrained to the nominal Lambda mass, chi2 < 30")
  .note(:candidate_selection, "for multiple candidates, choose the combination minimising the combined secondary-vertex plus kinematic chi2")

# Event selection chain
event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vr        10.0    # Vr < 10 cm
    Vz        20.0    # |Vz| < 20 cm
    nChrp     "==2"   # exactly two positively charged tracks
    nChrn     "==2"   # exactly two negatively charged tracks
    nNet      "==0"   # net charge zero
  }
  # no photon requirement -> no select_photon step
  .pid(method: :probability) {
    prob_cut 0.001                                   # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]        # p+ / anti-p against pi and K
    identify :kaon,   against: [:proton, :pion]      # K+ / K- against p and pi
    identify :pion,   against: [:kaon, :proton]      # pi+ / pi- against K and p
    nprp "==1"    # one p+  (Lambda -> p+ pi-)
    npim "==1"    # one pi- (Lambda -> p+ pi-)
    nkp  "==1"    # one K+  (phi -> K+ K-)
    nkm  "==1"    # one K-  (phi -> K+ K-)
  }
  # Reconstruct Lambda -> p+ pi- with a secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Treat the anti-Lambda as missing and partially reconstruct the Lambda K+ K- system
  .partial_miss([3]) {                            # recID 3 = anti-Lambda0 (daughters expanded automatically)
    best_combination_by_mass :Lambda, 1.1157      # Lambda mass window [1.112, 1.120]
    best_combination_by_mass :phi,    1.0195      # phi mass window [1.01, 1.03]
    require_recoil_mass 1.112, 1.120              # missing (anti-Lambda) mass constrained to nominal Lambda mass
  }

# Generate the algorithm for the decay card and run it on data + MC
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + [incMC_4180, exMC_signal])