# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Seven c.m. energy points of the 4.600-4.699 GeV scan (total ~4.4 fb^-1)
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV, 586.9 pb^-1
data_4610 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV
data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

# Matching inclusive MC samples at the same seven energy points
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c-
#   single-tag channel Lambda_c+ -> Lambda pi+ pi0, Lambda -> p pi-, pi0 -> gamma gamma
#   other side anti-Lambda_c- -> anti-p K+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-  PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Lambda0 pi+ pi0  PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-  PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal-MC events generated at every energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_LcLcbar"
  config.events        = 500000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "LcST"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.660]})   # representative c.m. energy; ECMS set per point at run time
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # charged track quality cuts
                  cos_theta  0.93       # |cos(theta)| < 0.93
                  Vz         10.0       # |Vz| < 10 cm
                  Vr         1.0        # Vr < 1 cm
                }
               .select_photon {         # photon selection
                  tdc_emc_start     0   # EMC time window 0 - 700 ns
                  tdc_emc_end       14
                  energyThreshold_b 0.025   # E > 25 MeV (barrel)
                  energyThreshold_e 0.050   # (endcap)
                  angle_to_track    10.0    # >= 10 deg from any charged track
                  nGam              ">=2"   # at least two photons to form pi0
                }
               .pid(method: :probability) {   # probability-method PID
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]   # L(p) > L(K), L(p) > L(pi)
                  identify :pion,   against: [:kaon]          # pi vs K
                  nprp ">=1"
                  npip ">=1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit -> pi0
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)  # gamma-gamma mass window
                  chi2_cut 25
                  npi0     ">=1"
                }
               .secondary_vertex_fit([:prp, :pim]) {   # secondary vertex -> Lambda
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Final nominal fit: form Lambda_c+ from Lambda pi+ pi0 (3C intent)
               .kinematic_fit([:Lambda, :pip, :pi0]) {
                  nominal
                  invariant_mass_of(:Lambda).within(1.111, 1.121)                               # Lambda mass window
                  invariant_mass_of(:Lambda).constrain_to_nominal_mass_of(:Lambda)              # m(Lambda) constraint
                  invariant_mass_of(:Lambda, :pip, :pi0).constrain_to_nominal_mass_of(:"Lambda_c+")  # m(Lambda_c+) constraint
                  constrain_four_momentum                                                        # beam / recoil-mass constraint
                  chi2_cut 200
                }

# Inexpressible BOSS-side procedures preserved for downstream use
my_algorithm
  .note(:background_veto, "Sigma0 background rejected by vetoing M(Lambda gamma) in [1.179, 1.203] GeV")
  .note(:lambda_vertex_selection, "secondary-vertex fit for Lambda -> p pi- requires chi2 < 100 and decay length > 2 sigma")
  .note(:delta_e_selection, "Lambda_c+ candidates required to satisfy DeltaE in (-0.03, 0.02) GeV; best candidate per event chosen by minimum |DeltaE|")
  .note(:mbc_signal_region, "M_BC signal region determined separately for each of the seven c.m. energy points")
  .note(:ecms_per_point, "algorithm runs over seven c.m. energy points (4.600-4.699 GeV); ECMS is set per point")
  .note(:kinematic_fit_constraints, "intended 3C fit constraining m(Lambda), m(Lambda_c+) and the recoil mass to m(Lambda_c+); implemented via the Lambda/Lambda_c+ mass constraints plus the beam four-momentum constraint")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points + exMC_signal)