# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data at the six center-of-mass energies of the 4.60-4.70 GeV scan (BOSS 706)
data_4610 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV

# Matching inclusive MC sample at each energy point
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

# Decay card for the signal process e+e- -> Lambda_c+ K-,
# Lambda_c+ -> Xi0 K+, Xi0 -> Lambda pi0, Lambda -> p pi-, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ K- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Xi0 K+ PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500,000-event exclusive MC for the signal mode, one sample per energy point
data_points = [data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lc_to_xi0k"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "LcToXi0K"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.66]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:ecms_scan, "analysis spans six energy points (4611.86-4698.82 MeV); the shared BOSS algorithm declares a single representative ECMS constant")

event_selection = Selection.new
    .select_track {          # charged track selection
        cos_theta 0.93       # |cos(theta)| < 0.93
        Vz        10.0       # |Vz| < 10 cm
        Vr        1.0        # Vr < 1 cm
        nChrp     ">=2"      # at least four charged tracks ...
        nChrn     ">=2"      # ... (2 positive and 2 negative)
        nNet      "==0"      # net charge zero
    }
    .select_photon {         # photon selection
        tdc_emc_start     0      # EMC TDC start
        tdc_emc_end       14     # EMC TDC end
        angle_to_track    10.0   # min angle to nearest charged track (deg)
        energyThreshold_b 0.025  # barrel energy threshold (GeV)
        energyThreshold_e 0.050  # endcap energy threshold (GeV)
        nGam              ">=2"  # at least two photons
    }
    .pid(method: :probability) {                    # probability-based PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # one proton (against K/pi)
        identify :kaon,   against: [:pion]          # one K+ (against pi)
        identify :pion,   against: [:kaon]          # at least one pi- (against K)
        nprp ">=1"
        nkp  ">=1"
        npim ">=1"
    }
    # remove overlaps between proton/charged-particle and kaon/charged-particle hypotheses
    .remove([:prp <= :chrgp, :kp <= :chrgp, :km <= :chrgn])
    # treat the remaining charged tracks as pions for the Lambda daughters
    .assign({:chrgp => :pip, :chrgn => :pim})
    # secondary vertex fit: Lambda -> p pi-, chosen by minimal mass difference
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # 1C Kalman fit: pi0 -> gamma gamma, invariant mass constrained to nominal pi0 mass
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"          # at least one pi0 candidate
    }
    # build Xi0 from Lambda pi0 and Lambda_c+ from Xi0 K+;
    # 4C fit with Xi0 mass consistency (Lambda mass already fixed by the SV fit)
    .kinematic_fit([:Lambda, :pi0, :kp, :km]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:Lambda, :pi0).constrain_to_nominal_mass_of(:Xi0)
        chi2_cut 200
    }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# The final DeltaE and MBC windows are applied in the ROOT analysis, not in this BOSS selection chain.
root_files = my_algorithm.execute_on(
  [data_4610, data_4620, data_4640, data_4660, data_4680, data_4700,
   incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700] + exMCs_signal
)