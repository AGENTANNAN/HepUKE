# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# Decay card: psi(2S) -> gamma chi_c0, chi_c0 -> Lambda anti-Lambda omega
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(2S) -> gamma chi_c1, chi_c1 -> Lambda anti-Lambda omega
decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(2S) -> gamma chi_c2, chi_c2 -> Lambda anti-Lambda omega
decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (200k events each) for the three chi_cJ signals
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c0_LLbar_omega"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c1_LLbar_omega"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c2_LLbar_omega"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# The three chi_cJ states share an identical final state
# (gamma Lambda anti-Lambda pi+ pi- pi0) and identical selection criteria,
# so a single Algorithm / Selection chain is used for all of them.
alg_name = "GammaChiCJLambdaLambdaBarOmega"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) center-of-mass energy (GeV)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                 # charged-track selection
    cos_theta   0.93              # |cos(theta)| < 0.93
    Vz          10.0              # |Vz| < 10 cm
    Vr          10.0              # Vr < 1 cm (10 mm) in the transverse plane
    nChrp       ">=3"             # at least 3 positively charged tracks
    nChrn       ">=3"             # at least 3 negatively charged tracks
    nNet        "==0"             # net charge zero
  }
  .select_photon {                # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0        # minimum angle to nearest charged track (degrees)
    energyThreshold_b 0.025       # 25 MeV barrel threshold
    energyThreshold_e 0.050       # 50 MeV endcap threshold
    nGam              ">=3"       # at least 3 photons
  }
  .pid(method: :probability) {    # probability-based PID
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]  # identify p and p-bar vs K and pi
    nprp ">=1"                    # at least one proton
    nprm ">=1"                    # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # drop identified (anti-)protons from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # treat remaining charged tracks as pions
  .secondary_vertex_fit([:prp, :pim]) {       # reconstruct Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # reconstruct anti-Lambda -> p-bar pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct pi0 -> gamma gamma (1C mass fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                # at least one pi0 candidate
  }
  # Nominal 5C fit: gamma Lambda anti-Lambda pi+ pi- pi0 with 4-momentum conservation
  # and the pi0 mass constraint (pi0 enters as a fitted particle of nominal mass).
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 30
  }
  # Competing hypothesis psi(2S) -> Lambda anti-Lambda omega (no radiative photon):
  # no chi2_cut and no nominal -> only the chi2_5c/4c value is stored for a ROOT-level veto.
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :pi0]) {
    constrain_four_momentum
  }

# Attach the (shared-topology) decay card and render the algorithm
my_algorithm.with_decay_card(decay_card_chi_c1).apply(event_selection)

# Execute on real data, inclusive MC and the three signal exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_chi_c0, exMC_chi_c1, exMC_chi_c2])