# Core DSL classes and dependencies are loaded automatically at execution.
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample at 3.097 GeV

# Decay card for the signal process J/psi -> gamma eta_c, eta_c -> Lambda anti-Lambda (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c        PHSP;
    Enddecay

    Decay eta_c
    1.0000 Lambda0 anti-Lambda0   PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-             HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal mode (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etac_llbar"
  config.related_dataset = jpsi_data        # associate with the 3.097 GeV real data
  config.events          = 100000           # number of events to generate
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaCLLbar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                     # Charged track selection
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        100.0                   # |Vz| < 100 cm
    Vr        10.0                    # Vr < 10 cm
    nChrp     ">=2"                   # At least 2 positively charged tracks
    nChrn     ">=2"                   # At least 2 negatively charged tracks
    nNet      "==0"                   # Net charge zero
  }
  .select_photon {                    # Photon selection
    tdc_emc_start     0               # EMC TDC start
    tdc_emc_end       14              # EMC TDC end
    angle_to_track    10.0            # Min angle to nearest charged track (degrees)
    energyThreshold_b 0.025           # 25 MeV threshold in the barrel region
    energyThreshold_e 0.050           # 50 MeV threshold in the endcap region
    nGam              ">=1"           # At least one photon
  }
  .pid(method: :probability) {        # PID by the probability method
    prob_cut 0.001                    # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]  # identify p+ and p- (charge-conjugation shorthand)
    nprp ">=1"                        # At least one proton
    nprm ">=1"                        # At least one anti-proton
  }
  .remove([:prp <= :chrgp])           # Remove identified protons from positive charged list
  .remove([:prm <= :chrgn])           # Remove identified anti-protons from negative charged list
  .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks -> pi+ / pi-
  .secondary_vertex_fit([:prp, :pim]) {   # Reconstruct Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {   # Reconstruct Lambda_bar -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {  # 4C kinematic fit to gamma Lambda anti-Lambda
    nominal                   # use corrected four-momenta from this fit
    constrain_four_momentum   # constrain total four-momentum to CMS energy
    chi2_cut 200              # loose chi-square cut; tight cut applied in ROOT
  }

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])