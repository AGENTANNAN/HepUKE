# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# J/psi at sqrt(s) = 3.097 GeV
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 3.097 GeV real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card: J/psi -> Lambda anti-Lambda; each Lambda decays 50/50 hadronically
# (p pi) and semileptonically (p mu nu); conjugate modes given for anti-Lambda
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    0.5000 p+ pi-             PHSP;
    0.5000 p+ mu- anti-nu_mu  PHSP;
    Enddecay

    Decay anti-Lambda0
    0.5000 anti-p- pi+        PHSP;
    0.5000 anti-p- mu+ nu_mu  PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for J/psi -> Lambda anti-Lambda
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_llbar_dt"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiLLbarDT"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                     # Charged-track selection (N_track = 4)
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"                   # exactly two positive tracks
    nChrn     "==2"                   # exactly two negative tracks
    nNet      "==0"
  }
  .pid(method: :probability) {        # PID: probability method, prob > 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # lepton (muon) vs electron
    identify :proton, against: [:kaon, :pion]   # protons against kaons/pions
    nprp ">=1"                        # at least one p
    nprm ">=1"                        # at least one anti-p
  }
  .secondary_vertex_fit([:prp, :pim]) {          # build Lambda from p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # build anti-Lambda from anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar]) {       # 4C fit on the Lambda anti-Lambda pair
    nominal
    constrain_four_momentum
    chi2_cut 200                                 # loose; optimal tight cut applied in ROOT
  }

# BOSS-side procedures without a dedicated DSL primitive, kept for the systematic step
my_algorithm
  .note(:tag_candidate_selection, "Tag-side Lambda chosen as the candidate minimising |DeltaE| within [-17, 13] MeV and satisfying M_BC in [1.089, 1.143] GeV/c^2, tag secondary-vertex chi2 < 100 and decay length > 2 sigma from the IP; these tag-side quantities are evaluated at DT-candidate selection.")
  .note(:muon_selection, "Signal muon required with L_mu > 0.001 and L_mu > L_e (either charge, charge-conjugate channels combined); handled by the high-momentum lepton identification of the PID step.")
  .note(:background_veto, "The fully-hadronic Lambda anti-Lambda hypothesis is vetoed by requiring the 4C kinematic-fit chi2 > 20: the semileptonic signal, with its undetected neutrino, does not satisfy the p pi / anti-p pi hypothesis. Only an upper chi2_cut is expressible in the fit block, so the lower-bound veto is applied on the stored chi2.")
  .note(:signal_recoil_mass, "Signal side additionally requires M_recoil(anti-Lambda p) > 0.170 GeV/c^2.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])