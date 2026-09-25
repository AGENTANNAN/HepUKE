# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
data_4260  = DatasetManager.real_data.find("703_4260")      # psi(4260) data @ 4.260 GeV (4.260 GeV point of the 3.90-4.60 GeV scan)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")   # corresponding inclusive MC sample

# Decay card — Mode I: e+e- -> K_S0 K+ pi- pi0
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K+ pi- pi0  PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode II: e+e- -> K_S0 K+ pi- eta
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K+ pi- eta  PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC sample for each of the two modes
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_KSKPiPi0"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_KSKPiEta"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# K_S0 criteria applied inside the secondary-vertex fit that the DSL cannot express formally
ks0_note = "K_S0 candidates required |M(pi+pi-) - M_K_S0| < 12 MeV/c^2 and a decay-length significance " \
           "> 2 sigma, with the best K_S0 taken as the candidate having the smallest secondary-vertex " \
           "chi^2. These SV-fit criteria are not expressible in secondary_vertex_fit " \
           "(by_minimizing_mass_difference is used to build the K_S0)."

alg_name_modeI = "KSKPiPi0"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .note(:ks0_selection, ks0_note)

alg_name_modeII = "KSKPiEta"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:ks0_selection, ks0_note)

# Shared event-selection chain (identical for both decay modes)
event_selection = Selection.new
    .select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        10.0      # |Vz| < 10 cm
        Vr        1.0       # Vr < 1 cm
        nChrp     ">=2"     # at least two positive tracks
        nChrn     ">=2"     # at least two negative tracks
        nNet      "==0"     # net charge zero
    }
    .select_photon {
        tdc_emc_start     0        # TDC window 0-14
        tdc_emc_end       14
        energyThreshold_b 0.025    # 25 MeV barrel
        energyThreshold_e 0.050    # 50 MeV endcap
        angle_to_track    20.0     # at least 20 deg from any charged track
        nGam              ">=2"    # at least two photons
    }
    .pid(method: :probability) {
        prob_cut 0.001                               # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]    # kaons vs pions/protons
        nkp "==1"                                    # one K+
        nkm "==1"                                    # one K-
    }
    .remove([:kp <= :chrgp, :km <= :chrgn])          # drop identified kaons from the generic charged lists
    .assign({:chrgp => :pip, :chrgn => :pim})        # remaining positives -> pi+, negatives -> pi-
    .secondary_vertex_fit([:pip, :pim]) {            # secondary vertex fit pi+ pi- -> K_S0
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:K_S0, :kp, :pim, :gamma, :gamma]) {   # nominal 4C fit to K_S0 K+ pi- gamma gamma
        nominal
        constrain_four_momentum
        chi2_cut 60
    }

# Both decay modes share the same selection chain (one Selection, applied to both algorithms)
alg_modeI.with_decay_card(decay_card_modeI).apply(event_selection)
alg_modeII.with_decay_card(decay_card_modeII).apply(event_selection)

# Post-fit K_S0 / pi0 / eta mass windows and the final χ² optimisation are applied in the ROOT analysis stage
root_files_modeI  = alg_modeI.execute_on([data_4260, incMC_4260, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([data_4260, incMC_4260, exMC_modeII])