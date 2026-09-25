# =============================================================================
# ψ(3770) single-tag D branching-fraction measurement
#   Mode I   : D+ → K_S0 K_S0 K+
#   Mode II  : D+ → K_S0 K_S0 π+   (90% via K*(892)+ → K_S0 π+, 10% direct 3-body)
#   Mode III : D0 → K_S0 K_S0
#   Mode IV  : D0 → K_S0 K_S0 K_S0
# with K_S0 → π+π-
# =============================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 2.93 fb^-1 ψ(3770) data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

### Decay cards (EvtGen format) ###
# Mode I : D+ → K_S0 K_S0 K+
decay_card_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K_S0 K_S0 K+ PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II : D+ → K_S0 K_S0 π+  (0.9 via K*(892)+, 0.1 direct 3-body)
decay_card_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    0.9000 K_S0 K*(892)+ PHSP;
    0.1000 K_S0 K_S0 pi+ PHSP;
    Enddecay

    Decay K*(892)+
    1.0000 K_S0 pi+ VSS;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode III : D0 → K_S0 K_S0
decay_card_modeIII = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV : D0 → K_S0 K_S0 K_S0
decay_card_modeIV = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC — 300k events for each of the four modes ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_KsKsK"
  config.related_dataset = data_3773
  config.events          = 300000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_KsKsPi"
  config.related_dataset = data_3773
  config.events          = 300000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0_KsKs"
  config.related_dataset = data_3773
  config.events          = 300000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0_KsKsKs"
  config.related_dataset = data_3773
  config.events          = 300000
  config.decay_card      = decay_card_modeIV
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ----------------- Mode I : D+ → K_S0 K_S0 K+ -----------------
alg_modeI = Algorithm.new("DpToKsKsK")
alg_modeI.set_header(["DpToKsKsKAlg/DpToKsKsK.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:ks0_selection,
               "each K_S0 candidate required to satisfy |M(pi+pi-) - M(K_S0)| < 12 MeV/c^2 and "
               "flight-length significance L/sigma_L > 2; K_S0 daughter tracks use a looser "
               "|Vz| < 20 cm than the prompt-track |Vz| < 10 cm cut")

sel_modeI = Selection.new
sel_modeI.select_track {                       # charged track selection
             cos_theta 0.93
             Vz 10.0
             Vr 1.0
             nChrp "==3"                       # 2 K_S0 daughters (pi+) + 1 prompt K+
             nChrn "==2"                       # 2 K_S0 daughters (pi-)
             nNet  "==1"                       # D+ net charge
           }
         .pid(method: :probability) {          # prompt K+ identified against pions (dE/dx & TOF, CL_K vs CL_pi)
             prob_cut 0.001
             identify :kaon, against: [:pion]
             nkp "==1"                         # one prompt K+
           }
         .remove([:kp <= :chrgp])              # keep the prompt K+ out of the K_S0 daughter pool
         .assign({:chrgp => :pip, :chrgn => :pim})   # K_S0 daughters used without PID
         .secondary_vertex_fit([:pip, :pim]) {       # first K_S0, constrained to common vertex
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
         .secondary_vertex_fit([:pip, :pim]) {       # second K_S0
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
         .kinematic_fit([:K_S0, :K_S0, :kp]) {       # 4C fit to the D+ candidate
             nominal
             constrain_four_momentum
             chi2_cut 200
           }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_3773, incMC_3773, exMC_modeI])

# ----------------- Mode II : D+ → K_S0 K_S0 π+ -----------------
alg_modeII = Algorithm.new("DpToKsKsPi")
alg_modeII.set_header(["DpToKsKsPiAlg/DpToKsKsPi.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:ks0_selection,
                "each K_S0 candidate required to satisfy |M(pi+pi-) - M(K_S0)| < 12 MeV/c^2 and "
                "flight-length significance L/sigma_L > 2; K_S0 daughter tracks use a looser "
                "|Vz| < 20 cm than the prompt-track |Vz| < 10 cm cut")

sel_modeII = Selection.new
sel_modeII.select_track {                      # charged track selection
              cos_theta 0.93
              Vz 10.0
              Vr 1.0
              nChrp "==3"                      # 2 K_S0 daughters (pi+) + 1 prompt pi+
              nChrn "==2"                      # 2 K_S0 daughters (pi-)
              nNet  "==1"                      # D+ net charge
            }
          .pid(method: :probability) {         # prompt π+ identified against kaons (CL_pi vs CL_K)
              prob_cut 0.001
              identify :pion, against: [:kaon]
            }
          .secondary_vertex_fit([:pip, :pim]) {      # first K_S0, constrained to common vertex
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
          .secondary_vertex_fit([:pip, :pim]) {      # second K_S0
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
          .kinematic_fit([:K_S0, :K_S0, :pip]) {     # 4C fit to the D+ candidate
              nominal
              constrain_four_momentum
              chi2_cut 200
            }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_3773, incMC_3773, exMC_modeII])

# ----------------- Mode III : D0 → K_S0 K_S0 -----------------
alg_modeIII = Algorithm.new("D0ToKsKs")
alg_modeIII.set_header(["D0ToKsKsAlg/D0ToKsKs.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .note(:ks0_selection,
                 "each K_S0 candidate required to satisfy |M(pi+pi-) - M(K_S0)| < 12 MeV/c^2 and "
                 "flight-length significance L/sigma_L > 2; K_S0 daughter tracks use a looser "
                 "|Vz| < 20 cm than the prompt-track |Vz| < 10 cm cut")

sel_modeIII = Selection.new
sel_modeIII.select_track {                     # charged track selection
               cos_theta 0.93
               Vz 10.0
               Vr 1.0
               nChrp "==2"                     # 2 K_S0 daughters (pi+)
               nChrn "==2"                     # 2 K_S0 daughters (pi-)
               nNet  "==0"                     # D0 net charge
             }
           .assign({:chrgp => :pip, :chrgn => :pim})   # K_S0 daughters used without PID
           .secondary_vertex_fit([:pip, :pim]) {       # first K_S0
               build_virtual_particle(:K_S0).by_minimizing_mass_difference
               remove_used_particle_from_candidate_list
             }
           .secondary_vertex_fit([:pip, :pim]) {       # second K_S0
               build_virtual_particle(:K_S0).by_minimizing_mass_difference
               remove_used_particle_from_candidate_list
             }
           .kinematic_fit([:K_S0, :K_S0]) {            # 4C fit to the D0 candidate
               nominal
               constrain_four_momentum
               chi2_cut 200
             }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)
alg_modeIII.execute_on([data_3773, incMC_3773, exMC_modeIII])

# ----------------- Mode IV : D0 → K_S0 K_S0 K_S0 -----------------
alg_modeIV = Algorithm.new("D0ToKsKsKs")
alg_modeIV.set_header(["D0ToKsKsKsAlg/D0ToKsKsKs.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:ks0_selection,
                "each K_S0 candidate required to satisfy |M(pi+pi-) - M(K_S0)| < 12 MeV/c^2 and "
                "flight-length significance L/sigma_L > 2; K_S0 daughter tracks use a looser "
                "|Vz| < 20 cm than the prompt-track |Vz| < 10 cm cut")

sel_modeIV = Selection.new
sel_modeIV.select_track {                      # charged track selection
              cos_theta 0.93
              Vz 10.0
              Vr 1.0
              nChrp "==3"                      # 3 K_S0 daughters (pi+)
              nChrn "==3"                      # 3 K_S0 daughters (pi-)
              nNet  "==0"                      # D0 net charge
            }
          .assign({:chrgp => :pip, :chrgn => :pim})   # K_S0 daughters used without PID
          .secondary_vertex_fit([:pip, :pim]) {       # first K_S0
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
          .secondary_vertex_fit([:pip, :pim]) {       # second K_S0
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
          .secondary_vertex_fit([:pip, :pim]) {       # third K_S0
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
          .kinematic_fit([:K_S0, :K_S0, :K_S0]) {     # 4C fit to the D0 candidate
              nominal
              constrain_four_momentum
              chi2_cut 200
            }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(sel_modeIV)
alg_modeIV.execute_on([data_3773, incMC_3773, exMC_modeIV])