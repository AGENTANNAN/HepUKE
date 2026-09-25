# BESIII: branching fractions for D+ -> K_S0 K_S0 K+, D+ -> K_S0 K_S0 pi+,
# D0 -> K_S0 K_S0 and D0 -> K_S0 K_S0 K_S0 (2.93 fb^-1 at sqrt(s) = 3.773 GeV)
#
# Single-tag method at the psi(3770) peak: a D (or anti-D) is reconstructed in
# the signal mode and the recoil side is not reconstructed.  The signal decays
# are:
#   Mode I   D+  -> K_S0 K_S0 K+
#   Mode II  D+  -> K_S0 K_S0 pi+
#   Mode III D0  -> K_S0 K_S0
#   Mode IV  D0  -> K_S0 K_S0 K_S0
# Each mode has its own final-state multiplicity and its own Delta E window, so
# each gets its own Algorithm/Selection chain (Rule T1).  The D yields are
# extracted from fits to M_BC in 2D/3D K_S0 signal and sideband regions; the
# sideband subtraction and the peaking-background normalisation are ROOT-level
# operations.

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")    # 2.93 fb^-1 at sqrt(s) = 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773") # 10x data: DDbar, ISR psi(3686)/Jpsi, qqbar, Bhabha, dimuon, ditau

### Decay cards for the exclusive signal MC ###
# Mode I: D+ -> K_S0 K_S0 K+ (phase space)
decay_card_dp_ksksk = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-                                           PHSP;
    Enddecay

    Decay D+
    1.0000  K_S0  K_S0  K+                                   PHSP;
    Enddecay

    Decay D-
    1.0000  anything                                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                         PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: D+ -> K_S0 K_S0 pi+.  90% proceeds through the intermediate
# K*(892)+ -> K_S0 pi+ with the other K_S0 recoiling, 10% is the direct
# three-body phase-space decay.
decay_card_dp_kskspi = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-                                           PHSP;
    Enddecay

    Decay D+
    0.9000  K_S0  K*+                                        PHSP;
    0.1000  K_S0  K_S0  pi+                                  PHSP;
    Enddecay

    Decay K*+
    1.0000  K_S0  pi+                                        PHSP;
    Enddecay

    Decay D-
    1.0000  anything                                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                         PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: D0 -> K_S0 K_S0 (phase space)
decay_card_d0_ksks = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D0  anti-D0                                      PHSP;
    Enddecay

    Decay D0
    1.0000  K_S0  K_S0                                       PHSP;
    Enddecay

    Decay anti-D0
    1.0000  anything                                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                         PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV: D0 -> K_S0 K_S0 K_S0 (phase space)
decay_card_d0_ksksks = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D0  anti-D0                                      PHSP;
    Enddecay

    Decay D0
    1.0000  K_S0  K_S0  K_S0                                 PHSP;
    Enddecay

    Decay anti-D0
    1.0000  anything                                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                         PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_ksksk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_KsKsK"
  config.related_dataset = psi3770_data
  config.events          = 300000
  config.decay_card      = decay_card_dp_ksksk
  config.cross_section   = :default
end

exMC_dp_kskspi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_KsKsPi"
  config.related_dataset = psi3770_data
  config.events          = 300000
  config.decay_card      = decay_card_dp_kskspi
  config.cross_section   = :default
end

exMC_d0_ksks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0_KsKs"
  config.related_dataset = psi3770_data
  config.events          = 300000
  config.decay_card      = decay_card_d0_ksks
  config.cross_section   = :default
end

exMC_d0_ksksks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0_KsKsKs"
  config.related_dataset = psi3770_data
  config.events          = 300000
  config.decay_card      = decay_card_d0_ksksks
  config.cross_section   = :default
end

datasets = [psi3770_data, psi3770_incMC]

# BOSS-side requirements shared by all four single-tag channels.
track_pid_ks0_note = "Charged tracks: |cos(theta)| < 0.93; good tracks not used for the " \
                     "K_S0 reconstruction must originate within V_xy < 1.0 cm and " \
                     "V_z < 10.0 cm.  Charged kaon/pion separation uses the combined dE/dx " \
                     "and TOF confidence levels: a track is identified as a kaon (pion) if " \
                     "CL_K > CL_pi (CL_pi > CL_K).  K_S0 -> pi+ pi-: both daughter pions must " \
                     "satisfy V_z < 20.0 cm, are used without PID requirements, are " \
                     "constrained to a common vertex, and the candidate is accepted if " \
                     "|M(pi+ pi-) - M(K_S0)^PDG| < 12 MeV/c^2 and L/sigma_L > 2.  Candidate " \
                     "selection uses Delta E = E_D - E_beam and " \
                     "M_BC = sqrt(E_beam^2/c^4 - |p_D|^2/c^2); when several candidates pass " \
                     "the selection requirements, only the one with the minimum |Delta E| is " \
                     "kept."

ks_sideband_note = "The K_S0 mesons define a 2D (M(pi+pi-)_1, M(pi+pi-)_2) plane: the signal " \
                   "region has both pairs in the K_S0 signal window, sideband 1 has one pair " \
                   "in the sideband, sideband 2 has both.  The K_S0 sideband is " \
                   "0.020 < |M(pi+ pi-) - M(K_S0)^PDG| < 0.044 GeV/c^2.  The D yield is " \
                   "extracted from fits to M_BC in each region (MC signal shape convoluted " \
                   "with a Gaussian with free resolution-difference parameters, ARGUS " \
                   "combinatorial background with a 1.8865 GeV/c^2 endpoint) and the net yield " \
                   "is N_net = N_sig - N_sb1/2 + N_sb2/4 - N_other, where N_other is the " \
                   "peaking background (D+ -> K_S0 K_L0 K+, D- -> K_S0 X) taken from the " \
                   "inclusive MC.  The combinatorial pi+pi- background in the M(pi+pi-) " \
                   "distribution is assumed flat, giving a background ratio of 0.5 between " \
                   "the K_S0 signal and sideband regions.  Applied in the ROOT analysis."

ksksks_sideband_note = "The three K_S0 mesons define a 3D (M(pi+pi-)_1, M(pi+pi-)_2, " \
                       "M(pi+pi-)_3) plane: the signal region has all three pairs in the K_S0 " \
                       "signal window, sideband i has i of the three pairs in the K_S0 " \
                       "sideband and the rest in the signal region " \
                       "(0.020 < |M(pi+pi-) - M(K_S0)^PDG| < 0.044 GeV/c^2).  The net yield is " \
                       "N_net = N_sig - N_sb1/2 + N_sb2/4 - N_sb3/8 - N_other, with the " \
                       "peaking background dominated by D0 -> K_S0 K_S0 K_L0 versus " \
                       "anti-D0 -> K_S0 X taken from the inclusive MC; the sideband 3 " \
                       "contribution is negligible since few events survive.  M_BC fits use an " \
                       "MC-derived signal shape convoluted with a Gaussian and an ARGUS " \
                       "background with a 1.8865 GeV/c^2 endpoint.  Applied in the ROOT " \
                       "analysis."

########################################################################
# Mode I: D+ -> K_S0 K_S0 K+   (2 K_S0 -> 2 pi+ 2 pi- , plus one prompt K+)
########################################################################
alg_name_I = "DpToKsKsK"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .note(:common_track_pid_ks0_selection, track_pid_ks0_note)
     .note(:delta_e_requirement,
           "Mode-dependent Delta E window from fits to the data and MC Delta E distributions " \
           "at +-3 sigma: data (-17, +19) MeV, MC (-16, +16) MeV.")
     .note(:ks_sideband_method, ks_sideband_note)

sel_I = Selection.new
sel_I.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr         1.0
        nChrp     "==3"     # K+ plus the two pi+ from the K_S0 pair decays
        nChrn     "==2"     # the two pi- from the K_S0 pair decays
        nNet      "==1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion]     # the prompt K+
      }
      .remove([:kp <= :chrgp])               # the remaining tracks are the K_S0 daughters
      .assign({ :chrgp => :pip, :chrgn => :pim })
      .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:K_S0, :K_S0, :kp]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_I.with_decay_card(decay_card_dp_ksksk).apply(sel_I)
alg_I.execute_on(datasets + [exMC_dp_ksksk])

########################################################################
# Mode II: D+ -> K_S0 K_S0 pi+   (2 K_S0 plus one prompt pi+)
########################################################################
alg_name_II = "DpToKsKsPi"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .note(:common_track_pid_ks0_selection, track_pid_ks0_note)
      .note(:delta_e_requirement,
            "Mode-dependent Delta E window from fits to the data and MC Delta E " \
            "distributions at +-3 sigma: data (-17, +17) MeV, MC (-17, +16) MeV.")
      .note(:ks_sideband_method, ks_sideband_note)
      .note(:mc_modeling,
            "Signal MC is a mixed sample: 90% D+ -> K_S0 K*(892)+ with " \
            "K*(892)+ -> K_S0 pi+, 10% direct three-body phase-space decay; the efficiency " \
            "difference from the a0(980) / f0(980) sub-resonances is taken as a 1.0% model " \
            "uncertainty.")

sel_II = Selection.new
sel_II.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr         1.0
         nChrp     "==3"    # prompt pi+ plus the two pi+ from the K_S0 pair decays
         nChrn     "==2"    # the two pi- from the K_S0 pair decays
         nNet      "==1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon]    # the prompt pi+; the K_S0 pions need no PID
       }
       .assign({ :chrgp => :pip, :chrgn => :pim })
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .kinematic_fit([:K_S0, :K_S0, :pip]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }

alg_II.with_decay_card(decay_card_dp_kskspi).apply(sel_II)
alg_II.execute_on(datasets + [exMC_dp_kskspi])

########################################################################
# Mode III: D0 -> K_S0 K_S0   (all tracks are K_S0 daughters)
########################################################################
alg_name_III = "D0ToKsKs"
alg_III = Algorithm.new(alg_name_III)
alg_III.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .note(:common_track_pid_ks0_selection, track_pid_ks0_note)
       .note(:delta_e_requirement,
             "Mode-dependent Delta E window from fits to the data and MC Delta E " \
             "distributions at +-3 sigma: data (-19, +17) MeV, MC (-17, +14) MeV.")
       .note(:ks_sideband_method, ks_sideband_note)
       .note(:two_body_topology,
             "This two-body mode has no residual peaking background (N_other = 0); the net " \
             "yield follows the same 2D sideband subtraction as the other channels.")

sel_III = Selection.new
sel_III.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr         1.0
          nChrp     "==2"
          nChrn     "==2"
          nNet      "==0"
        }
        .assign({ :chrgp => :pip, :chrgn => :pim })
        .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        .kinematic_fit([:K_S0, :K_S0]) {
          nominal
          constrain_four_momentum
          chi2_cut 200
        }

alg_III.with_decay_card(decay_card_d0_ksks).apply(sel_III)
alg_III.execute_on(datasets + [exMC_d0_ksks])

########################################################################
# Mode IV: D0 -> K_S0 K_S0 K_S0   (all tracks are K_S0 daughters)
########################################################################
alg_name_IV = "D0ToKsKsKs"
alg_IV = Algorithm.new(alg_name_IV)
alg_IV.set_header(["#{alg_name_IV}Alg/#{alg_name_IV}.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .note(:common_track_pid_ks0_selection, track_pid_ks0_note)
      .note(:delta_e_requirement,
            "Mode-dependent Delta E window from fits to the data and MC Delta E " \
            "distributions at +-3 sigma: data (-14, +16) MeV, MC (-13, +13) MeV.")
      .note(:ks_sideband_method, ksksks_sideband_note)

sel_IV = Selection.new
sel_IV.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr         1.0
         nChrp     "==3"
         nChrn     "==3"
         nNet      "==0"
       }
       .assign({ :chrgp => :pip, :chrgn => :pim })
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .kinematic_fit([:K_S0, :K_S0, :K_S0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }

alg_IV.with_decay_card(decay_card_d0_ksksks).apply(sel_IV)
alg_IV.execute_on(datasets + [exMC_d0_ksksks])
