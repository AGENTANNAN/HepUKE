# =============================================================================
# BOSS event-selection specification
#   e+e- -> Ds+ Ds1(2536)- ,  Ds1(2536)- -> D*0bar K-
#                            -> (D0bar pi0) K- , D0bar -> K+ pi- , pi0 -> gamma gamma
#   e+e- -> Ds+ Ds2*(2573)- , Ds2*(2573)- -> D0bar K- , D0bar -> K+ pi-
#   both with Ds+ -> K- K+ pi+
# 15 BESIII scan points, 4.530 - 4.946 GeV. Real data only (no inclusive MC here).
# =============================================================================

### ------------------------------ Datasets ------------------------------ ###
data_4530 = DatasetManager.real_data.find("703_4530")   # 4.530 GeV
data_4575 = DatasetManager.real_data.find("703_4575")   # 4.575 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.610 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.620 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.640 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.660 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.680 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.700 GeV
data_4740 = DatasetManager.real_data.find("707_4740")   # 4.740 GeV
data_4750 = DatasetManager.real_data.find("707_4750")   # 4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")   # 4.780 GeV
data_4840 = DatasetManager.real_data.find("707_4840")   # 4.840 GeV
data_4914 = DatasetManager.real_data.find("707_4914")   # 4.914 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # 4.946 GeV

# All 15 scan points analysed together (no inclusive MC sample is used)
scan_data = [data_4530, data_4575, data_4600, data_4610, data_4620,
             data_4640, data_4660, data_4680, data_4700, data_4740,
             data_4750, data_4780, data_4840, data_4914, data_4946]

### ----------------------------- Decay cards ----------------------------- ###
# Mode I: e+e- -> Ds+ Ds1(2536)-
decay_card_ds1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D_s+  D_s1-                         PHSP;
    Enddecay

    Decay D_s+
    1.0000  K+  K-  pi+                         PHSP;
    Enddecay

    Decay D_s1-
    1.0000  anti-D*0  K-                        PHSP;
    Enddecay

    Decay anti-D*0
    1.0000  anti-D0  pi0                        PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                        PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: e+e- -> Ds+ Ds2*(2573)-
decay_card_ds2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D_s+  D_s2*-                        PHSP;
    Enddecay

    Decay D_s+
    1.0000  K+  K-  pi+                         PHSP;
    Enddecay

    Decay D_s2*-
    1.0000  anti-D0  K-                         PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                            PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------- Exclusive MC samples ------------------------ ###
# 200k-event signal MC per mode, generated at every one of the 15 energy points
exMCs_ds1 = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "exmc_dsp_ds1_2536"   # auto-suffixed per energy point
  config.events      = 200_000
  config.decay_card  = decay_card_ds1
  config.cross_section = :default
end
exMCs_ds1.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

exMCs_ds2 = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "exmc_dsp_ds2_2573"   # auto-suffixed per energy point
  config.events      = 200_000
  config.decay_card  = decay_card_ds2
  config.cross_section = :default
end
exMCs_ds2.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### ------------------------------ Algorithms ----------------------------- ###
alg_name_ds1 = "Ds1_2536"
alg_ds1 = Algorithm.new(alg_name_ds1)
alg_ds1.set_header(["#{alg_name_ds1}Alg/#{alg_name_ds1}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})   # representative scan energy; set per point when running

alg_name_ds2 = "Ds2_2573"
alg_ds2 = Algorithm.new(alg_name_ds2)
alg_ds2.set_header(["#{alg_name_ds2}Alg/#{alg_name_ds2}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})   # representative scan energy; set per point when running

### ----------------- Shared event selection (both modes) ---------------- ###
event_selection_common = Selection.new
  .select_track {                       # charged-track quality cuts
      cos_theta   0.93                  # |cos(theta)| < 0.93
      Vz          10.0                  # |Vz| < 10 cm
      Vr          1.0                   # Vr < 1 cm
      nChrp       ">=2"                 # at least two positive tracks
      nChrn       ">=2"                 # at least two negative tracks
  }
  # no photon selection is applied in this analysis
  .pid(method: :probability) {          # K / pi separation, no lepton identification
      prob_cut    0.001                 # PID probability > 0.001
      identify :kaon, against: [:pion]  # K+ and K- versus pions
      identify :pion, against: [:kaon]  # pi+ versus kaons
      nkm  ">=2"                        # at least two K-
      nkp  ">=1"                        # at least one K+
      npip ">=1"                        # at least one pi+
  }
  # Ds+ -> K- K+ pi+ : common vertex for the three tracks
  .secondary_vertex_fit([:km, :kp, :pip]) {
      build_virtual_particle(:Ds_p).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Ds+ -> K_S0 K+ sub-mode : K_S0 -> pi+ pi- from a secondary vertex
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }

### --------------------- Mode I : Ds1(2536)- selection ------------------- ###
ds1_selection = event_selection_common.dup
  .kinematic_fit([:Ds_p, :km, :D_star0_bar]) {
      nominal                                              # nominal fit (corrected 4-momenta saved)
      miss_track_of :D_star0_bar                           # D*0bar is not reconstructed (missing mass constraint)
      invariant_mass_of(:Ds_p).constrain_to_nominal_mass_of(:Ds_p)   # Ds+ mass constraint
      # additional intermediate windows of the Ds+ sub-modes
      invariant_mass_of(:kp, :km).within(1.004, 1.034)     # Ds+ -> phi pi+   : M(K+ K-)
      invariant_mass_of(:km, :pip).within(0.832, 0.928)    # Ds+ -> K*0bar K+ : M(K- pi+)
      constrain_four_momentum                              # 4C energy-momentum constraint
      chi2_cut 200                                         # loose BOSS-level cut (tightened in ROOT)
  }

### -------------------- Mode II : Ds2*(2573)- selection ------------------ ###
ds2_selection = event_selection_common.dup
  .kinematic_fit([:Ds_p, :km, :D0_bar]) {
      nominal
      miss_track_of :D0_bar                                # D0bar is not reconstructed (missing mass constraint)
      invariant_mass_of(:Ds_p).constrain_to_nominal_mass_of(:Ds_p)   # Ds+ mass constraint
      invariant_mass_of(:kp, :km).within(1.004, 1.034)     # Ds+ -> phi pi+   : M(K+ K-)
      invariant_mass_of(:km, :pip).within(0.832, 0.928)    # Ds+ -> K*0bar K+ : M(K- pi+)
      constrain_four_momentum                              # 4C energy-momentum constraint
      chi2_cut 200                                         # loose BOSS-level cut (tightened in ROOT)
  }

### ----------- Inexpressible BOSS-side procedures (notes) --------------- ###
alg_ds1
  .note(:ds_plus_vertex_chi2,
        "Ds+ -> K- K+ pi+ candidates are built from a common vertex; only combinations with
         vertex chi2 < 100 and |M(K- K+ pi+) - m_Ds+| < 8 MeV/c^2 are retained.
         The combinatorial loop iterates over the two K- assignments (one K- to Ds+, one to Ds1-).")
  .note(:recoil_mass_window,
        "|RM(Ds+ K-) - m(D*0bar)| < 9 MeV/c^2 for the Ds+ -> K- K+ pi+ sub-mode;
         tightened to < 7 MeV/c^2 for the Ds+ -> K_S0 K+ sub-mode.")
  .note(:intermediate_mass_windows,
        "Ds+ -> phi pi+ : M(K+ K-) in [1.004, 1.034] GeV/c^2;
         Ds+ -> K*0bar K+ : M(K- pi+) in [0.832, 0.928] GeV/c^2;
         Ds+ -> K_S0 K+ : |M(K_S0 K+) - m_Ds+| < 8 MeV/c^2 with K_S0 from a secondary vertex.")
  .note(:helicity_angle,
        "|cos(theta)| > 0.4 for the K+ from phi (Ds+ -> phi pi+) and > 0.52 for the K- from K*0bar
         (Ds+ -> K*0bar K+); applied at ROOT level on the helicity angles.")
  .note(:branching_fraction_combination,
        "Ds+ -> K- K+ pi+, phi pi+, K*0bar K+ and K_S0 K+ sub-modes are combined according to
         their branching fractions when forming the Ds+ yield.")
  .note(:ecms_scan,
        "The ECMS constant must be set to the CMS energy of each of the 15 scan points
         (4.530 - 4.946 GeV) before running; the value declared above is representative.")

alg_ds2
  .note(:ds_plus_vertex_chi2,
        "Ds+ -> K- K+ pi+ candidates are built from a common vertex; only combinations with
         vertex chi2 < 100 and |M(K- K+ pi+) - m_Ds+| < 8 MeV/c^2 are retained.
         The combinatorial loop iterates over the two K- assignments (one K- to Ds+, one to Ds2*-).")
  .note(:recoil_mass_window,
        "|RM(Ds+ K-) - m(D0bar)| < 11 MeV/c^2 for the Ds+ -> K- K+ pi+ sub-mode;
         tightened to < 9 MeV/c^2 for the Ds+ -> K_S0 K+ sub-mode.")
  .note(:intermediate_mass_windows,
        "Ds+ -> phi pi+ : M(K+ K-) in [1.004, 1.034] GeV/c^2;
         Ds+ -> K*0bar K+ : M(K- pi+) in [0.832, 0.928] GeV/c^2;
         Ds+ -> K_S0 K+ : |M(K_S0 K+) - m_Ds+| < 8 MeV/c^2 with K_S0 from a secondary vertex.")
  .note(:helicity_angle,
        "|cos(theta)| > 0.4 for the K+ from phi (Ds+ -> phi pi+) and > 0.52 for the K- from K*0bar
         (Ds+ -> K*0bar K+); applied at ROOT level on the helicity angles.")
  .note(:branching_fraction_combination,
        "Ds+ -> K- K+ pi+, phi pi+, K*0bar K+ and K_S0 K+ sub-modes are combined according to
         their branching fractions when forming the Ds+ yield.")
  .note(:ecms_scan,
        "The ECMS constant must be set to the CMS energy of each of the 15 scan points
         (4.530 - 4.946 GeV) before running; the value declared above is representative.")

### ------------------------------ Execution ------------------------------ ###
alg_ds1.with_decay_card(decay_card_ds1).apply(ds1_selection)
alg_ds2.with_decay_card(decay_card_ds2).apply(ds2_selection)

root_files_ds1 = alg_ds1.execute_on(scan_data + exMCs_ds1)
root_files_ds2 = alg_ds2.execute_on(scan_data + exMCs_ds2)