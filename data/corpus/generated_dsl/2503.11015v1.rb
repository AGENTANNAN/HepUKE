# ============================================================================
#  Search for X(1^-+) in  e+e- -> gamma Ds+ Ds1(2536)-
#    Ds1(2536)- -> anti-D*0 K- ,  anti-D*0 -> anti-D0 pi0
#    Mode 1: Ds+ -> K+ K- pi+   (phi / K*(892) windows + helicity cut)
#    Mode 2: Ds+ -> K_S0 K+     (K_S0 -> pi+ pi-)
#  12 energy points, 4.612 - 4.951 GeV  (5.8 fb^-1)
# ============================================================================

### ------------------------------- Datasets ------------------------------- ###
# 12 scan points: BOSS 706 (4610-4700 MeV) and BOSS 707 (4740-4946 MeV)
scan_names = %w[
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
real_data = scan_names.map { |n| DatasetManager.real_data.find(n) }     # real data
inc_mc    = scan_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # inclusive MC

### ------------------------------- Decay cards ------------------------------- ###
# Mode 1: Ds+ -> K+ K- pi+
decay_card_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma Ds+ Ds1- PHSP;
  Enddecay

  Decay Ds+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay Ds1-
  1.000 anti-D*0 K- PHSP;
  Enddecay

  Decay anti-D*0
  1.000 anti-D0 pi0 PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2: Ds+ -> K_S0 K+ , K_S0 -> pi+ pi-
decay_card_mode2 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma Ds+ Ds1- PHSP;
  Enddecay

  Decay Ds+
  1.000 K_S0 K+ PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay Ds1-
  1.000 anti-D*0 K- PHSP;
  Enddecay

  Decay anti-D*0
  1.000 anti-D0 pi0 PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### -------------------- Exclusive signal MC (one per mode) -------------------- ###
# 500k events for each Ds+ decay mode, generated at every energy point
exMC_mode1 = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_X1pp_DsDs1_mode1"
  config.events        = 500_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_X1pp_DsDs1_mode2"
  config.events        = 500_000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

### ------------------------------- Algorithms ------------------------------- ###
# --- Mode 1: Ds+ -> K+ K- pi+ ---
alg_name_mode1 = "X1ppDsDs1Mode1"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 4.612]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_scan, "Analysis runs over 12 scan points from 4.612 to 4.951 GeV; ECMS "
                           "must be set to the measured centre-of-mass energy of the "
                           "individual point each time the algorithm is applied.")
         .note(:phi_mass_window, "Mode 1: the K+K- pair from Ds+ -> K+K-pi+ is required to "
                                 "lie in the phi(1020) window 1.019 +/- 0.010 GeV.")
         .note(:kstar_mass_window, "Mode 1: the K-pi+ pair from Ds+ -> K+K-pi+ is required to "
                                   "lie in the K*(892) window 0.892 +/- 0.050 GeV.")
         .note(:helicity_angle, "Mode 1: the Ds+ -> K+K-pi+ helicity angle is required to "
                                "satisfy |cos(theta_hel)| > 0.5.")

# --- Mode 2: Ds+ -> K_S0 K+ ---
alg_name_mode2 = "X1ppDsDs1Mode2"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 4.612]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_scan, "Analysis runs over 12 scan points from 4.612 to 4.951 GeV; ECMS "
                           "must be set to the measured centre-of-mass energy of the "
                           "individual point each time the algorithm is applied.")
         .note(:ks_mass_window, "Mode 2: the K_S0 -> pi+pi- candidate is required to have an "
                                "invariant mass in (487, 511) MeV/c^2 after the secondary "
                                "vertex fit.")
         .note(:ks_decay_length, "Mode 2: the K_S0 decay length is required to exceed twice "
                                 "the vertex resolution.")

### ------------------------------- Event selection ------------------------------- ###
# Common preselection shared by the two Ds+ modes
common_selection = Selection.new
    .select_track {                 # charged-track quality cuts
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        1.0               # Vr < 1 cm
        nChrp     ">=2"
        nChrn     ">=2"
        nNet      "==0"
    }
    .select_photon {                # photon selection
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0      # > 10 deg from the nearest charged track
        energyThreshold_b 0.025     # 25 MeV in the barrel
        energyThreshold_e 0.050     # 50 MeV in the endcap
        nGam              ">=2"     # at least two photons
    }

# ---- Mode 1: Ds+ -> K+ K- pi+ ----
# Final state K+ K- pi+ (Ds+) plus K- (Ds1-): require one K+ and one K-; the two
# remaining tracks are identified as pions (one pi+, one pi-).
sel_mode1 = common_selection.dup
    .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]   # K+ and K-
        identify :pion, against: [:kaon, :proton]   # remaining tracks as pions
        nkp  "==1"
        nkm  "==1"
        npip "==1"
        npim "==1"
    }
    # Partial reconstruction of gamma, Ds+(K+K-pi+) and the K- from Ds1-
    # (rec ids: 1=gamma, 2=Ds+, 4=K+_Ds+, 5=K-_Ds+, 6=pi+_Ds+, 8=K-_Ds1-)
    .partial_rec([1, 2, 4, 5, 6, 8]) {
        best_combination_by_mass :Ds_plus, 1.968   # Ds+ mass constraint 1.968 GeV
        require_recoil_mass 2.000, 2.020           # recoil of gamma Ds+ K- = D*0 window
    }

# ---- Mode 2: Ds+ -> K_S0 K+ , K_S0 -> pi+ pi- ----
sel_mode2 = common_selection.dup
    .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp "==1"                                   # exactly one K+ from Ds+
    }
    .secondary_vertex_fit([:pip, :pim]) {           # K_S0 -> pi+ pi-
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .assign({:chrgn => :km})                        # remaining negative track is the K- from Ds1-
    # Partial reconstruction of gamma, Ds+(K_S0 K+) and the K- from Ds1-
    # (rec ids: 1=gamma, 2=Ds+, 4=K_S0, 5=K+_Ds+, 6=pi+_K_S0, 7=pi-_K_S0, 9=K-_Ds1-)
    .partial_rec([1, 2, 4, 5, 6, 7, 9]) {
        best_combination_by_mass :Ds_plus, 1.968   # Ds+ mass constraint 1.968 GeV
        require_recoil_mass 2.000, 2.020           # recoil of gamma Ds+ K- = D*0 window
    }

### ------------------------------- Wire up ------------------------------- ###
alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)

root_files_mode1 = alg_mode1.execute_on(real_data + inc_mc + exMC_mode1)
root_files_mode2 = alg_mode2.execute_on(real_data + inc_mc + exMC_mode2)