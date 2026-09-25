# =============================================================================
#  BESIII : search for the weak decay J/psi -> D_s^{(*)-} e+ nu_e  (+ c.c.)
#           at sqrt(s) = 3.097 GeV.  BOSS-side (dataset prep + event selection).
# =============================================================================

### ----------------------------------------------------------------- Datasets
jpsi_data  = DatasetManager.real_data.find("708_3097")      # 2.25e8 J/psi events @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # matching inclusive MC
# psi(4040) @ 4.009 GeV (478 pb-1): D_s-efficiency control sample, NOT part of
# the signal selection (only used to measure D_s reconstruction efficiency).
psi4040_data = DatasetManager.real_data.find("703_4009")

### ---------------------------------------------------- Decay cards (8 modes)
# The description asks for the dedicated weak-interaction c->s generator; its
# EvtGen identifier is not exposed by this DSL, so PHSP is used as a stand-in.
# Both charge-conjugate branches are written in each top decay.

# ---- direct D_s channels -----------------------------------------------------
decay_card_ds_kkpi = <<~DECAYCARD
    Decay J/psi
    0.5 D_s-  e+  nu_e        PHSP;
    0.5 D_s+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s-
    1.0 K+  K-  pi-           PHSP;
    Enddecay
    Decay D_s+
    1.0 K-  K+  pi+           PHSP;
    Enddecay
    End
DECAYCARD

decay_card_ds_kkpipi0 = <<~DECAYCARD
    Decay J/psi
    0.5 D_s-  e+  nu_e        PHSP;
    0.5 D_s+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s-
    1.0 K+  K-  pi-  pi0      PHSP;
    Enddecay
    Decay D_s+
    1.0 K-  K+  pi+  pi0      PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma           PHSP;
    Enddecay
    End
DECAYCARD

decay_card_ds_ksk = <<~DECAYCARD
    Decay J/psi
    0.5 D_s-  e+  nu_e        PHSP;
    0.5 D_s+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s-
    1.0 K_S0  K-              PHSP;
    Enddecay
    Decay D_s+
    1.0 K_S0  K+              PHSP;
    Enddecay
    Decay K_S0
    1.0 pi+  pi-              VSS;
    Enddecay
    End
DECAYCARD

decay_card_ds_kskpipi = <<~DECAYCARD
    Decay J/psi
    0.5 D_s-  e+  nu_e        PHSP;
    0.5 D_s+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s-
    1.0 K_S0  K+  pi-  pi-    PHSP;
    Enddecay
    Decay D_s+
    1.0 K_S0  K-  pi+  pi+    PHSP;
    Enddecay
    Decay K_S0
    1.0 pi+  pi-              VSS;
    Enddecay
    End
DECAYCARD

# ---- D_s* -> D_s gamma channels ---------------------------------------------
decay_card_dsstar_kkpi = <<~DECAYCARD
    Decay J/psi
    0.5 D_s*-  e+  nu_e        PHSP;
    0.5 D_s*+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s*-
    1.0 D_s-  gamma            PHSP;
    Enddecay
    Decay D_s*+
    1.0 D_s+  gamma            PHSP;
    Enddecay
    Decay D_s-
    1.0 K+  K-  pi-            PHSP;
    Enddecay
    Decay D_s+
    1.0 K-  K+  pi+            PHSP;
    Enddecay
    End
DECAYCARD

decay_card_dsstar_kkpipi0 = <<~DECAYCARD
    Decay J/psi
    0.5 D_s*-  e+  nu_e        PHSP;
    0.5 D_s*+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s*-
    1.0 D_s-  gamma            PHSP;
    Enddecay
    Decay D_s*+
    1.0 D_s+  gamma            PHSP;
    Enddecay
    Decay D_s-
    1.0 K+  K-  pi-  pi0       PHSP;
    Enddecay
    Decay D_s+
    1.0 K-  K+  pi+  pi0       PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma            PHSP;
    Enddecay
    End
DECAYCARD

decay_card_dsstar_ksk = <<~DECAYCARD
    Decay J/psi
    0.5 D_s*-  e+  nu_e        PHSP;
    0.5 D_s*+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s*-
    1.0 D_s-  gamma            PHSP;
    Enddecay
    Decay D_s*+
    1.0 D_s+  gamma            PHSP;
    Enddecay
    Decay D_s-
    1.0 K_S0  K-               PHSP;
    Enddecay
    Decay D_s+
    1.0 K_S0  K+               PHSP;
    Enddecay
    Decay K_S0
    1.0 pi+  pi-               VSS;
    Enddecay
    End
DECAYCARD

decay_card_dsstar_kskpipi = <<~DECAYCARD
    Decay J/psi
    0.5 D_s*-  e+  nu_e        PHSP;
    0.5 D_s*+  e-  anti-nu_e   PHSP;
    Enddecay
    Decay D_s*-
    1.0 D_s-  gamma            PHSP;
    Enddecay
    Decay D_s*+
    1.0 D_s+  gamma            PHSP;
    Enddecay
    Decay D_s-
    1.0 K_S0  K+  pi-  pi-     PHSP;
    Enddecay
    Decay D_s+
    1.0 K_S0  K-  pi+  pi+     PHSP;
    Enddecay
    Decay K_S0
    1.0 pi+  pi-               VSS;
    Enddecay
    End
DECAYCARD

### ------------------------------------------- Exclusive signal MC (100k each)
exMC_ds_kkpi        = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_ds_kkpi";        c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_ds_kkpi
  c.cross_section = :default
end
exMC_ds_kkpipi0     = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_ds_kkpipi0";     c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_ds_kkpipi0
  c.cross_section = :default
end
exMC_ds_ksk         = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_ds_ksk";         c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_ds_ksk
  c.cross_section = :default
end
exMC_ds_kskpipi     = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_ds_kskpipi";     c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_ds_kskpipi
  c.cross_section = :default
end
exMC_dsstar_kkpi    = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_dsstar_kkpi";    c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_dsstar_kkpi
  c.cross_section = :default
end
exMC_dsstar_kkpipi0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_dsstar_kkpipi0"; c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_dsstar_kkpipi0
  c.cross_section = :default
end
exMC_dsstar_ksk     = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_dsstar_ksk";     c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_dsstar_ksk
  c.cross_section = :default
end
exMC_dsstar_kskpipi = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_jpsi_dsstar_kskpipi"; c.related_dataset = jpsi_data
  c.events = 100_000;                          c.decay_card = decay_card_dsstar_kskpipi
  c.cross_section = :default
end

### ============================================================ Event selection
# Common charged-track quality: |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm.
# Common photon quality: E>25 MeV (barrel) / >50 MeV (endcap), >20 deg from any
# track, EMC timing 0-700 ns (tdc 0-14).

# ------------------------------------------------- Mode A1: D_s- -> K+K-pi-
alg_ds_kkpi = Algorithm.new("DsKKpi")
alg_ds_kkpi.set_header(["DsKKpiAlg/DsKKpi.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .note(:missing_momentum, "neutrino undetected -> no 4C fit; missing four-momentum
                 reconstructed, |p_miss|>50 MeV applied, U_miss=E_miss-|p_miss| stored as the
                 signal variable (unbinned fit in -0.2<U_miss<0.2 GeV done offline)")
           .note(:extra_neutral_energy, "total extra EMC neutral energy required below
                 0.15-0.30 GeV depending on mode and channel")
           .note(:ds_mass_window, "D_s mass window applied at +-3 sigma")
           .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05 applied on top of
                 the probability PID (P(e)>0.001, P(e)>P(K), P(e)>P(pi))")

sel_ds_kkpi = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"      # 4 tracks total (2+,2-)
     nChrn     "==2"
     nNet      "==0"
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6   # exactly one positron
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"          # one positron
     nkp  "==1"          # final multiplicity 1K+ 1K- 1pi-
     nkm  "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])            # remove the identified positron from the charged list
  .partial_miss([3]) {                # recID 3 = nu_e; everything else reconstructed from recoil
     best_combination_by_mass :D_s, 1.9685
   }
alg_ds_kkpi.with_decay_card(decay_card_ds_kkpi).apply(sel_ds_kkpi)

# ------------------------------- Mode A2: D_s- -> K+K-pi-pi0  (pi0 -> gamma gamma)
alg_ds_kkpipi0 = Algorithm.new("DsKKpiPi0")
alg_ds_kkpipi0.set_header(["DsKKpiPi0Alg/DsKKpiPi0.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                    U_miss stored as signal variable (offline fit in -0.2<U_miss<0.2 GeV)")
              .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
              .note(:ds_mass_window, "D_s mass window at +-3 sigma")
              .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_ds_kkpipi0 = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
   }
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end   14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     angle_to_track    20.0
     nGam ">=2"                 # >=2 photons for the pi0
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkp  "==1"                 # 1K+ 1K- 1pi-
     nkm  "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])
  .kalman_kinematic_fit([:gamma, :gamma]) {           # pi0 mass-constrained gammagamma fit
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 100                                     # chi2 < 100
     npi0 ">=1"
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_ds_kkpipi0.with_decay_card(decay_card_ds_kkpipi0).apply(sel_ds_kkpipi0)

# ---------------------------------------- Mode A3: D_s- -> K_S0 K-  (K_S0 -> pi+pi-)
alg_ds_ksk = Algorithm.new("DsKsK")
alg_ds_ksk.set_header(["DsKsKAlg/DsKsK.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .note(:background_veto, "K_S0 from secondary vertex fit requires decay length >2 sigma
                and 0.487<M(pi+pi-)<0.511 GeV/c^2")
          .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                U_miss stored as signal variable")
          .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
          .note(:ds_mass_window, "D_s mass window at +-3 sigma")
          .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_ds_ksk = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkm  "==1"                  # final multiplicity 1K- 1pi+ 1pi-
     npip "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])
  .secondary_vertex_fit([:pip, :pim]) {        # reconstruct K_S0
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_ds_ksk.with_decay_card(decay_card_ds_ksk).apply(sel_ds_ksk)

# ------------------------------- Mode A4: D_s- -> K_S0 K+ pi- pi-  (6 tracks)
alg_ds_kskpipi = Algorithm.new("DsKsKPiPi")
alg_ds_kskpipi.set_header(["DsKsKPiPiAlg/DsKsKPiPi.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .note(:background_veto, "K_S0 secondary vertex fit: decay length >2 sigma and
                    0.487<M(pi+pi-)<0.511 GeV/c^2")
              .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                    U_miss stored as signal variable")
              .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
              .note(:ds_mass_window, "D_s mass window at +-3 sigma")
              .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_ds_kskpipi = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==3"              # 6 tracks total (3+,3-)
     nChrn     "==3"
     nNet      "==0"
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkp  "==1"                   # final multiplicity 1K+ 1pi+ 3pi-
     npip "==1"
     npim "==3"
   }
  .remove([:lp <= :chrgp])
  .secondary_vertex_fit([:pip, :pim]) {
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_ds_kskpipi.with_decay_card(decay_card_ds_kskpipi).apply(sel_ds_kskpipi)

# =============================================================================
#  D_s* -> D_s gamma channels (add >=1 photon, and Delta-M = m(D_s gamma)-m(D_s))
# =============================================================================

# ---------------------------------------- Mode B1: D_s*-, D_s- -> K+K-pi-
alg_dsstar_kkpi = Algorithm.new("DsStarKKpi")
alg_dsstar_kkpi.set_header(["DsStarKKpiAlg/DsStarKKpi.h"])
               .set_constant({"ECMS" => [:double, 3.097]})
               .note(:dsstar_delta_m, "D_s*: D_s combined with a photon requiring
                     0.125 < DeltaM = m(D_s gamma)-m(D_s) < 0.150 GeV/c^2")
               .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                     U_miss stored as signal variable")
               .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
               .note(:ds_mass_window, "D_s mass window at +-3 sigma")
               .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_dsstar_kkpi = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
   }
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end   14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     angle_to_track    20.0
     nGam ">=1"                  # the D_s* photon
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkp  "==1"
     nkm  "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_dsstar_kkpi.with_decay_card(decay_card_dsstar_kkpi).apply(sel_dsstar_kkpi)

# ------------------------------ Mode B2: D_s*-, D_s- -> K+K-pi-pi0 (>=3 photons)
alg_dsstar_kkpipi0 = Algorithm.new("DsStarKKpiPi0")
alg_dsstar_kkpipi0.set_header(["DsStarKKpiPi0Alg/DsStarKKpiPi0.h"])
                  .set_constant({"ECMS" => [:double, 3.097]})
                  .note(:dsstar_delta_m, "D_s*: DeltaM = m(D_s gamma)-m(D_s) in
                        [0.125, 0.150] GeV/c^2")
                  .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                        U_miss stored as signal variable")
                  .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
                  .note(:ds_mass_window, "D_s mass window at +-3 sigma")
                  .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_dsstar_kkpipi0 = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
   }
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end   14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     angle_to_track    20.0
     nGam ">=3"                  # D_s* photon + 2 photons from pi0
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkp  "==1"
     nkm  "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 100
     npi0 ">=1"
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_dsstar_kkpipi0.with_decay_card(decay_card_dsstar_kkpipi0).apply(sel_dsstar_kkpipi0)

# --------------------------------- Mode B3: D_s*-, D_s- -> K_S0 K-
alg_dsstar_ksk = Algorithm.new("DsStarKsK")
alg_dsstar_ksk.set_header(["DsStarKsKAlg/DsStarKsK.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .note(:dsstar_delta_m, "D_s*: DeltaM = m(D_s gamma)-m(D_s) in
                    [0.125, 0.150] GeV/c^2")
              .note(:background_veto, "K_S0 secondary vertex fit: decay length >2 sigma and
                    0.487<M(pi+pi-)<0.511 GeV/c^2")
              .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                    U_miss stored as signal variable")
              .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
              .note(:ds_mass_window, "D_s mass window at +-3 sigma")
              .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_dsstar_ksk = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
   }
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end   14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     angle_to_track    20.0
     nGam ">=1"
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkm  "==1"
     npip "==1"
     npim "==1"
   }
  .remove([:lp <= :chrgp])
  .secondary_vertex_fit([:pip, :pim]) {
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_dsstar_ksk.with_decay_card(decay_card_dsstar_ksk).apply(sel_dsstar_ksk)

# -------------------- Mode B4: D_s*-, D_s- -> K_S0 K+ pi- pi-  (6 tracks, >=1 photon)
alg_dsstar_kskpipi = Algorithm.new("DsStarKsKPiPi")
alg_dsstar_kskpipi.set_header(["DsStarKsKPiPiAlg/DsStarKsKPiPi.h"])
                  .set_constant({"ECMS" => [:double, 3.097]})
                  .note(:dsstar_delta_m, "D_s*: DeltaM = m(D_s gamma)-m(D_s) in
                        [0.125, 0.150] GeV/c^2")
                  .note(:background_veto, "K_S0 secondary vertex fit: decay length >2 sigma and
                        0.487<M(pi+pi-)<0.511 GeV/c^2")
                  .note(:missing_momentum, "neutrino undetected -> no 4C fit; |p_miss|>50 MeV,
                        U_miss stored as signal variable")
                  .note(:extra_neutral_energy, "total extra EMC neutral energy below 0.15-0.30 GeV")
                  .note(:ds_mass_window, "D_s mass window at +-3 sigma")
                  .note(:pid_correction_method, "positron E/p window 0.80<E/p<1.05")
sel_dsstar_kskpipi = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==3"
     nChrn     "==3"
     nNet      "==0"
   }
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end   14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     angle_to_track    20.0
     nGam ">=1"
   }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
     identify :kaon, against: [:pion, :proton]
     identify :pion, against: [:kaon, :proton]
     nlp  "==1"
     nkp  "==1"
     npip "==1"
     npim "==3"
   }
  .remove([:lp <= :chrgp])
  .secondary_vertex_fit([:pip, :pim]) {
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
   }
  .partial_miss([3]) {
     best_combination_by_mass :D_s, 1.9685
   }
alg_dsstar_kskpipi.with_decay_card(decay_card_dsstar_kskpipi).apply(sel_dsstar_kskpipi)

### ------------------------------------------------------------ Execute on data
# The four D_s modes are also run on the psi(4040) control sample to measure the
# D_s reconstruction efficiency (it is not part of the signal sample).

alg_ds_kkpi.with_decay_card(decay_card_ds_kkpi).apply(sel_ds_kkpi)
root_files = []
root_files += alg_ds_kkpi.execute_on([jpsi_data, jpsi_incMC, psi4040_data, exMC_ds_kkpi])
root_files += alg_ds_kkpipi0.execute_on([jpsi_data, jpsi_incMC, psi4040_data, exMC_ds_kkpipi0])
root_files += alg_ds_ksk.execute_on([jpsi_data, jpsi_incMC, psi4040_data, exMC_ds_ksk])
root_files += alg_ds_kskpipi.execute_on([jpsi_data, jpsi_incMC, psi4040_data, exMC_ds_kskpipi])
root_files += alg_dsstar_kkpi.execute_on([jpsi_data, jpsi_incMC, exMC_dsstar_kkpi])
root_files += alg_dsstar_kkpipi0.execute_on([jpsi_data, jpsi_incMC, exMC_dsstar_kkpipi0])
root_files += alg_dsstar_ksk.execute_on([jpsi_data, jpsi_incMC, exMC_dsstar_ksk])
root_files += alg_dsstar_kskpipi.execute_on([jpsi_data, jpsi_incMC, exMC_dsstar_kskpipi])