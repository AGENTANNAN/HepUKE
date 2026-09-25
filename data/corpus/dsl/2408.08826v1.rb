### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # (10087 +/- 44) x 10^6 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ------------------------------------------------------------------------------
# Decay cards
# ------------------------------------------------------------------------------
# Mode I: J/psi -> gamma D0, D0 -> K- pi+
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma   D0                     JPE;
    Enddecay

    Decay D0
    1.0000  K-   pi+                        PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma D0, D0 -> K- pi+ pi0
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma   D0                     JPE;
    Enddecay

    Decay D0
    1.0000  K-   pi+   pi0                  PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma   gamma                   PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: J/psi -> gamma D0, D0 -> K- pi+ pi+ pi-
decay_card_modeIII = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma   D0                     JPE;
    Enddecay

    Decay D0
    1.0000  K-   pi+   pi+   pi-            PHSP;
    Enddecay

    End
DECAYCARD

exMC_jpsi_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaD0_Kpi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_jpsi_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaD0_Kpipi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_jpsi_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaD0_Kpipipi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

# ==============================================================================
# Common event selection template
# ==============================================================================
def build_common_selection
  Selection.new
    .select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
     }
    .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start     0
       tdc_emc_end       14
     }
end

# ==============================================================================
# Mode I: J/psi -> gamma D0, D0 -> K- pi+
# ==============================================================================
alg_modeI = Algorithm.new("JpsiGammaD0_Kpi")
alg_modeI.set_header(["JpsiGammaD0Alg/JpsiGammaD0.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = build_common_selection
  .pid(method: :probability) {
     prob_cut 0.001
     identify :kaon, against: [:pion]
     identify :pion, against: [:kaon]
     nkp "==1"
     nkm "==0"
     npip "==1"
     npim "==0"
   }
  # Nominal 4C kinematic fit: J/psi -> K- pi+ gamma
  .kinematic_fit([:km, :pip, :gamma]) {
     nominal
     constrain_four_momentum
     chi2_cut 200                # loose in BOSS; tight cut (30) in ROOT
   }
  # Alternative hypothesis for pi+pi-gamma background veto
  .kinematic_fit([:pim, :pip, :gamma]) {
     constrain_four_momentum
   }
  # Alternative hypothesis for K+K-gamma background veto
  .kinematic_fit([:km, :kp, :gamma]) {
     constrain_four_momentum
   }
  # Alternative hypothesis for pi+pi-gammagamma background veto
  .kinematic_fit([:pim, :pip, :gamma, :gamma]) {
     constrain_four_momentum
   }

alg_modeI
  .note(:radiative_photon,
        "Radiative photon with maximum energy selected; N_extra_trk == 0 required in ROOT.")
  .note(:E_over_p_veto,
        "E/p < 0.8 applied for charged pion candidate in ROOT to suppress radiative Bhabha background.")
  .note(:momentum_cuts,
        "p_K < 1.25 GeV/c, p_pi < 1.35 GeV/c applied in ROOT for K/pi mis-ID suppression.")
  .note(:alternative_fit_veto,
        "chi2(pi+pi-gamma) > 85, chi2(K+K-gamma) > 50, chi2(pi+pi-gammagamma) > 35 applied in ROOT.")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_modeI])

# ==============================================================================
# Mode II: J/psi -> gamma D0, D0 -> K- pi+ pi0
# ==============================================================================
alg_modeII = Algorithm.new("JpsiGammaD0_Kpipi0")
alg_modeII.set_header(["JpsiGammaD0Alg/JpsiGammaD0.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = build_common_selection
  .pid(method: :probability) {
     prob_cut 0.001
     identify :kaon, against: [:pion]
     identify :pion, against: [:kaon]
     nkp "==1"
     nkm "==0"
     npip "==1"
     npim "==0"
   }
  # Select pi0 candidates via kinematic fit
  .kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 200                # loose; tight cut (20) in ROOT
   }
  # Nominal 5C kinematic fit: J/psi -> K- pi+ gamma gamma gamma with pi0 mass constraint
  .kinematic_fit([:km, :pip, :gamma, :gamma, :gamma]) {
     nominal
     constrain_four_momentum
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 200                # loose in BOSS; tight cut (30) in ROOT
   }
  # Alternative hypothesis for pi+pi-pi0 gamma background veto
  .kinematic_fit([:pim, :pip, :gamma, :gamma, :gamma]) {
     constrain_four_momentum
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
   }
  # Alternative hypothesis for K+K-pi0 gamma background veto
  .kinematic_fit([:km, :kp, :gamma, :gamma, :gamma]) {
     constrain_four_momentum
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
   }
  # Alternative hypothesis for pi+pi-pi0 gamma gamma background veto
  .kinematic_fit([:pim, :pip, :gamma, :gamma, :gamma, :gamma]) {
     constrain_four_momentum
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
   }

alg_modeII
  .note(:radiative_photon,
        "Radiative photon with maximum energy selected; pi0 mass window [0.115, 0.150] GeV/c^2 in ROOT.")
  .note(:pi0_1C_fit,
        "1C kinematic fit constraining gamma gamma to pi0 mass; chi2 < 20 in ROOT.")
  .note(:Ks_veto,
        "M_recoil(K-pi+) window [0.32, 0.66] GeV/c^2 to veto K_S^0 -> pi0 pi0 background (ROOT).")
  .note(:alternative_fit_veto,
        "chi2(pi+pi-pi0 gamma) > 40, chi2(K+K-pi0 gamma) > 40, chi2(pi+pi-pi0 gamma gamma) > 95 in ROOT.")
  .note(:omega_veto,
        "M(gamma pi0) window [0.68, 0.92] GeV/c^2 to veto omega background (ROOT).")
  .note(:K_decay_veto,
        "M(pi+ pi0) window [0.45, 0.53] GeV/c^2 to veto K+ -> pi+ pi0 background (ROOT).")
  .note(:MVA_KL_veto,
        "MVA (BDTG) on shower shape to veto K_L^0 background, applied in ROOT.")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_modeII])

# ==============================================================================
# Mode III: J/psi -> gamma D0, D0 -> K- pi+ pi+ pi-
# ==============================================================================
alg_modeIII = Algorithm.new("JpsiGammaD0_Kpipipi")
alg_modeIII.set_header(["JpsiGammaD0Alg/JpsiGammaD0.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeIII = build_common_selection
  .pid(method: :probability) {
     prob_cut 0.001
     identify :kaon, against: [:pion]
     identify :pion, against: [:kaon]
     nkp "==1"
     nkm "==0"
     npip "==2"
     npim "==1"
   }
  # Nominal 4C kinematic fit: J/psi -> K- pi+ pi+ pi- gamma
  .kinematic_fit([:km, :pip, :pip, :pim, :gamma]) {
     nominal
     constrain_four_momentum
     chi2_cut 200                # loose in BOSS; tight cut (30) in ROOT
   }
  # Alternative hypothesis for pi+pi-pi+pi-gamma background veto
  .kinematic_fit([:pim, :pip, :pip, :pim, :gamma]) {
     constrain_four_momentum
   }

alg_modeIII
  .note(:radiative_photon,
        "Radiative photon with maximum energy selected; N_extra_trk == 0 required in ROOT.")
  .note(:Ks_veto,
        "M(pi+pi-) window [0.46, 0.55] GeV/c^2 to veto K_S^0 background (ROOT).")
  .note(:alternative_fit_veto,
        "chi2(pi+pi-pi+pi-gamma) > 430 applied in ROOT.")
  .note(:K_decay_veto,
        "M(pi+pi+pi-) window [0.47, 0.52] and M_recoil(K-pi+pi-) window [0.42, 0.52] GeV/c^2 for K decay veto (ROOT).")
  .note(:Dalitz_veto,
        "M(e+e-gamma) window [0.09, 0.19] GeV/c^2 for pi0 Dalitz decay veto (ROOT).")

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)
root_files_modeIII = alg_modeIII.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_modeIII])