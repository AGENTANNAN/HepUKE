# BESIII: "Studies of Charmonium at BESIII" (R.-G. Ping) — proceedings summarizing
# six independent charmonium analyses. Each signal process has a different final
# state and selection, so each gets its own Algorithm + Selection chain.
#
# Data: psi(3686) 106 M, J/psi 1.31 B, psi(3770) 2.9 fb^-1, 42 pb^-1 continuum at 3.65 GeV.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")        # psi(3770), 2.92 fb^-1
data_3686  = DatasetManager.real_data.find("709_3686")        # psi(3686), 106 M events
data_3097  = DatasetManager.real_data.find("708_3097")        # J/psi, 1.31 B events
data_3650  = DatasetManager.real_data.find("709_3650")        # 3.65 GeV continuum (QED background)

incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")

# =============================================================================
# ANALYSIS 1 — Search for psi(3770) exclusive baryonic decays
# =============================================================================
decay_card_LLpipi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi+ pi-       PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_LLpi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi0           PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_LLeta = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 eta           PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_SigmaSigma = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma+ anti-Sigma-                 PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0                            PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- pi0                       PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                       PHSP;
  Enddecay
  End
DECAYCARD

decay_card_Sigma0Sigma0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma0 anti-Sigma0                 PHSP;
  Enddecay
  Decay Sigma0
  1.0000 Lambda0 gamma                      PHSP;
  Enddecay
  Decay anti-Sigma0
  1.0000 anti-Lambda0 gamma                 PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_XiXi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi- anti-Xi+                       PHSP;
  Enddecay
  Decay Xi-
  1.0000 Lambda0 pi-                        PHSP;
  Enddecay
  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+                   PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_Xi0Xi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi0 anti-Xi0                       PHSP;
  Enddecay
  Decay Xi0
  1.0000 Lambda0 pi0                        PHSP;
  Enddecay
  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0                   PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-                             PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+                        PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

exMC_modes = [
  ["baryonic_LLpipi",          decay_card_LLpipi],
  ["baryonic_LLpi0",           decay_card_LLpi0],
  ["baryonic_LLeta",           decay_card_LLeta],
  ["baryonic_SigmaSigma",      decay_card_SigmaSigma],
  ["baryonic_Sigma0Sigma0",    decay_card_Sigma0Sigma0],
  ["baryonic_XiXi",            decay_card_XiXi],
  ["baryonic_Xi0Xi0",          decay_card_Xi0Xi0],
].map do |name, card|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = name
    config.related_dataset = data_3773
    config.events          = 200000
    config.decay_card      = card
    config.cross_section   = :default
  end
end

# Common track + photon selection for all baryonic channels
baryonic_common = Selection.new
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz        10.0        # closest approach to IP < 10 cm along the beam
    Vr        1.0         # closest approach to IP < 1 cm transverse
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0    # >= 20 deg from any charged track
    energyThreshold_b 0.025   # E > 25 MeV in the barrel
    energyThreshold_e 0.050   # E > 50 MeV in the end cap
  }

# --- Mode: psi(3770) -> Lambda anti-Lambda pi+ pi- ---
alg_LLpipi = Algorithm.new("BaryonicLLpipi")
alg_LLpipi.set_header(["BaryonicLLpipiAlg/BaryonicLLpipi.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
sel_LLpipi = baryonic_common.dup
  .select_track { nChrp "==3"; nChrn "==3" }       # p, pi+ (Lambda_bar), pi+  /  pi-, pbar, pi-
  .select_photon  { nGam "==0" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])        # keep protons for the Lambda vertex fits
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .assign({ chrgp: :pip, chrgn: :pim })            # remaining tracks are the two pions
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_LLpipi.note(:baryon_vertex_fit,
                "Lambda0 -> p pi- and anti-Lambda0 -> anti-p pi+ are reconstructed with a " \
                "secondary vertex fit from oppositely charged p and pi tracks; the " \
                "|m(p pi) - m(Lambda)| < 0.02 GeV/c^2 requirement is applied downstream in ROOT.")
             .note(:qed_background,
                   "QED backgrounds are estimated from data taken at sqrt(s) = 3.542, 3.554, " \
                   "3.561, 3.600 and 3.650 GeV; ISR backgrounds (e+e- -> gamma psi(3686), " \
                   "gamma psi(3770)) are estimated with MC simulation and subtracted before " \
                   "the 90% C.L. upper limits are extracted.")
             .with_decay_card(decay_card_LLpipi).apply(sel_LLpipi)
alg_LLpipi.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[0..0])

# --- Mode: psi(3770) -> Lambda anti-Lambda pi0 ---
alg_LLpi0 = Algorithm.new("BaryonicLLpi0")
alg_LLpi0.set_header(["BaryonicLLpi0Alg/BaryonicLLpi0.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
sel_LLpi0 = baryonic_common.dup
  .select_track { nChrp "==2"; nChrn "==2" }       # p pi- / pbar pi+
  .select_photon  { nGam ">=2" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_LLpi0.note(:pi0_reconstruction,
               "pi0 candidates are built from photon pairs with 0.115 < m(gamma gamma) < " \
               "0.150 GeV/c^2 and enter the fit as a single mass-constrained participant.")
        .with_decay_card(decay_card_LLpi0).apply(sel_LLpi0)
alg_LLpi0.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[1..1])

# --- Mode: psi(3770) -> Lambda anti-Lambda eta ---
alg_LLeta = Algorithm.new("BaryonicLLeta")
alg_LLeta.set_header(["BaryonicLLetaAlg/BaryonicLLeta.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
sel_LLeta = baryonic_common.dup
  .select_track { nChrp "==2"; nChrn "==2" }
  .select_photon  { nGam ">=2" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # 0.50 < m(gg) < 0.58
    chi2_cut 25
    neta "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_LLeta.note(:eta_reconstruction,
               "eta candidates are built from photon pairs with 0.50 < m(gamma gamma) < " \
               "0.58 GeV/c^2 and enter the fit as a single mass-constrained participant.")
        .with_decay_card(decay_card_LLeta).apply(sel_LLeta)
alg_LLeta.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[2..2])

# --- Mode: psi(3770) -> Sigma+ anti-Sigma-, Sigma+ -> p pi0 ---
alg_SigmaSigma = Algorithm.new("BaryonicSigmaSigma")
alg_SigmaSigma.set_header(["BaryonicSigmaSigmaAlg/BaryonicSigmaSigma.h"])
              .set_constant({ "ECMS" => [:double, 3.773] })
sel_SigmaSigma = baryonic_common.dup
  .select_track { nChrp "==1"; nChrn "==1" }       # p / pbar
  .select_photon  { nGam ">=4" }                   # 2 pi0 -> 4 photons
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_SigmaSigma.with_decay_card(decay_card_SigmaSigma).apply(sel_SigmaSigma)
alg_SigmaSigma.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[3..3])

# --- Mode: psi(3770) -> Sigma0 anti-Sigma0, Sigma0 -> Lambda gamma ---
alg_Sigma0Sigma0 = Algorithm.new("BaryonicSigma0Sigma0")
alg_Sigma0Sigma0.set_header(["BaryonicSigma0Sigma0Alg/BaryonicSigma0Sigma0.h"])
                .set_constant({ "ECMS" => [:double, 3.773] })
sel_Sigma0Sigma0 = baryonic_common.dup
  .select_track { nChrp "==2"; nChrn "==2" }       # p pi- / pbar pi+
  .select_photon  { nGam ">=2" }                   # one photon per Sigma0
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_Sigma0Sigma0.with_decay_card(decay_card_Sigma0Sigma0).apply(sel_Sigma0Sigma0)
alg_Sigma0Sigma0.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[4..4])

# --- Mode: psi(3770) -> Xi- anti-Xi+, Xi- -> Lambda pi- ---
alg_XiXi = Algorithm.new("BaryonicXiXi")
alg_XiXi.set_header(["BaryonicXiXiAlg/BaryonicXiXi.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
sel_XiXi = baryonic_common.dup
  .select_track { nChrp "==3"; nChrn "==3" }       # p pi- pi- / pbar pi+ pi+
  .select_photon  { nGam "==0" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_XiXi.note(:cascade_daughters,
              "The extra charged pion of Xi- -> Lambda pi- (anti-Xi+ -> anti-Lambda pi+) is " \
              "taken from the remaining tracks after the Lambda / anti-Lambda secondary vertex " \
              "fits; the intermediate Xi mass window is applied downstream in ROOT.")
      .with_decay_card(decay_card_XiXi).apply(sel_XiXi)
alg_XiXi.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[5..5])

# --- Mode: psi(3770) -> Xi0 anti-Xi0, Xi0 -> Lambda pi0 ---
alg_Xi0Xi0 = Algorithm.new("BaryonicXi0Xi0")
alg_Xi0Xi0.set_header(["BaryonicXi0Xi0Alg/BaryonicXi0Xi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
sel_Xi0Xi0 = baryonic_common.dup
  .select_track { nChrp "==2"; nChrn "==2" }
  .select_photon  { nGam ">=4" }                   # 2 pi0 -> 4 photons
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"; nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_Xi0Xi0.with_decay_card(decay_card_Xi0Xi0).apply(sel_Xi0Xi0)
alg_Xi0Xi0.execute_on([data_3773, incMC_3773, data_3650] + exMC_modes[6..6])

# =============================================================================
# ANALYSIS 2 — psi(3770) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+- pi-+
# =============================================================================
decay_card_etac2S = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma eta_c(2S)                    PHSP;
  Enddecay
  Decay eta_c(2S)
  1.0000 K_S0 K+ pi-                        PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay
  End
DECAYCARD

exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_to_gamma_etac2S"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

alg_etac2S = Algorithm.new("Psi3770ToGammaEtac2S")
alg_etac2S.set_header(["Psi3770ToGammaEtac2SAlg/Psi3770ToGammaEtac2S.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
sel_etac2S = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"     # K_S0 daughter pi+ plus K+ / pi+
    nChrn     ">=2"     # K_S0 daughter pi- plus K- / pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"   # the energetic radiative photon
  }
  .assign({ chrgp: :pip, chrgn: :pim })   # K_S0 daughters are treated as pions (no PID)
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]    # the single charged kaon of K_S0 K+- pi-+
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_etac2S.note(:kS0_reconstruction,
                "K_S0 -> pi+ pi- is reconstructed from pairs of oppositely charged tracks " \
                "with a vertex-constrained fit and 0.487 < m(pi+ pi-) < 0.511 GeV/c^2; the " \
                "K_S0 enters the kinematic fit as one virtual particle.")
          .note(:event_topology,
                "Exactly one K_S0, one charged kaon and one charged pion are required in the " \
                "event (charge-conjugation inclusive: K_S0 K+ pi- or K_S0 K- pi+); the " \
                "kaon/pion assignment is resolved by the kinematic fit and the residual " \
                "multiplicity requirement is applied downstream in ROOT.")
          .note(:upper_limit,
                "No significant eta_c(2S) signal is observed. The 90% C.L. upper limits are " \
                "B(psi(3770) -> gamma eta_c(2S) -> gamma K_S0 K+- pi-+) < 5.6e-6 and " \
                "B(psi(3770) -> gamma eta_c(2S)) < 2.0e-3.")
          .with_decay_card(decay_card_etac2S).apply(sel_etac2S)
alg_etac2S.execute_on([data_3773, incMC_3773, exMC_etac2S])

# =============================================================================
# ANALYSIS 3 — psi(3770) -> gamma chi_c1 (and search for gamma chi_c2),
#              chi_cJ -> gamma J/psi, J/psi -> l+ l-
# =============================================================================
decay_card_chic1_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1                       PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi                        PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e-                              PHOTOS PHSP;
  Enddecay
  End
DECAYCARD

decay_card_chic1_mumu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1                       PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi                        PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-                            PHOTOS PHSP;
  Enddecay
  End
DECAYCARD

decay_card_chic2_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c2                       PHSP;
  Enddecay
  Decay chi_c2
  1.0000 gamma J/psi                        PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e-                              PHOTOS PHSP;
  Enddecay
  End
DECAYCARD

exMC_chic1_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_to_gamma_chic1_to_gamma_ee"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_chic1_ee
  config.cross_section   = :default
end

exMC_chic1_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_to_gamma_chic1_to_gamma_mumu"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_chic1_mumu
  config.cross_section   = :default
end

exMC_chic2_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_to_gamma_chic2_to_gamma_ee"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_chic2_ee
  config.cross_section   = :default
end

# --- Electron channel: J/psi -> e+ e- ---
alg_chic_ee = Algorithm.new("Psi3770ToGammaChicEE")
alg_chic_ee.set_header(["Psi3770ToGammaChicEEAlg/Psi3770ToGammaChicEE.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
sel_chic_ee = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"     # e+
    nChrn     "==1"     # e-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              "==2"   # one radiative photon from psi(3770), one from chi_cJ -> gamma J/psi
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"; nlm "==1"
  }
  .kalman_kinematic_fit([:lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 25
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_chic_ee.note(:electron_pid,
                 "Electron identification requires L_e/(L_e + L_pi + L_K) > 0.8 from combined " \
                 "dE/dx and TOF likelihoods together with E/p in [0.8, 1.2]; bremsstrahlung " \
                 "recovery adds EMC clusters within 5 deg of the track.")
            .note(:results,
                  "B(psi(3770) -> gamma chi_c1) = (2.48 +- 0.15 +- 0.23)e-3; the 90% C.L. " \
                  "upper limit for chi_c2 is B(psi(3770) -> gamma chi_c2) < 0.64e-3. The " \
                  "chi_cJ invariant mass is formed from the energetic photon and the J/psi.")
            .with_decay_card(decay_card_chic1_ee).apply(sel_chic_ee)
alg_chic_ee.execute_on([data_3773, incMC_3773, exMC_chic1_ee, exMC_chic2_ee])

# --- Muon channel: J/psi -> mu+ mu- ---
alg_chic_mumu = Algorithm.new("Psi3770ToGammaChicMuMu")
alg_chic_mumu.set_header(["Psi3770ToGammaChicMuMuAlg/Psi3770ToGammaChicMuMu.h"])
             .set_constant({ "ECMS" => [:double, 3.773] })
sel_chic_mumu = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              "==2"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"; nlm "==1"
  }
  .kalman_kinematic_fit([:lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 25
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_chic_mumu.note(:muon_pid,
                   "The muon channel uses the standard BESIII muon identification from MUC " \
                   "penetration depth combined with the high-momentum lepton selector.")
             .with_decay_card(decay_card_chic1_mumu).apply(sel_chic_mumu)
alg_chic_mumu.execute_on([data_3773, incMC_3773, exMC_chic1_mumu])

# =============================================================================
# ANALYSIS 4 — Isospin-violating chi_c0,2 -> pi0 eta_c, eta_c -> K_S0 K+- pi-+
#              in psi(3686) -> gamma chi_c0,2
# =============================================================================
decay_card_chic0_etac = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0                       PHSP;
  Enddecay
  Decay chi_c0
  1.0000 pi0 eta_c                          PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi-                        PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_chic2_etac = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2                       PHSP;
  Enddecay
  Decay chi_c2
  1.0000 pi0 eta_c                          PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi-                        PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

exMC_chic0_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic0_to_pi0_etac"
  config.related_dataset = data_3686
  config.events          = 200000
  config.decay_card      = decay_card_chic0_etac
  config.cross_section   = :default
end

exMC_chic2_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic2_to_pi0_etac"
  config.related_dataset = data_3686
  config.events          = 200000
  config.decay_card      = decay_card_chic2_etac
  config.cross_section   = :default
end

# --- Mode: chi_c0 -> pi0 eta_c ---
alg_chic0_etac = Algorithm.new("Chic0ToPi0Etac")
alg_chic0_etac.set_header(["Chic0ToPi0EtacAlg/Chic0ToPi0Etac.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })
sel_chic0_etac = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"   # 2 photons from pi0 + the radiative photon
  }
  .assign({ chrgp: :pip, chrgn: :pim })
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==1"
  }
  .kinematic_fit([:gamma, :pi0, :K_S0, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_chic0_etac.note(:kS0_and_eta_c,
                    "K_S0 -> pi+ pi- from oppositely charged track pairs with a vertex fit and " \
                    "0.487 < m(pi+ pi-) < 0.511 GeV/c^2; the eta_c is reconstructed through " \
                    "eta_c -> K_S0 K+- pi-+ and its mass window is applied downstream.")
               .note(:upper_limit,
                     "No statistically significant signal is observed: " \
                     "B(chi_c0 -> pi0 eta_c) < 1.6e-3 at the 90% C.L.")
               .with_decay_card(decay_card_chic0_etac).apply(sel_chic0_etac)
alg_chic0_etac.execute_on([data_3686, incMC_3686, exMC_chic0_etac])

# --- Mode: chi_c2 -> pi0 eta_c ---
alg_chic2_etac = Algorithm.new("Chic2ToPi0Etac")
alg_chic2_etac.set_header(["Chic2ToPi0EtacAlg/Chic2ToPi0Etac.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })
sel_chic2_etac = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .assign({ chrgp: :pip, chrgn: :pim })
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==1"
  }
  .kinematic_fit([:gamma, :pi0, :K_S0, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_chic2_etac.note(:upper_limit,
                    "No statistically significant signal is observed: " \
                    "B(chi_c2 -> pi0 eta_c) < 3.2e-3 at the 90% C.L.")
               .with_decay_card(decay_card_chic2_etac).apply(sel_chic2_etac)
alg_chic2_etac.execute_on([data_3686, incMC_3686, exMC_chic2_etac])

# =============================================================================
# ANALYSIS 5 — C-violation searches psi(3686) -> pi+ pi- J/psi,
#              J/psi -> gamma gamma and J/psi -> gamma phi, phi -> K+ K-
# =============================================================================
decay_card_jpsi_gammagamma = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi                      PHSP;
  Enddecay
  Decay J/psi
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

decay_card_jpsi_gammaphi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi                      PHSP;
  Enddecay
  Decay J/psi
  1.0000 gamma phi                          PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K-                              VSS;
  Enddecay
  End
DECAYCARD

exMC_gammagamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pipi_jpsi_to_gammagamma"
  config.related_dataset = data_3686
  config.events          = 200000
  config.decay_card      = decay_card_jpsi_gammagamma
  config.cross_section   = :default
end

exMC_gammaphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pipi_jpsi_to_gammaphi"
  config.related_dataset = data_3686
  config.events          = 200000
  config.decay_card      = decay_card_jpsi_gammaphi
  config.cross_section   = :default
end

# --- Channel: psi(3686) -> pi+ pi- J/psi, J/psi -> gamma gamma ---
alg_gammagamma = Algorithm.new("JpsiToGammaGamma")
alg_gammagamma.set_header(["JpsiToGammaGammaAlg/JpsiToGammaGamma.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })
sel_gammagamma = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"     # pi+
    nChrn     "==1"     # pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              "==2"   # both photons from J/psi -> gamma gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"; npim "==1"
  }
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_gammagamma.note(:signal_extraction,
                    "The J/psi candidates are searched for in the recoil mass " \
                    "M_{pi+ pi-}^{rec}; the global fit to that distribution gives the " \
                    "90% C.L. upper limit B(J/psi -> gamma gamma) < 2.7e-7, one order of " \
                    "magnitude more stringent than the previous limit.")
              .with_decay_card(decay_card_jpsi_gammagamma).apply(sel_gammagamma)
alg_gammagamma.execute_on([data_3686, incMC_3686, exMC_gammagamma])

# --- Channel: psi(3686) -> pi+ pi- J/psi, J/psi -> gamma phi, phi -> K+ K- ---
alg_gammaphi = Algorithm.new("JpsiToGammaPhi")
alg_gammaphi.set_header(["JpsiToGammaPhiAlg/JpsiToGammaPhi.h"])
            .set_constant({ "ECMS" => [:double, 3.686] })
sel_gammaphi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"     # pi+ and K+
    nChrn     "==2"     # pi- and K-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"   # the radiative photon of J/psi -> gamma phi
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"; nkm "==1"
    npip "==1"; npim "==1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_gammaphi.note(:signal_extraction,
                  "The phi candidates are searched for in the K+ K- invariant mass " \
                  "distribution; no significant signal is observed and the 90% C.L. upper " \
                  "limit is B(J/psi -> gamma phi) < 1.4e-6, the first limit for this channel.")
            .with_decay_card(decay_card_jpsi_gammaphi).apply(sel_gammaphi)
alg_gammaphi.execute_on([data_3686, incMC_3686, exMC_gammaphi])

# =============================================================================
# ANALYSIS 6 — OZI-suppressed decay J/psi -> pi0 phi, phi -> K+ K-, pi0 -> gamma gamma
# =============================================================================
decay_card_jpsi_pi0phi = <<~DECAYCARD
  Decay J/psi
  1.0000 pi0 phi                            PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K-                              VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay
  End
DECAYCARD

exMC_pi0phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_pi0_phi"
  config.related_dataset = data_3097
  config.events          = 500000
  config.decay_card      = decay_card_jpsi_pi0phi
  config.cross_section   = :default
end

alg_pi0phi = Algorithm.new("JpsiToPi0Phi")
alg_pi0phi.set_header(["JpsiToPi0PhiAlg/JpsiToPi0Phi.h"])
         .set_constant({ "ECMS" => [:double, 3.097] })
sel_pi0phi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"     # K+
    nChrn     "==1"     # K-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              "==2"   # pi0 -> gamma gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"; nkm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # 0.115 < m(gg) < 0.150
    chi2_cut 25
    npi0 "==1"
  }
  .kinematic_fit([:kp, :km, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0phi.note(:kaon_pid,
               "Kaons are identified from combined dE/dx and TOF likelihoods with " \
               "L_K > L_pi; photon candidates must satisfy the EMC timing requirement and be " \
               "at least 20 deg from any charged track.")
         .note(:pi0_selection,
               "pi0 -> gamma gamma candidates require 0.115 < m(gamma gamma) < 0.150 GeV/c^2; " \
               "candidates with both photons in the end cap are rejected, and the chi2 of the " \
               "1C mass-constrained Kalman fit is required to be small.")
         .note(:signal_extraction,
               "A structure around 1.02 GeV/c^2 is observed in the K+ K- invariant mass " \
               "spectrum and interpreted as interference of J/psi -> phi pi0 with other " \
               "processes decaying to the same final state. The fit yields two solutions: " \
               "[2.94 +- 0.16 +- 0.16]e-6 and [1.24 +- 0.33 +- 0.30]e-6. Sideband " \
               "subtraction is applied to the M(K+ K-) spectrum before fitting.")
         .with_decay_card(decay_card_jpsi_pi0phi).apply(sel_pi0phi)
alg_pi0phi.execute_on([data_3097, incMC_3097, exMC_pi0phi])
