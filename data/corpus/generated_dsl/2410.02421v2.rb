# =====================================================================
# BOSS DSL — search for lepton-number-violating D_s+ -> h- h0 e+ e+
# Single-tag method: tag (recoil) side D_s- -> K+ K- pi-
# c.m. energies 4.128 - 4.226 GeV (8 data sets, 7.33 fb^-1)
# =====================================================================

### --------------------------- Datasets --------------------------- ###
# Eight real data sets (4.128 / 4.157 / 4.178 / 4.189 / 4.199 / 4.209 / 4.219 / 4.226 GeV)
data_4130 = DatasetManager.real_data.find("705_4130")   # 4.128 GeV
data_4160 = DatasetManager.real_data.find("705_4160")   # 4.157 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

# Matching inclusive MC samples
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

data_all  = [data_4130, data_4160, data_4180, data_4190,
             data_4200, data_4210, data_4220, data_4230]
incMC_all = [incMC_4130, incMC_4160, incMC_4180, incMC_4190,
             incMC_4200, incMC_4210, incMC_4220, incMC_4230]

### ------------------------- Decay cards -------------------------- ###
# Mode 1: D_s+ -> phi pi- e+ e+ , phi -> K+ K-
decay_card_phi_pi_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 phi pi- e+ e+                          PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                                  VSS;
    Enddecay

    End
DECAYCARD

# Mode 2: D_s+ -> phi K- e+ e+ , phi -> K+ K-
decay_card_phi_k_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 phi K- e+ e+                           PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                                  VSS;
    Enddecay

    End
DECAYCARD

# Mode 3: D_s+ -> K_S0 pi- e+ e+ , K_S0 -> pi+ pi-
decay_card_ks_pi_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 K_S0 pi- e+ e+                         PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                                PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4: D_s+ -> K_S0 K- e+ e+ , K_S0 -> pi+ pi-
decay_card_ks_k_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 K_S0 K- e+ e+                          PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                                PHSP;
    Enddecay

    End
DECAYCARD

# Mode 5: D_s+ -> pi- pi0 e+ e+ , pi0 -> gamma gamma
decay_card_pi_pi0_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 pi- pi0 e+ e+                          PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

# Mode 6: D_s+ -> K- pi0 e+ e+ , pi0 -> gamma gamma
decay_card_k_pi0_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s-                              PHSP;
    Enddecay

    Decay D_s+
    1.0000 K- pi0 e+ e+                           PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-                              PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

### ----------------- Exclusive MC (100k evt @ 4.178 GeV) ---------- ###
exMC_phi_pi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_phi_pi_ee"
  config.related_dataset = data_4180           # 4.178 GeV
  config.events          = 100_000
  config.decay_card      = decay_card_phi_pi_ee
  config.cross_section   = :default
end

exMC_phi_k_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_phi_k_ee"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_phi_k_ee
  config.cross_section   = :default
end

exMC_ks_pi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_ks_pi_ee"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_ks_pi_ee
  config.cross_section   = :default
end

exMC_ks_k_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_ks_k_ee"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_ks_k_ee
  config.cross_section   = :default
end

exMC_pi_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_pi_pi0_ee"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_pi_pi0_ee
  config.cross_section   = :default
end

exMC_k_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_lnv_k_pi0_ee"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_k_pi0_ee
  config.cross_section   = :default
end

### =====================================================================
### Mode 1 — D_s+ -> phi pi- e+ e+   (phi -> K+ K-)
### =====================================================================
alg_phi_pi_ee = TagAnalysis.new("DsToPhiPiEE")
alg_phi_pi_ee.set_header(["DsToPhiPiEEAlg/DsToPhiPiEE.h"])
             .set_constant({"ECMS" => [:double, 4.178]})
             .with_decay_card(decay_card_phi_pi_ee)

# Tag (recoil) side: D_s- -> K+ K- pi-  (single tag; charm = -1 pins D_s-)
alg_phi_pi_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# Signal side: h- = pi-, h0 = phi (-> K+ K-), two positrons; net charge +1
alg_phi_pi_ee.signal_side do |s|
  s.charged(kp: 1, km: 1, pim: 1, ep: 2)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_phi_pi_ee.fit do |f|
  f.constrain_four_momentum                                    # 4C fit (tag + signal)
  f.invariant_mass_of(:kp, :km).between(1.00, 1.05)            # phi mass window
  f.chi2_cut 200
end

alg_phi_pi_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm; DTagTool's own track-quality flags are used as-is")
alg_phi_pi_ee.note(:pid_correction_method,
  "the two positrons must pass the electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut; this tighter selector replaces the fixed SimplePIDSvc threshold used by the tag-layer charged(ep: ...) key")
alg_phi_pi_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV")

alg_phi_pi_ee.apply

### =====================================================================
### Mode 2 — D_s+ -> phi K- e+ e+   (phi -> K+ K-)
### =====================================================================
alg_phi_k_ee = TagAnalysis.new("DsToPhiKEE")
alg_phi_k_ee.set_header(["DsToPhiKEEAlg/DsToPhiKEE.h"])
            .set_constant({"ECMS" => [:double, 4.178]})
            .with_decay_card(decay_card_phi_k_ee)

alg_phi_k_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# h- = K-, h0 = phi (-> K+ K-), two positrons
alg_phi_k_ee.signal_side do |s|
  s.charged(kp: 1, km: 2, ep: 2)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_phi_k_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :km).between(1.00, 1.05)            # phi mass window
  f.chi2_cut 200
end

alg_phi_k_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm")
alg_phi_k_ee.note(:pid_correction_method,
  "electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut; in this mode only ONE positron is required to pass it")
alg_phi_k_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV; the M_rec-DeltaM signal region (Punzi-FOM optimised per energy point) is not applied to this mode")

alg_phi_k_ee.apply

### =====================================================================
### Mode 3 — D_s+ -> K_S0 pi- e+ e+   (K_S0 -> pi+ pi-)
### =====================================================================
alg_ks_pi_ee = TagAnalysis.new("DsToKsPiEE")
alg_ks_pi_ee.set_header(["DsToKsPiEEAlg/DsToKsPiEE.h"])
            .set_constant({"ECMS" => [:double, 4.178]})
            .with_decay_card(decay_card_ks_pi_ee)

alg_ks_pi_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# h- = pi-, h0 = K_S0 (-> pi+ pi-), two positrons
alg_ks_pi_ee.signal_side do |s|
  s.charged(pip: 1, pim: 2, ep: 2)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_ks_pi_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)        # K_S0 mass window
  f.chi2_cut 200
end

alg_ks_pi_ee.note(:ks_reconstruction,
  "K_S0 rebuilt by a secondary-vertex fit of the pi+ pi- pair requiring |Vz| < 20 cm and L/sigma_L > 2; the mass window M(pi+ pi-) in [0.487, 0.511] GeV/c^2 is carried in the 4C fit via invariant_mass_of(:pip,:pim).between")
alg_ks_pi_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm")
alg_ks_pi_ee.note(:pid_correction_method,
  "the two positrons must pass the electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut")
alg_ks_pi_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV")

alg_ks_pi_ee.apply

### =====================================================================
### Mode 4 — D_s+ -> K_S0 K- e+ e+   (K_S0 -> pi+ pi-)
### =====================================================================
alg_ks_k_ee = TagAnalysis.new("DsToKsKEE")
alg_ks_k_ee.set_header(["DsToKsKEEAlg/DsToKsKEE.h"])
           .set_constant({"ECMS" => [:double, 4.178]})
           .with_decay_card(decay_card_ks_k_ee)

alg_ks_k_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# h- = K-, h0 = K_S0 (-> pi+ pi-), two positrons
alg_ks_k_ee.signal_side do |s|
  s.charged(pip: 1, pim: 1, km: 1, ep: 2)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_ks_k_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)        # K_S0 mass window
  f.chi2_cut 200
end

alg_ks_k_ee.note(:ks_reconstruction,
  "K_S0 rebuilt by a secondary-vertex fit of the pi+ pi- pair requiring |Vz| < 20 cm and L/sigma_L > 2; the mass window M(pi+ pi-) in [0.487, 0.511] GeV/c^2 is carried in the 4C fit via invariant_mass_of(:pip,:pim).between")
alg_ks_k_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm")
alg_ks_k_ee.note(:pid_correction_method,
  "the two positrons must pass the electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut")
alg_ks_k_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV")

alg_ks_k_ee.apply

### =====================================================================
### Mode 5 — D_s+ -> pi- pi0 e+ e+   (pi0 -> gamma gamma)
### =====================================================================
alg_pi_pi0_ee = TagAnalysis.new("DsToPiPi0EE")
alg_pi_pi0_ee.set_header(["DsToPiPi0EEAlg/DsToPiPi0EE.h"])
             .set_constant({"ECMS" => [:double, 4.178]})
             .with_decay_card(decay_card_pi_pi0_ee)

alg_pi_pi0_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# h- = pi-, h0 = pi0 (-> gamma gamma), two positrons
alg_pi_pi0_ee.signal_side do |s|
  s.charged(pim: 1, ep: 2)
  s.require_charge 1
  s.photons 2                                  # at least two signal photons
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_pi_pi0_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 1C pi0 constraint
  f.chi2_cut 200
end

alg_pi_pi0_ee.note(:pi0_reconstruction,
  "pi0 rebuilt by a 1C Kalman fit of the two signal photons to the pi0 mass with chi2 < 50 and M(gamma gamma) in [0.115, 0.150] GeV/c^2; the mass constraint is carried in the 4C fit via invariant_mass_of(:gamma,:gamma).constrain_to_nominal_mass_of(:pi0)")
alg_pi_pi0_ee.note(:photon_selection,
  "signal photons: E > 25 MeV (barrel) / > 50 MeV (endcap), TDC 0-700 ns, > 10 deg from any charged track; the barrel/endcap energy split and the TDC window are handled by DTagTool isGoodShower")
alg_pi_pi0_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm")
alg_pi_pi0_ee.note(:pid_correction_method,
  "the two positrons must pass the electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut")
alg_pi_pi0_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV")

alg_pi_pi0_ee.apply

### =====================================================================
### Mode 6 — D_s+ -> K- pi0 e+ e+   (pi0 -> gamma gamma)
### =====================================================================
alg_k_pi0_ee = TagAnalysis.new("DsToKPi0EE")
alg_k_pi0_ee.set_header(["DsToKPi0EEAlg/DsToKPi0EE.h"])
            .set_constant({"ECMS" => [:double, 4.178]})
            .with_decay_card(decay_card_k_pi0_ee)

alg_k_pi0_ee.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# h- = K-, h0 = pi0 (-> gamma gamma), two positrons
alg_k_pi0_ee.signal_side do |s|
  s.charged(km: 1, ep: 2)
  s.require_charge 1
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_k_pi0_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_k_pi0_ee.note(:pi0_reconstruction,
  "pi0 rebuilt by a 1C Kalman fit of the two signal photons to the pi0 mass with chi2 < 50 and M(gamma gamma) in [0.115, 0.150] GeV/c^2; the mass constraint is carried in the 4C fit via invariant_mass_of(:gamma,:gamma).constrain_to_nominal_mass_of(:pi0)")
alg_k_pi0_ee.note(:photon_selection,
  "signal photons: E > 25 MeV (barrel) / > 50 MeV (endcap), TDC 0-700 ns, > 10 deg from any charged track; the barrel/endcap energy split and the TDC window are handled by DTagTool isGoodShower")
alg_k_pi0_ee.note(:track_quality,
  "signal-side charged tracks must satisfy |cos(theta)| < 0.93, |Vz| < 10 cm, Vxy < 1 cm")
alg_k_pi0_ee.note(:pid_correction_method,
  "the two positrons must pass the electron criterion L_e/(L_e+L_pi+L_K) > 0.8 including an E/p cut")
alg_k_pi0_ee.note(:background_veto,
  "cos(theta)(e+ e-) < 0.95; cos(theta)(e+ pi-) < 0.98; L_pi/sigma_L < 3 (K_S0 veto); L/sigma_L > 2 (fake-K_S0 veto); E(pi0) > 0.17 GeV")

alg_k_pi0_ee.apply

### --------------------------- Execution --------------------------- ###
# Every mode runs on the eight data points + their inclusive MC + its own 100k signal MC.
# NOTE: the measured per-run beam energy and boost come from MeasuredEcmsSvc; the ECMS
# constant only fixes the reference/spread (signal MC sits at 4.178 GeV).
root_files_phi_pi_ee = alg_phi_pi_ee.execute_on(data_all + incMC_all + [exMC_phi_pi_ee])
root_files_phi_k_ee  = alg_phi_k_ee.execute_on(data_all + incMC_all + [exMC_phi_k_ee])
root_files_ks_pi_ee  = alg_ks_pi_ee.execute_on(data_all + incMC_all + [exMC_ks_pi_ee])
root_files_ks_k_ee   = alg_ks_k_ee.execute_on(data_all + incMC_all + [exMC_ks_k_ee])
root_files_pi_pi0_ee = alg_pi_pi0_ee.execute_on(data_all + incMC_all + [exMC_pi_pi0_ee])
root_files_k_pi0_ee  = alg_k_pi0_ee.execute_on(data_all + incMC_all + [exMC_k_pi0_ee])