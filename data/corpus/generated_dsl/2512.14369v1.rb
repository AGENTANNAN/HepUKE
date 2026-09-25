# ============================================================
# ψ(3686) → γ χ_c1, χ_c1 → φφη / φφη' / φ K+K− η  (BOSS part)
# ============================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# --- Decay cards (EvtGen syntax) ---
# Mode 1: χ_c1 → φ φ η,  η → γγ
decay_card_phipheta = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi phi eta PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2 (I): χ_c1 → φ φ η',  η' → π+ π- γ
decay_card_etap_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi phi eta' PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta'
    1.000 pi+ pi- gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3 (II): χ_c1 → φ φ η',  η' → π+ π- η,  η → γγ
decay_card_etap_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi phi eta' PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4: χ_c1 → φ K+ K- η,  η → γγ
decay_card_phiKKeta = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi K+ K- eta PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples (1,000,000 events each) ---
exMC_phipheta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_phipheta"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_phipheta
    config.cross_section   = :default
end

exMC_etap_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_phiphi_etap_modeI"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_etap_modeI
    config.cross_section   = :default
end

exMC_etap_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_phiphi_etap_modeII"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_etap_modeII
    config.cross_section   = :default
end

exMC_phiKKeta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_phiKKeta"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_phiKKeta
    config.cross_section   = :default
end

# ============================================================
# Event selection (BOSS)
# ============================================================

# ------------------------- Mode 1: φφη -------------------------
alg_phipheta = Algorithm.new("Chic1PhiPhiEta")
alg_phipheta.set_header(["Chic1PhiPhiEtaAlg/Chic1PhiPhiEta.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

sel_phipheta = Selection.new
  .select_track {                # 2 positive + 2 negative tracks (all kaons), net charge 0
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==2"
      nChrn     "==2"
      nNet      "==0"
  }
  .select_photon {               # ≥3 photons (barrel 25 MeV / endcap 50 MeV, TDC 0–14, >10° to tracks)
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
      nGam               ">=3"
  }
  .pid(method: :probability) {   # K/π separation; all tracks must be kaons (no pions)
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=2"
      nkm ">=2"
  }
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # γγ → η mass constraint
      chi2_cut 200              # loose; tight χ²<38 applied at ROOT
  }
alg_phipheta.with_decay_card(decay_card_phipheta).apply(sel_phipheta)

# ------------------ Mode 2 (I): φφη', η'→π+π−γ ------------------
alg_etap_modeI = Algorithm.new("Chic1PhiPhiEtapModeI")
alg_etap_modeI.set_header(["Chic1PhiPhiEtapModeIAlg/Chic1PhiPhiEtapModeI.h"])
              .set_constant({"ECMS" => [:double, 3.686]})

sel_etap_modeI = Selection.new
  .select_track {                # 3 positive + 3 negative tracks (2K+2K− + π+π−), net charge 0
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==3"
      nChrn     "==3"
      nNet      "==0"
  }
  .select_photon {               # ≥2 photons
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
      nGam               ">=2"
  }
  .pid(method: :probability) {   # K/π separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp ">=2"
      nkm ">=2"
      npip ">=1"
      npim ">=1"
  }
  .kinematic_fit([:kp, :km, :kp, :km, :pip, :pim, :gamma]) {
      nominal
      constrain_four_momentum
      chi2_cut 200              # loose; tight cut at ROOT
  }
alg_etap_modeI.with_decay_card(decay_card_etap_modeI).apply(sel_etap_modeI)

# -------------- Mode 3 (II): φφη', η'→π+π−η, η→γγ --------------
alg_etap_modeII = Algorithm.new("Chic1PhiPhiEtapModeII")
alg_etap_modeII.set_header(["Chic1PhiPhiEtapModeIIAlg/Chic1PhiPhiEtapModeII.h"])
               .set_constant({"ECMS" => [:double, 3.686]})

sel_etap_modeII = Selection.new
  .select_track {                # 3 positive + 3 negative tracks, net charge 0
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==3"
      nChrn     "==3"
      nNet      "==0"
  }
  .select_photon {               # ≥3 photons
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
      nGam               ">=3"
  }
  .pid(method: :probability) {   # K/π separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp ">=2"
      nkm ">=2"
      npip ">=1"
      npim ">=1"
  }
  .kinematic_fit([:kp, :km, :kp, :km, :pip, :pim, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # γγ → η mass constraint
      chi2_cut 200              # loose; tight cut at ROOT
  }
alg_etap_modeII.with_decay_card(decay_card_etap_modeII).apply(sel_etap_modeII)

# --------------------- Mode 4: φ K+K− η ---------------------
alg_phiKKeta = Algorithm.new("Chic1PhiKKeta")
alg_phiKKeta.set_header(["Chic1PhiKKetaAlg/Chic1PhiKKeta.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

sel_phiKKeta = Selection.new
  .select_track {                # 2 positive + 2 negative tracks (K+K− from φ and K+K− direct), net charge 0
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==2"
      nChrn     "==2"
      nNet      "==0"
  }
  .select_photon {               # ≥3 photons
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
      nGam               ">=3"
  }
  .pid(method: :probability) {   # K/π separation; no pion requirement
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=2"
      nkm ">=2"
  }
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # γγ → η mass constraint
      chi2_cut 200              # loose; tight cut at ROOT
  }
alg_phiKKeta.with_decay_card(decay_card_phiKKeta).apply(sel_phiKKeta)

# ============================================================
# Execute on real data + inclusive MC + exclusive MC
# ============================================================
root_files_phipheta = alg_phipheta.execute_on([psip_data, psip_incMC, exMC_phipheta])
root_files_etap_modeI = alg_etap_modeI.execute_on([psip_data, psip_incMC, exMC_etap_modeI])
root_files_etap_modeII = alg_etap_modeII.execute_on([psip_data, psip_incMC, exMC_etap_modeII])
root_files_phiKKeta = alg_phiKKeta.execute_on([psip_data, psip_incMC, exMC_phiKKeta])

# NOTE: final φ / η / η' mass windows, the π0 veto (Mode I) and the J/ψ veto
# (φ K+K− η) are applied at ROOT level — out of scope of the BOSS spec.