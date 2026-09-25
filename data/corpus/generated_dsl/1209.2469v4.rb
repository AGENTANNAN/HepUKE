# =====================================================================
# Dataset preparation  (J/psi(3097) resonance)
# =====================================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# ---------------------------------------------------------------------
# Decay cards (EvtGen format) — seven exclusive-MC modes
# ---------------------------------------------------------------------
# (1) J/psi -> phi eta, eta left undecayed (invisible)
decay_card_eta_inv = <<~DECAYCARD
  Decay J/psi
  1.0 phi eta PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# (2) J/psi -> phi eta', eta' left undecayed (invisible)
decay_card_etap_inv = <<~DECAYCARD
  Decay J/psi
  1.0 phi eta' PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# (3) J/psi -> phi eta, eta -> gamma gamma   (normalisation channel)
decay_card_eta_gg = <<~DECAYCARD
  Decay J/psi
  1.0 phi eta PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (4) J/psi -> phi eta', eta' -> gamma gamma  (normalisation channel)
decay_card_etap_gg = <<~DECAYCARD
  Decay J/psi
  1.0 phi eta' PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay eta'
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (5) background: J/psi -> gamma eta_c, eta_c -> K+ pi- K_L0
decay_card_bkg_etac = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ pi- K_L0 PHSP;
  Enddecay

  End
DECAYCARD

# (6) background: J/psi -> phi K_L0 K_L0
decay_card_bkg_phiKLKL = <<~DECAYCARD
  Decay J/psi
  1.0 phi K_L0 K_L0 PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# (7) background: J/psi -> phi f_0(980), f_0(980) -> K_L0 K_L0
decay_card_bkg_phif0 = <<~DECAYCARD
  Decay J/psi
  1.0 phi f_0 PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay f_0
  1.0 K_L0 K_L0 PHSP;
  Enddecay

  End
DECAYCARD

# ---------------------------------------------------------------------
# Exclusive MC samples — 200k events each, matched to the J/psi data
# ---------------------------------------------------------------------
exMC_eta_inv = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_invisible"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_eta_inv
  config.cross_section   = :default
end

exMC_etap_inv = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_etap_invisible"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_etap_inv
  config.cross_section   = :default
end

exMC_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_eta_gg
  config.cross_section   = :default
end

exMC_etap_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_etap_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_etap_gg
  config.cross_section   = :default
end

exMC_bkg_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac_KpiKL"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_etac
  config.cross_section   = :default
end

exMC_bkg_phiKLKL = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_KLKL"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_phiKLKL
  config.cross_section   = :default
end

exMC_bkg_phif0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_f0_KLKL"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_phif0
  config.cross_section   = :default
end

# =====================================================================
# Event selection (BOSS) — four independent signal chains.
# Charged-track, photon and PID criteria are identical for all chains.
# =====================================================================

# -------- Algorithms --------
alg_eta_inv = Algorithm.new("JpsiPhiEtaInv")
alg_eta_inv.set_header(["JpsiPhiEtaInvAlg/JpsiPhiEtaInv.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:photon_transition_region,
                 "photon energy threshold in the EMC barrel/endcap transition region " \
                 "(0.80 < |cos(theta)| < 0.86) is 100 MeV; select_photon only exposes " \
                 "energyThreshold_b (barrel) and energyThreshold_e (endcap), so this " \
                 "transition-region threshold is applied in the generated BOSS code")

alg_etap_inv = Algorithm.new("JpsiPhiEtaPrimeInv")
alg_etap_inv.set_header(["JpsiPhiEtaPrimeInvAlg/JpsiPhiEtaPrimeInv.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:photon_transition_region,
                  "photon energy threshold in the EMC barrel/endcap transition region " \
                  "(0.80 < |cos(theta)| < 0.86) is 100 MeV; applied in the generated BOSS code")

alg_eta_gg = Algorithm.new("JpsiPhiEtaGG")
alg_eta_gg.set_header(["JpsiPhiEtaGGAlg/JpsiPhiEtaGG.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:photon_transition_region,
                "photon energy threshold in the EMC barrel/endcap transition region " \
                "(0.80 < |cos(theta)| < 0.86) is 100 MeV; applied in the generated BOSS code")

alg_etap_gg = Algorithm.new("JpsiPhiEtaPrimeGG")
alg_etap_gg.set_header(["JpsiPhiEtaPrimeGGAlg/JpsiPhiEtaPrimeGG.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:photon_transition_region,
                 "photon energy threshold in the EMC barrel/endcap transition region " \
                 "(0.80 < |cos(theta)| < 0.86) is 100 MeV; applied in the generated BOSS code")

# -------- Common bases --------
# Common base for the two invisible (eta/eta' -> invisible) chains:
# two charged tracks (1 K+ + 1 K-), no photon-count requirement.
base_invisible = Selection.new
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93
    Vz        10.0          # |Vz| < 10 cm
    Vr        1.0           # Vr < 1 cm
    nChrp     "==1"         # exactly one positive track
    nChrn     "==1"         # exactly one negative track
    nNet      "==0"         # net charge zero
  }
  .select_photon {
    tdc_emc_start     0     # TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0  # isolated by > 10 deg from all charged tracks
    energyThreshold_b 0.025 # 25 MeV (barrel)
    energyThreshold_e 0.050 # 50 MeV (endcap)
  }

# Common base for the two gamma-gamma (eta/eta' -> gamma gamma) chains:
# same track / photon quality cuts, plus at least two photons.
base_gg = Selection.new
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
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2" # at least two photons for the gamma-gamma modes
  }

# -------- Chain A: J/psi -> phi eta, eta invisible --------
sel_eta_inv = base_invisible.dup
  .pid(method: :probability) {
    prob_cut 0.001                                # PID probability > 0.001
    identify :kaon, against: [:pion, :proton]     # K+ and K- (charge-conjugation shorthand)
    nkp "==1"                                     # exactly one K+
    nkm "==1"                                     # exactly one K-
  }
  # Reconstruct only the phi -> K+ K- tag (recIDs: 1 = phi, 3 = K+, 4 = K-);
  # the invisible eta is inferred from the recoil against the phi.
  .partial_rec([1, 3, 4]) {
    require_recoil_mass 0.495, 0.601              # eta mass window (+-3 sigma, sigma = 17.8 MeV)
  }

# -------- Chain B: J/psi -> phi eta', eta' invisible --------
sel_etap_inv = base_invisible.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  }
  .partial_rec([1, 3, 4]) {
    require_recoil_mass 0.80, 1.20                # eta' recoil-mass window
  }

# -------- Chain C: J/psi -> phi eta, eta -> gamma gamma --------
sel_eta_gg = base_gg.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  }
  # 4C kinematic fit to J/psi -> K+ K- gamma gamma (loose chi2 cut kept in BOSS;
  # optimal chi2 < 90 for eta applied later in ROOT)
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.99, 1.10)       # phi mass window
    invariant_mass_of(:gamma, :gamma).within(0.35, 0.75) # eta mass window
    chi2_cut 200
  }

# -------- Chain D: J/psi -> phi eta', eta' -> gamma gamma --------
sel_etap_gg = base_gg.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  }
  # 4C kinematic fit to J/psi -> K+ K- gamma gamma (loose chi2 cut kept in BOSS;
  # optimal chi2 < 40 for eta' applied later in ROOT)
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.99, 1.10)       # phi mass window
    invariant_mass_of(:gamma, :gamma).within(0.75, 1.15) # eta' mass window
    chi2_cut 200
  }

# -------- Attach decay cards and render the BOSS algorithms --------
alg_eta_inv.with_decay_card(decay_card_eta_inv).apply(sel_eta_inv)
alg_etap_inv.with_decay_card(decay_card_etap_inv).apply(sel_etap_inv)
alg_eta_gg.with_decay_card(decay_card_eta_gg).apply(sel_eta_gg)
alg_etap_gg.with_decay_card(decay_card_etap_gg).apply(sel_etap_gg)

# -------- Execute on real data, inclusive MC and all exclusive MC --------
datasets = [jpsi_data, jpsi_incMC,
            exMC_eta_inv, exMC_etap_inv, exMC_eta_gg, exMC_etap_gg,
            exMC_bkg_etac, exMC_bkg_phiKLKL, exMC_bkg_phif0]

root_files_eta_inv  = alg_eta_inv.execute_on(datasets)
root_files_etap_inv = alg_etap_inv.execute_on(datasets)
root_files_eta_gg   = alg_eta_gg.execute_on(datasets)
root_files_etap_gg  = alg_etap_gg.execute_on(datasets)