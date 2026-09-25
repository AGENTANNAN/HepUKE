# DSL: 2503.11015v1 — Search for X(1-+) via e+e- → γ Ds+ Ds1-(2536)
# BESIII: 5.8 fb-1 at 12 energy points √s = 4.612-4.951 GeV
# Ordinary analysis — Partial reconstruction (D_s1- → D̄*0 K- via recoil)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ── Datasets: 12 energy points ────────────────────────────────────────
# BOSS 706-707 datasets covering 4.612-4.951 GeV
ds_4610 = DatasetManager.real_data.find("706_4610")  # 4611.86 MeV
ds_4620 = DatasetManager.real_data.find("706_4620")  # 4628.00 MeV
ds_4640 = DatasetManager.real_data.find("706_4640")  # 4640.91 MeV
ds_4660 = DatasetManager.real_data.find("706_4660")  # 4661.24 MeV
ds_4680 = DatasetManager.real_data.find("706_4680")  # 4681.92 MeV
ds_4700 = DatasetManager.real_data.find("706_4700")  # 4698.82 MeV
ds_4740 = DatasetManager.real_data.find("707_4740")  # 4739.70 MeV
ds_4750 = DatasetManager.real_data.find("707_4750")  # 4750.05 MeV
ds_4780 = DatasetManager.real_data.find("707_4780")  # 4780.54 MeV
ds_4840 = DatasetManager.real_data.find("707_4840")  # 4843.07 MeV
ds_4914 = DatasetManager.real_data.find("707_4914")  # 4918.02 MeV
ds_4946 = DatasetManager.real_data.find("707_4946")  # 4950.93 MeV

all_data = [ds_4610, ds_4620, ds_4640, ds_4660, ds_4680, ds_4700,
            ds_4740, ds_4750, ds_4780, ds_4840, ds_4914, ds_4946]

# ═══════════════════════════════════════════════════════════════════════
# Mode 1: Ds+ → K+ K- π+ (with φ/K* intermediate cuts)
# ═══════════════════════════════════════════════════════════════════════
decay_card_m1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma D_s+ D_s1(2536)- PHSP;
  Enddecay

  Decay D_s+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s1(2536)-
  1.000 anti-D*0 K- PHSP;
  Enddecay

  Decay anti-D*0
  1.000 anti-D0 pi0 PHSP;
  Enddecay
DECAYCARD

exMC_m1 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_X_gamma_Ds_Ds1_KKpi"
  config.events        = 500_000
  config.decay_card    = decay_card_m1
  config.cross_section = :default
end

# ── Algorithm ─────────────────────────────────────────────────────────
alg_m1 = Algorithm.new("X1m_Ds_Ds1_KKpi")
alg_m1.set_header(["X1m_Ds_Ds1_KKpiAlg/X1m_Ds_Ds1_KKpi.h"])
    .set_constant({ "ECMS" => [:double, 4.612] })   # will be overridden per-execute_on dataset
    .with_decay_card(decay_card_m1)

event_selection_m1 = Selection.new
  # Track selection
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  end
  # Photon selection: at least 2 photons (Ds1 transition + D*0 decay)
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # PID: identify kaons and pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==1"
  end
  # Remove identified kaons, assign remainder as pions
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  end
  # Partial reconstruction: reconstruct Ds+ (KKπ) + K- from Ds1 decay;
  # infer Ds1-(2536) and D*0 from recoil mass
  # Ds+ mass window constraint and RM(γ Ds+ K-) ~ m_D*0
  .partial_rec([1, 2, 3, 4, 5, 6]) do
    # recIDs: 1=gamma, 2=Ds+, 3=Ds1-, 4=anti-D*0, 5=pi0, 6=K- (from Ds1)
    best_combination_by_mass :D_s, 1.968    # Ds+ nominal mass (GeV)
    require_recoil_mass 2.000, 2.020         # D*0 recoil mass window ~2.007 GeV
  end

alg_m1.with_decay_card(decay_card_m1).apply(event_selection_m1)
alg_m1.execute_on(all_data + [exMC_m1].flatten)

alg_m1.note(:multi_energy_fit, "Three energy sets (I: 4.612-4.699, II: 4.740-4.805, III: 4.843-4.951) for optimized cuts per Table II")
alg_m1.note(:ds_reconstruction, "Ds+ → K+K-π+ with φ(→K+K-) and K*(→K±π∓) intermediate mass windows and helicity angle cuts")
alg_m1.note(:partial_reco_detail, "D_s1- → D̄*0 K-: D̄*0 not fully reconstructed; RM(γ Ds+ K-) = m_D*0 constraint applied; Ds+ mass constrained to 1.968 GeV")
alg_m1.note(:background_veto, "Veto on Ds+ → K_S0K+ (π0 veto) and Ds*±Ds∓ backgrounds")
alg_m1.note(:x_signal_search, "X(1-+) search via γ Ds+ Ds1- cross section enhancement; ROOT-stage fit to γ Ds+ mass spectrum")
alg_m1.note(:helicity_angles, "φ → K+K- helicity angle |cosθ_hel| > 0.5; K*(892) mass window ±50 MeV; helicity angle cut applied in ROOT")

# ═══════════════════════════════════════════════════════════════════════
# Mode 2: Ds+ → K_S0 K+ (K_S0 → π+ π-)
# ═══════════════════════════════════════════════════════════════════════
decay_card_m2 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma D_s+ D_s1(2536)- PHSP;
  Enddecay

  Decay D_s+
  1.000 K_S0 K+ PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay D_s1(2536)-
  1.000 anti-D*0 K- PHSP;
  Enddecay

  Decay anti-D*0
  1.000 anti-D0 pi0 PHSP;
  Enddecay
DECAYCARD

exMC_m2 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_X_gamma_Ds_Ds1_KsK"
  config.events        = 500_000
  config.decay_card    = decay_card_m2
  config.cross_section = :default
end

# ── Algorithm ─────────────────────────────────────────────────────────
alg_m2 = Algorithm.new("X1m_Ds_Ds1_KsK")
alg_m2.set_header(["X1m_Ds_Ds1_KsKAlg/X1m_Ds_Ds1_KsK.h"])
    .set_constant({ "ECMS" => [:double, 4.612] })
    .with_decay_card(decay_card_m2)

event_selection_m2 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # PID: identify kaons
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"
  end
  .remove([:kp <= :chrgp])
  # K_S0 → π+π- secondary vertex fit
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Remaining track is K- from Ds1 decay
  .assign({ chrgn: :km })
  # Partial reconstruction: reconstruct Ds+(K_S0 K+) + K-(from Ds1)
  # and transition photon; infer Ds1- via recoil mass
  .partial_rec([1, 2, 3, 4, 5, 6]) do
    best_combination_by_mass :D_s, 1.968
    require_recoil_mass 2.000, 2.020
  end

alg_m2.with_decay_card(decay_card_m2).apply(event_selection_m2)
alg_m2.execute_on(all_data + [exMC_m2].flatten)

alg_m2.note(:ks_reconstruction, "K_S0 → π+π- with secondary vertex fit; invariant mass in (487,511) MeV/c² (±3σ); decay length > 2× vertex resolution")
alg_m2.note(:ds_mass_window, "Ds+ mass from KKπ/K_S0K combinations constrained to 1.968 GeV; φ mass window ±10 MeV/c²; K* mass window ±50 MeV/c²")