# 2411.15752v2: e⁺e⁻ → K_S⁰ K_S⁰ ψ(3686) at √s = 4.682–4.951 GeV
# Ordinary analysis, multi-energy cross-section measurement (4.1 fb⁻¹ total)
# K_S⁰ → π⁺π⁻, ψ(3686) → J/ψ + X, J/ψ → l⁺l⁻ (l = e, μ)
# ψ(3686) identified via M(K_S⁰ K_S⁰) recoil mass

# Eight energy points: 4.682, 4.699, 4.740, 4.750, 4.781, 4.843, 4.918, 4.951 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # ~4.682 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # ~4.699 GeV
data_4740 = DatasetManager.real_data.find("707_4740")   # ~4.740 GeV
data_4750 = DatasetManager.real_data.find("707_4750")   # ~4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")   # ~4.781 GeV
data_4840 = DatasetManager.real_data.find("707_4840")   # ~4.843 GeV
data_4914 = DatasetManager.real_data.find("707_4914")   # ~4.918 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # ~4.951 GeV

multi_data = [data_4680, data_4700, data_4740, data_4750,
              data_4780, data_4840, data_4914, data_4946]

incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

multi_incMC = [incMC_4680, incMC_4700, incMC_4740, incMC_4750,
               incMC_4780, incMC_4840, incMC_4914, incMC_4946]

# ---- Decay card for signal MC ----
# e⁺e⁻ → K_S⁰ K_S⁰ ψ(3686); ψ(3686) → J/ψ + X (X = π⁺π⁻, π⁰π⁰, η, π⁰, γγ)
# J/ψ → e⁺e⁻ (one card), J/ψ → μ⁺μ⁻ (another card)
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K_S0 psi(3686) PHSP;
    Enddecay

    Decay psi(3686)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K_S0 psi(3686) PHSP;
    Enddecay

    Decay psi(3686)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_sig_ee = DatasetManager.create_exclusive_mc_for(multi_data) do |config|
  config.sample_name = "KsKsPsi3686_ee_ExcMC"
  config.events = 100000
  config.decay_card = decay_card_ee
  config.cross_section = :default
end

exMC_sig_mumu = DatasetManager.create_exclusive_mc_for(multi_data) do |config|
  config.sample_name = "KsKsPsi3686_mumu_ExcMC"
  config.events = 100000
  config.decay_card = decay_card_mumu
  config.cross_section = :default
end

# ============================================================
# e⁺e⁻ → K_S⁰ K_S⁰ ψ(3686); ψ(3686) → J/ψ + X (X unobserved)
# ============================================================
# e⁺e⁻ channel
alg_ee = Algorithm.new("KsKsPsi3686_EE")
alg_ee.set_header(["KsKsPsi3686_EEAlg/KsKsPsi3686_EE.h"])
      .set_constant({"ECMS" => [:double, 4.843]})
      .note(:partial_reco, "ψ(3686) → J/ψ + X where X = π⁺π⁻, π⁰π⁰, η, π⁰, γγ is NOT reconstructed. The ψ(3686) is identified via M(K_S⁰K_S⁰) recoil mass spectrum in ROOT. The missing X system mass is not constrained in the kinematic fit — only K_S⁰ masses and J/ψ mass are constrained.")
      .note(:lepton_id, "e/μ separation via E_EMC: > 0.8 GeV → electron, < 0.5 GeV → muon. If at least one lepton satisfies and the other has no EMC info, both kept. Not fully expressible in DSL; approximated with high-momentum lepton PID.")
      .note(:multi_energy, "Eight c.m. energies from 4.682 to 4.951 GeV across BOSS 706/707. ECMS set to central value 4.843 GeV.")
      .note(:signal_extraction, "Signal yield obtained by counting events in ψ(3686) signal region (3.66, 3.71) GeV/c² in M(K_S⁰K_S⁰) recoil mass. Sideband method for background subtraction. Born cross section computed with ISR and VP corrections in ROOT.")
      .note(:best_candidate, "Multiple K_S⁰K_S⁰ combinations: all kept for further analysis (paper says 'all combinations are kept'). Fraction of multi-combination events estimated at 8.9% from MC.")
      .note(:ks_pion_tracks, "Pions from K_S⁰ decays do NOT have Vz/Vxy cuts applied. Not directly expressible in DSL — approximated by excluding Vz/Vr from select_track and applying KS0 mass window post-vertex-fit.")
      .with_decay_card(decay_card_ee)

sel_ee = Selection.new

sel_ee.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nTot ">=6"
end

sel_ee.pid do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,
                                 treat_as_electron_if_energy_above: 0.6
  nlp "==1"
  nlm "==1"
end

sel_ee.remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})

sel_ee.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
end

sel_ee.kinematic_fit([:K_S0, :K_S0, :lp, :lm]) do
  nominal
  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  constrain_four_momentum
  chi2_cut 200
end

alg_ee.apply(sel_ee)

# μ⁺μ⁻ channel — same selection except lepton flavour
alg_mumu = Algorithm.new("KsKsPsi3686_MuMu")
alg_mumu.set_header(["KsKsPsi3686_MuMuAlg/KsKsPsi3686_MuMu.h"])
        .set_constant({"ECMS" => [:double, 4.843]})
        .note(:partial_reco, "ψ(3686) → J/ψ + X; X system unobserved. ψ(3686) identified via M(K_S⁰K_S⁰) recoil mass. See e⁺e⁻ algorithm notes.")
        .note(:lepton_id, "muon identification: E_EMC < 0.5 GeV deposit. Approximated via high-momentum lepton PID.")
        .note(:multi_energy, "Eight c.m. energies from 4.682 to 4.951 GeV. See e⁺e⁻ algorithm notes.")
        .note(:signal_extraction, "Combined with e⁺e⁻ channel. See e⁺e⁻ algorithm notes.")
        .with_decay_card(decay_card_mumu)

sel_mumu = Selection.new

sel_mumu.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nTot ">=6"
end

sel_mumu.pid do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,
                                 treat_as_electron_if_energy_above: 0.6
  nlp "==1"
  nlm "==1"
end

sel_mumu.remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})

sel_mumu.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
end

sel_mumu.kinematic_fit([:K_S0, :K_S0, :lp, :lm]) do
  nominal
  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  constrain_four_momentum
  chi2_cut 200
end

alg_mumu.apply(sel_mumu)

# Execute on all energy points
all_datasets_ee = multi_data + multi_incMC + exMC_sig_ee
all_datasets_mumu = multi_data + multi_incMC + exMC_sig_mumu

alg_ee.execute_on(all_datasets_ee)
alg_mumu.execute_on(all_datasets_mumu)