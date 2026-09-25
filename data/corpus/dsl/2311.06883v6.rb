# DSL for 2311.06883v6: Evidence for Λc+ → pπ0 and Λc+ → pη
# TagAnalysis — single-tag (ST) with signal side reconstruction
# Double-tag method: tag anti-Λc- then reconstruct Λc+ → pπ0/pη on signal side
# 10 energy points: 4.600–4.843 GeV, 6.0 fb-1

# ============================================================
# Datasets
# ============================================================
data_4600 = DatasetManager.load_real_data.find("703_4600")
data_4610 = DatasetManager.load_real_data.find("706_4610")
data_4620 = DatasetManager.load_real_data.find("706_4620")
data_4640 = DatasetManager.load_real_data.find("706_4640")
data_4660 = DatasetManager.load_real_data.find("706_4660")
data_4680 = DatasetManager.load_real_data.find("706_4680")
data_4700 = DatasetManager.load_real_data.find("706_4700")
data_4740 = DatasetManager.load_real_data.find("707_4740")
data_4750 = DatasetManager.load_real_data.find("707_4750")
data_4840 = DatasetManager.load_real_data.find("707_4840")

inc_4600 = DatasetManager.load_inclusive_mc.find("703_4600")
inc_4610 = DatasetManager.load_inclusive_mc.find("706_4610")
inc_4620 = DatasetManager.load_inclusive_mc.find("706_4620")
inc_4640 = DatasetManager.load_inclusive_mc.find("706_4640")
inc_4660 = DatasetManager.load_inclusive_mc.find("706_4660")
inc_4680 = DatasetManager.load_inclusive_mc.find("706_4680")
inc_4700 = DatasetManager.load_inclusive_mc.find("706_4700")
inc_4740 = DatasetManager.load_inclusive_mc.find("707_4740")
inc_4750 = DatasetManager.load_inclusive_mc.find("707_4750")
inc_4840 = DatasetManager.load_inclusive_mc.find("707_4840")

# Signal MC: e+e- → Λc+Λc-, Λc- → tag modes, Λc+ → pπ0 or pη
sig_mc_ppi0 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card = <<~DECAY
    Decay Lambda_c+
      1.0 p+ pi0 PHSP;
    Enddecay
    Decay Lambda_c-
      1.0 p- K+ pi- PHSP;
    Enddecay
  DECAY
  c.related_dataset [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700, data_4740, data_4750, data_4840]
end

sig_mc_peta = DatasetManager.create_exclusive_mc do |c|
  c.decay_card = <<~DECAY
    Decay Lambda_c+
      1.0 p+ eta PHSP;
    Enddecay
    Decay Lambda_c-
      1.0 p- K+ pi- PHSP;
    Enddecay
  DECAY
  c.related_dataset [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700, data_4740, data_4750, data_4840]
end

# ============================================================
# TagAnalysis
# ============================================================
alg = TagAnalysis.new("Lc2pPi0pEta")

alg.set_header(["Lc2pPi0pEta/Lc2pPi0pEta.h"])
   .set_constant({ "ECMS" => [:double, 4.600] })

# Tag side: anti-Λc- reconstructed via 9 hadronic modes (Table I)
# All 9 modes verified against authoritative tag_modes.yaml
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,        # pK+π-
          :LambdacPtoKsP,          # pK0S
          :LambdacPtoKPiPi0P,      # pK+π-π0
          :LambdacPtoKsPi0P,       # pK0Sπ0
          :LambdacPtoKsPiPiP,      # pK0Sπ+π-
          :LambdacPtoLambdaPi,     # Λπ-
          :LambdacPtoLambdaPiPi0,  # Λπ-π0
          :LambdacPtoLambdaPiPiPi, # Λπ-π+π-
          :LambdacPtoPiSIGMA0LambdaGam  # Σπ (Σ0→Λγ)
  t.charm -1
end

# Signal side: proton + 2 photons (π0/η → γγ)
# π0 and η distinguished by M_γγ windows in ROOT stage
alg.signal_side do |s|
  s.charged prp: 1
  s.photons 2
  s.min_photon_angle 30.0       # >30° between photon and proton
  s.min_photon_energy 0.025     # E>25 MeV barrel or >50 MeV endcap
end

# Kinematic fit: 4-momentum conservation
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:photon_quality,
  "Photons additionally require lateral moment in (0.05, 0.40) and E3x3/E5x5 > 0.85. " \
  "Not expressible in DSL v1 — implemented via manual C++ in the algorithm.")
alg.note(:deltaE_signal_side,
  "Signal side ΔE requirement: -0.080 < ΔE_{p2γ} < 0.035 GeV. " \
  "Applied as loose pre-selection; optimal window determined in ROOT.")
alg.note(:Lambda_veto,
  "Veto events where M(p_signal, π-_from_tag) ∈ [1.111, 1.121] GeV/c² to suppress Λ background.")
alg.note(:omega_veto,
  "Veto events where M(π0_signal, π+π-_from_tag) ∈ [0.733, 0.833] GeV/c² to suppress ω background.")
alg.note(:M_gammagamma_windows,
  "π0 signal region: M_γγ ∈ [0.115, 0.150] GeV/c². " \
  "η signal region: M_γγ ∈ [0.490, 0.583] GeV/c². " \
  "Applied in ROOT stage; 2D simultaneous fit to M_BC^ST vs M_BC^{p2γ}.")
alg.note(:multiple_candidates,
  "When multiple candidates exist (~10% of events), the one with minimum |ΔE_{p2γ}| is kept.")
alg.note(:best_tag_selection,
  "ST candidate with minimum |ΔE| among all tag candidates is selected per event.")

alg.apply
alg.execute_on([data_4600, data_4610, data_4620, data_4640, data_4660,
                data_4680, data_4700, data_4740, data_4750, data_4840,
                inc_4600, inc_4610, inc_4620, inc_4640, inc_4660,
                inc_4680, inc_4700, inc_4740, inc_4750, inc_4840,
                sig_mc_ppi0, sig_mc_peta])