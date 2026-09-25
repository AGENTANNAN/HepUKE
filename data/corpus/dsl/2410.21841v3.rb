# 2410.21841v3: Search for Λ-Λbar oscillation in J/ψ → Λ Λbar
# Ordinary analysis at J/ψ (√s = 3.097 GeV), (10087±44)×10⁶ J/ψ events
# Λ → pπ⁻, Λbar → pbar π⁺

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Decay card for RS signal: J/ψ → Λ Λbar ----
decay_card_rs = <<~DECAYCARD
    Decay J/psi
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_rs = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_Lambda_Lambdabar_RS_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 10000000
  config.decay_card = decay_card_rs
  config.cross_section = :default
end

# ---- Decay card for WS signal: J/ψ → Λ Λ (oscillation) ----
decay_card_ws = <<~DECAYCARD
    Decay J/psi
    1.000 Lambda0 Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

exMC_ws = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_Lambda_Lambda_WS_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 1000000
  config.decay_card = decay_card_ws
  config.cross_section = :default
end

# ---- Peaking background: J/ψ → Λ Σbar⁰ + c.c. ----
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.000 Lambda0 anti-Sigma0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Sigma0
    1.000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_Lambda_Sigmabar0_bkg_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_card_bkg
  config.cross_section = :default
end

# ============================================================
# J/ψ → Λ Λbar (Right Sign) and J/ψ → Λ Λ + c.c. (Wrong Sign)
# ============================================================
alg = Algorithm.new("JpsiToLLbar")
alg.set_header(["JpsiToLLbarAlg/JpsiToLLbar.h"])
   .set_constant({"ECMS" => [:double, 3.097]})
   .note(:low_momentum_pion, "Tracks with p < 0.5 GeV/c are taken as pion candidates without PID; tracks with p ≥ 0.5 GeV/c use probability PID for proton vs kaon/pion. Momentum-dependent particle-type assignment not expressible in DSL; approximated with probability PID on all tracks.")
   .note(:mass_window, "Λ/Λbar mass window set to ±3σ (σ = 1.8 MeV) from signal MC; applied in ROOT")
   .note(:best_candidate, "If multiple ΛΛbar pairs survive, the one with minimum χ²_4C is retained; automatic in DSL via kinematic_fit best-combination selection")
   .note(:signal_extraction, "Simultaneous fit of M(pπ⁻) and M(pbar π⁺) distributions in ROOT; signal shape from MC convolved with Gaussian for resolution difference")
   .with_decay_card(decay_card_rs)

sel = Selection.new

sel.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot ">=4"
end

sel.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
  nprm ">=1"
end

sel.remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})

sel.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

sel.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

sel.kinematic_fit([:Lambda, :Lambda_bar]) do
  nominal
  constrain_four_momentum
  chi2_cut 50
end

alg.apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_rs, exMC_ws, exMC_bkg])