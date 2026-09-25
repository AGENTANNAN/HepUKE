# DSL for paper 2403.19256v2: Absolute BFs of 15 Ds+ hadronic decays
# with double-tag technique at √s = 4.128-4.226 GeV, 7.33 fb-1
# DT analysis: both Ds+ and Ds- reconstructed via hadronic tag modes
# 19 final states (many not in standard DTagAlg list; noted below)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Decay card for inclusive signal MC: e+e- → Ds*±Ds∓ with Ds → hadronic
# The transition photon/π0 from Ds*→Ds(γ,π0) is not reconstructed
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 D_s+ gamma PHSP;
    Enddecay

    Decay D_s+
    1.000 K+ K- pi+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 7 data sample groups: √s = 4.128+4.157 (merged), 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV
# BOSS 703/705 XYZ data with total 7.33 fb-1
scan_data = DatasetManager.real_data.where(cms_energy: {value: 4128..4226})
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 4128..4226})

scan_exMC = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_DsDs_dt"
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ============================================================
# TagAnalysis: double-tag Ds (both sides)
# ============================================================
alg = TagAnalysis.new("DsHadronicDTag")
alg.set_header(["DsHadronicDTagAlg/DsHadronicDTag.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })
   .note(:dt_method, "Double-tag technique at e+e- → Ds*±Ds∓: both Ds+ and Ds- reconstructed via hadronic decays. Transition γ/π0 from Ds*→Ds(γ,π0) not reconstructed. Signal events from Ds*±Ds∓ production (cross section ~20× larger than Ds+Ds- in this energy range)")
   .note(:nineteen_final_states, "15 decay modes (19 final states) of Ds+ measured: K_S0K+, K+K-π+, K_S0K+π0, K_S0K_S0π+, K+K-π+π0, K_S0K+π+π-, K_S0K-π+π+, π+π+π-, π+η(γγ), π+π0η(γγ), π+π+π-η(γγ), π+η'(π+π-η(γγ)/γρ), π+π0η'(π+π-η(γγ)/γρ), K_S0π+π0, K+π+π-. Many final states involve η/η'/ρ resonances not in standard DTagAlg mode list; only standard Ds hadronic tag modes are expressed in DSL.")
   .note(:st_yield_extraction, "ST yields from fits to M(Ds) distributions: MC-simulated shape convolved with Gaussian for signal, 2nd-order polynomial for background. Peaking backgrounds from D+→K_S0π+ (to K_S0K+) and Ds+→π+π+π-η_γγ (to π+η'_ππη) included. Efficiency-correction matrix C^ST accounts for crossfeed between signal modes")
   .note(:dt_yield_extraction, "DT yields from counting in signal region (|m̄ - m_Ds| < 15 MeV/c², |Δm| < 30 MeV/c²) with sideband subtraction (80 < |Δm| < 140 MeV/c²). Uniform background scale factor f derived from inclusive MC. Crossfeed matrix C^DT included in maximum-likelihood BF fit")
   .note(:ml_fit, "Maximum-likelihood fit to 266 ST yields + 2527 DT yields across 19 final states × 7 energy groups, extracting 15 BFs + 7 N(Ds+Ds-) values. BF constrained to be identical across different η/η' sub-modes per physics channel")
   .note(:recoil_mass, "ST selection uses recoil mass M_rec = √((√s - E*_Ds)²/c⁴ - |p*_Ds|²/c²). Tighter M_rec cut [2.10, 2.13] GeV/c² for π+π+π-η, π+π0η'_γρ, K_S0π+π0 modes; looser intervals per energy point for others. Best candidate: minimum |M_rec - m(Ds*+)|")
   .note(:reconstruction_details, "K_S0: π+π- pairs, |M(π+π-)-m_Ks| < 12 MeV/c², vertex χ² < 100, decay length > 2σ for K_S0K_S0π+ and K_S0π+π0. γ: E>25(50) MeV barrel(endcap), angle to track > 10°, 0 < t < 700 ns. π0/η_γγ: 1-C mass-constrained fit with χ² < 30, M(π0) in [115,150], M(η) in [490,580] MeV/c², ≥1 γ in barrel. η_3π: M(π+π-π0) in [530,560] MeV/c². η'_ππη: M(π+π-η) in [943,973] MeV/c². η'_γρ: M(γπ+π-) in [946,970] MeV/c². ρ0: M(π+π-) in [570,970] MeV/c². K+π+π- mode: M(π+π-) outside [487,511] MeV/c² (K_S0 veto)")
   .note(:cp_asymmetry, "CP asymmetry A_CP,i = (N_i/ε_i - N_ī/ε_ī)/(N_i/ε_i + N_ī/ε_ī) measured per mode from ST yields; all found compatible with zero")
   .with_decay_card(decay_card)

alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPi0, :DstoKsKsPi, :DstoKKPiPi0, :DstoKsKPiPi
  t.charm 1
end

alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPi0, :DstoKsKsPi, :DstoKKPiPi0, :DstoKsKPiPi
  t.charm -1
  t.rank_by :inv
end

alg.signal_side do |s|
  s.photons 0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on(scan_data + scan_incMC + scan_exMC)