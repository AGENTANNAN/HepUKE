# DSL for paper 2403.19091v1: Observation of D0→K_S0π-π0e+νe and D+→K_S0π+π-e+νe
# Tag-based analysis: double-tag method at ψ(3770), 2.93 fb-1
# Two separate TagAnalysis objects for D0 and D+ tag sides

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# Common dataset setup: ψ(3770), 2.93 fb-1
# ============================================================
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for signal MC: ψ(3770) → D0 D0bar (or D+ D-)
# D → K_S0 π π e+ ν_e via ISGW2 model
decay_card_D0_sig = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 pi- pi0 e+ nu_e ISGW2;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_Dp_sig = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 pi+ pi- e+ nu_e ISGW2;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC
exMC_D0_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0toKsPiPi0Enu"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_D0_sig
  config.cross_section = :default
end

exMC_Dp_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_DptoKsPiPiEnu"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_Dp_sig
  config.cross_section = :default
end

# ============================================================
# Algorithm 1: D0 tag → signal D0→K_S0π-π0e+νe
# ============================================================
alg_D0 = TagAnalysis.new("D0toKsPiPi0EnuTag")
alg_D0.set_header(["D0toKsPiPi0EnuTagAlg/D0toKsPiPi0EnuTag.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .note(:dt_method, "Double-tag method at ψ(3770): ST anti-D0 from 4 hadronic tag modes (K+π-, K+π-π0, K+π-π+π-, K+π-π+π-π0); signal D0→K_S0π-π0e+νe reconstructed from remaining tracks")
       .note(:ks_reconstruction, "K_S0 reconstructed from π+π- pairs in signal side: vertex-constrained secondary vertex fit, flight distance > 2σ, invariant mass within (0.486, 0.510) GeV/c²")
       .note(:pi0_reconstruction, "π0 → γγ in signal side: photon energies > 25(50) MeV barrel(endcap), invariant mass within (0.115, 0.150) GeV/c², mass-constrained kinematic fit applied")
       .note(:electron_pid, "Positron PID: L'(e) > 0.001 and L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8; EMC energy/momentum ratio cut: E/|p|c > 0.18×χ²(dE/dx) + 0.32")
       .note(:fsr_recovery, "Photons within 5° cone of positron with E > 50 MeV added back to positron 4-momentum")
       .note(:signal_selection_cuts, "D0 signal: π0 energy > 0.22 GeV, |cos θ(π0)| < 0.83; M(K_S0π-π0π+[e→π]) < 1.78 GeV/c² to suppress D0→K_S0π+π-π0 bkg")
       .note(:missing_mass_fit, "DT yield from 2D unbinned maximum-likelihood fit to M²_miss vs M(K_S0ππ) using RooNDKeysPdf signal and background shapes; M²_miss = (E_beam - ΣE_i)² - |Σp_i|² after 4C kinematic fit with D mass constraints")
       .note(:k1_1270_resonance, "Signal dominated by K1(1270) axial-vector resonance: mass 1.253 GeV/c², width 90 MeV, relativistic Breit-Wigner. BFs of K1(1270)→K_S0ππ sub-decays from PDG")
       .with_decay_card(decay_card_D0_sig)

alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
  t.charm -1
end

alg_D0.signal_side do |s|
  s.photons 2
  s.charged(pip: 1, pim: 2, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

alg_D0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0.apply
alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0_sig])

# ============================================================
# Algorithm 2: D+ tag → signal D+→K_S0π+π-e+νe
# ============================================================
alg_Dp = TagAnalysis.new("DptoKsPiPiEnuTag")
alg_Dp.set_header(["DptoKsPiPiEnuTagAlg/DptoKsPiPiEnuTag.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .note(:dt_method_Dp, "Double-tag method at ψ(3770): ST D- from 6 hadronic tag modes (K+π-π-, K+π-π-π0, K_S0π-, K_S0π-π0, K_S0π+π-π-, K+K-π-); signal D+→K_S0π+π-e+νe reconstructed from remaining tracks")
       .note(:ks_reconstruction, "K_S0 reconstructed from π+π- pairs in signal side: vertex-constrained secondary vertex fit, flight distance > 2σ, invariant mass within (0.486, 0.510) GeV/c²")
       .note(:electron_pid, "Positron PID: L'(e) > 0.001 and L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8; EMC energy/momentum ratio cut: E/|p|c > 0.18×χ²(dE/dx) + 0.32")
       .note(:fsr_recovery, "Photons within 5° cone of positron with E > 50 MeV added back to positron 4-momentum")
       .note(:signal_selection_cuts_Dp, "D+ signal: M(K_S0π+π-π+[e→π]) < 1.83 GeV/c²; cos θ(α) between e+ and π- < 0.95 (Dalitz veto); M(K_S0π+π-π+[e→π]π0) < 1.4 GeV/c²; cos θ(β) between missing momentum and most energetic unused shower < 0.88")
       .note(:missing_mass_fit, "DT yield from 2D unbinned maximum-likelihood fit to M²_miss vs M(K_S0ππ) using RooNDKeysPdf signal and background shapes; M²_miss = (E_beam - ΣE_i)² - |Σp_i|² after 4C kinematic fit with D mass constraints")
       .note(:k1_1270_resonance, "Signal dominated by K1(1270) axial-vector resonance: mass 1.253 GeV/c², width 90 MeV. BFs of K1(1270)→K_S0ππ sub-decays from Belle measurement")
       .with_decay_card(decay_card_Dp_sig)

alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_Dp.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 2, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Dp.apply
alg_Dp.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_sig])