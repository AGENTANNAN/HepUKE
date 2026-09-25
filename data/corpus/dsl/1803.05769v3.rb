# ============================================================
# BESIII Double-Tag Analysis of D0 -> All-Neutral Decays
# Paper: arXiv:1803.05769
# Measurement of Singly Cabibbo-Suppressed Decays
#   D0 -> pi0pi0pi0, pi0pi0eta, pi0etaeta, etaetaeta
# sqrt(s) = 3.773 GeV, L = 2.93 fb-1 at BESIII
# ============================================================

# --- Datasets ---
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ========================================================================
# Mode I: D0 -> pi0 pi0 pi0
# Signal: all-neutral, 3 pi0 from 6 photons
# Tag: anti-D0 -> K+pi-, K+pi-pi0, K+pi-pi-pi+
# ========================================================================

decay_card_3pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 pi0 pi0 pi0 PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0to3pi0_DT_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_3pi0
  config.cross_section = :default
end

alg_3pi0 = TagAnalysis.new("D0to3pi0_DT")
alg_3pi0.set_header(["D0to3pi0_DTAlg/D0to3pi0_DT.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_3pi0)

alg_3pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.83
end

alg_3pi0.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_3pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_3pi0
  .note(:tag_photon_limit, "DSL v1 TagAnalysis consumes 0..2 signal photons; analysis requires >=6 for 3pi0 (3x pi0->gammagamma)")
  .note(:pi0_constraints, "Three 1C pi0 mass constraints needed; only 1 expressible due to v1 photon limit")
  .note(:tag_deltaE_windows, "Mode-dependent ST DeltaE: Kpi (-0.027,0.025), KpiPi0 (-0.071,0.041), Kpi3pi (-0.025,0.022) GeV; best tag by min |DeltaE|")
  .note(:photon_selection, "E>25 MeV barrel(|cos_theta|<0.8), E>50 MeV endcap(0.86<|cos_theta|<0.92); EMC TDC [0,700]ns")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID from TOF+dE/dx, highest C.L. assigned")
  .note(:signal_deltaE, "D0 signal DeltaE window: (-0.115, 0.059) GeV; best signal candidate by min |DeltaE|")
  .note(:background_veto, "D0->4pi0 veto: chi2_4pi>20 when >=4 independent pi0; KS0->pi0pi0 veto: M(pi0pi0) in [445,535] MeV; cross-feed: chi2_pi0pi0eta>20; no other combo chi2<20")
  .note(:mbc_signal_fit, "M_BC unbinned ML fit: MC-convolved Gaussian signal + ARGUS comb. bkg + peaking bkg (D0->4pi0, KS0pi0, cross-feeds)")
  .note(:mc_model, "Efficiency from 3-body phase space MC (no intermediate resonances observed in Dalitz plot). ST efficiency from inclusive MC.")
  .apply

alg_3pi0.execute_on([psi3770_data, psi3770_incMC, exMC_3pi0])


# ========================================================================
# Mode II: D0 -> pi0 pi0 eta
# Signal: all-neutral, 2pi0 + 1eta from 6 photons
# ========================================================================

decay_card_2pi0eta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 pi0 pi0 eta PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_2pi0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0to2pi0eta_DT_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_2pi0eta
  config.cross_section = :default
end

alg_2pi0eta = TagAnalysis.new("D0to2pi0eta_DT")
alg_2pi0eta.set_header(["D0to2pi0eta_DTAlg/D0to2pi0eta_DT.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_2pi0eta)

alg_2pi0eta.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.83
end

alg_2pi0eta.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_2pi0eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_2pi0eta
  .note(:tag_photon_limit, "DSL v1 TagAnalysis consumes 0..2 signal photons; analysis requires >=6 for 2pi0+eta (4g+2g)")
  .note(:resonance_constraints, "2pi0 + 1eta 1C mass constraints needed; only 1 pi0 constraint expressible due to v1 photon limit")
  .note(:tag_deltaE_windows, "Mode-dependent ST DeltaE: Kpi (-0.027,0.025), KpiPi0 (-0.071,0.041), Kpi3pi (-0.025,0.022) GeV")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC TDC [0,700]ns")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID max C.L.")
  .note(:signal_deltaE, "D0 signal DeltaE window: (-0.088, 0.053) GeV; best candidate by min |DeltaE|")
  .note(:background_veto, "D0->4pi0 veto: chi2_4pi>20; KS0->pi0pi0 veto: M(pi0pi0) in [445,535] MeV; cross-feed: chi2_3pi0>20; no other combo chi2<20")
  .note(:mbc_signal_fit, "M_BC unbinned ML fit with peaking backgrounds (D0->4pi0, KS0eta, cross-feeds)")
  .note(:mc_model, "Efficiency from 3-body phase space MC (no intermediate resonances observed). ST efficiency from inclusive MC.")
  .apply

alg_2pi0eta.execute_on([psi3770_data, psi3770_incMC, exMC_2pi0eta])


# ========================================================================
# Mode III: D0 -> pi0 eta eta
# Dominant process: D0 -> a0(980)0 eta -> pi0 eta eta (significance 2.6sigma)
# Efficiency from D0->a0(980)0 eta MC with Flatte a0(980) parameters
# ========================================================================

decay_card_pi0etaeta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 pi0 eta eta PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_pi0etaeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0topi0etaeta_DT_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_pi0etaeta
  config.cross_section = :default
end

alg_pi0etaeta = TagAnalysis.new("D0topi0etaeta_DT")
alg_pi0etaeta.set_header(["D0topi0etaeta_DTAlg/D0topi0etaeta_DT.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_pi0etaeta)

alg_pi0etaeta.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.83
end

alg_pi0etaeta.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_pi0etaeta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_pi0etaeta
  .note(:tag_photon_limit, "DSL v1 TagAnalysis consumes 0..2 signal photons; analysis requires >=6 for pi0+2eta (2g+2g+2g)")
  .note(:resonance_constraints, "1pi0 + 2eta 1C mass constraints needed; only 1 pi0 constraint expressible due to v1 photon limit")
  .note(:intermediate_resonance, "Dominant process D0->a0(980)0 eta->pi0 eta eta; a0(980)0 significance 2.6sigma; Flatte params from Crystal Barrel")
  .note(:tag_deltaE_windows, "Mode-dependent ST DeltaE: Kpi (-0.027,0.025), KpiPi0 (-0.071,0.041), Kpi3pi (-0.025,0.022) GeV")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC TDC [0,700]ns")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID max C.L.")
  .note(:signal_deltaE, "D0 signal DeltaE window: (-0.061, 0.045) GeV; best candidate by min |DeltaE|")
  .note(:background_veto, "Cross-feed veto: chi2_3pi0>20 AND chi2_2pi0eta>20; no other combo chi2<20; D0->4pi0 veto not applied (no strong contamination)")
  .note(:mc_model, "Efficiency from D0->a0(980)0 eta->pi0 eta eta MC (a0(980)0 observed as dominant intermediate state). ST efficiency from inclusive MC.")
  .note(:mbc_signal_fit, "M_BC unbinned ML fit with peaking backgrounds. Dalitz plot fit for M(pi0 eta): a0(980)0 signal yield 21+-5 events.")
  .apply

alg_pi0etaeta.execute_on([psi3770_data, psi3770_incMC, exMC_pi0etaeta])


# ========================================================================
# Mode IV: D0 -> eta eta eta
# No significant signal observed; upper limit set at 90% C.L.
# ========================================================================

decay_card_3eta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 eta eta eta PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_3eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0to3eta_DT_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_3eta
  config.cross_section = :default
end

alg_3eta = TagAnalysis.new("D0to3eta_DT")
alg_3eta.set_header(["D0to3eta_DTAlg/D0to3eta_DT.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_3eta)

alg_3eta.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.83
end

alg_3eta.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_3eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_3eta
  .note(:tag_photon_limit, "DSL v1 TagAnalysis consumes 0..2 signal photons; analysis requires >=6 for 3eta (3x eta->gammagamma)")
  .note(:resonance_constraints, "3 eta 1C mass constraints needed; only 1 eta constraint expressible due to v1 photon limit")
  .note(:no_signal, "No significant D0->etaetaeta signal observed; upper limit B<1.3e-4 at 90% C.L.")
  .note(:tag_deltaE_windows, "Mode-dependent ST DeltaE: Kpi (-0.027,0.025), KpiPi0 (-0.071,0.041), Kpi3pi (-0.025,0.022) GeV")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC TDC [0,700]ns")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID max C.L.")
  .note(:signal_deltaE, "D0 signal DeltaE window: (-0.030, 0.028) GeV; best candidate by min |DeltaE|")
  .note(:mc_model, "Efficiency from 3-body phase space MC (small phase space, no intermediate resonance study). ST efficiency from inclusive MC.")
  .note(:mbc_signal_fit, "M_BC fit: MC-convolved Gaussian signal + ARGUS background; upper limit from normalized likelihood with systematic uncertainties incorporated")
  .apply

alg_3eta.execute_on([psi3770_data, psi3770_incMC, exMC_3eta])