# ============================================================
# BESIII Double-Tag Cross Section Measurement
# Paper: arXiv:1803.06293
# Measurement of e+e- -> DDbar Cross Sections at psi(3770)
# sqrt(s) = 3.773 GeV, L = 2.93 fb-1
# ============================================================

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ========================================================================
# D0 D0bar analysis — 3 x 3 tag mode combinations (9 DT pairs)
# D0 tag modes: K-pi+, K-pi+pi0, K-pi+pi+pi-
# Anti-D0 tag: charm -1
# ========================================================================

decay_card_D0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_D0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0D0bar_XS_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_D0
  config.cross_section = :default
end

alg_D0 = TagAnalysis.new("D0D0bar_XS")
alg_D0.set_header(["D0D0bar_XSAlg/D0D0bar_XS.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_D0)

# Tag side 1: D0 (charm +1) with 3 hadronic modes
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

# Tag side 2: anti-D0 (charm -1) with same 3 hadronic modes
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_D0.signal_side do |s|
  s.photons 0
end

alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_D0
  .note(:xs_measurement, "Cross section measurement: N_D0D0bar from DT/ST ratio with efficiency cancellation; quantum correlation corrections from Asner & Sun (2006)")
  .note(:tag_deltaE, "Mode-dependent asymmetric DeltaE: Kpi +-3sigma about mean; modes with pi0 extended to -4sigma low side. Best tag by min |DeltaE|")
  .note(:tag_mbc, "MBC signal region: D0 modes 1.858-1.874 GeV/c2. MBC fit: MC-convolved double-Gaussian signal + ARGUS bkg; tag yields from signal region minus ARGUS integral")
  .note(:dt_selection, "DT events: opposite net charge, opposite charm, no shared tracks; best DT by (MBC_D + MBC_Dbar)/2 closest to D mass")
  .note(:ks0_reco, "K_S0->pi+pi-: secondary vertex fit chi2<100, mass 487-511 MeV; KS0 pion daughters exempt from IP/PID requirements")
  .note(:pi0_reco, "pi0->gamma gamma: invariant mass 115-150 MeV; at least one photon in barrel; refit with pi0 mass constraint")
  .note(:photon_selection, "E>25 MeV barrel(|cos_theta|<0.8), E>50 MeV endcap(0.84<|cos_theta|<0.92); EMC TDC [0,700]ns; no track association")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID: TOF+dE/dx, highest probability hypothesis")
  .note(:dt_mbc_fit, "2D MBC fit: MC signal shape + 4-component bkg model (DTag-ARGUS product, mispartitioned continuum, ARGUSxARGUS)")
  .note(:cosmic_veto, "D0->Kpi cosmic/QED veto: TOF cosmic, e+e- PID, EMC e+e-, muon PID+MUC rejection")
  .note(:quantum_correlation, "D0D0bar C=-1 correlation correction: -0.2% overall effect; strong phase factors from HFLAV and CLEO inputs")
  .note(:peaking_bkg, "MC-determined peaking background subtraction (DCS modes, K/pi mis-ID, fake KS0); max 2.5% for D+->Kspi3pi")
  .apply

alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0])


# ========================================================================
# D+ D- analysis — 6 x 6 tag mode combinations (36 DT pairs)
# D+ tag modes: K-pi+pi+, K-pi+pi+pi0, Kspi+, Kspi+pi0, Kspi+pi+pi-, K+K-pi+
# ========================================================================

decay_card_Dp = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay
    Decay D+
    1.000 K- pi+ pi+ PHSP;
    Enddecay
    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_Dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "DpDm_XS_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card_Dp
  config.cross_section = :default
end

alg_Dp = TagAnalysis.new("DpDm_XS")
alg_Dp.set_header(["DpDm_XSAlg/DpDm_XS.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_Dp)

# Tag side 1: D+ (charm +1) with 6 hadronic modes
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm 1
end

# Tag side 2: D- (charm -1) with same 6 hadronic modes
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.rank_by :inv
end

alg_Dp.signal_side do |s|
  s.photons 0
end

alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:"D+")
  f.chi2_cut 200
end

alg_Dp
  .note(:xs_measurement, "Cross section measurement: N_DpDm from DT/ST ratio; yield-weighted mean N_DpDm = (8296+-31)x10^3, sigma(e+e-->D+D-) = (2.830+-0.011) nb (stat)")
  .note(:tag_deltaE, "Mode-dependent asymmetric DeltaE: +-3sigma about mean; modes with pi0 extended to -4sigma low side. Best tag by min |DeltaE|")
  .note(:tag_mbc, "MBC signal region: D+ modes 1.8628-1.8788 GeV/c2. MBC fit: MC-convolved double-Gaussian signal + ARGUS bkg")
  .note(:dt_selection, "DT events: opposite net charge, opposite charm, no shared tracks; best DT by (MBC_D + MBC_Dbar)/2 closest to D mass")
  .note(:ks0_reco, "K_S0->pi+pi- reconstructed with secondary vertex fit chi2<100, mass 487-511 MeV; KS0 daughters exempt from IP/PID requirements")
  .note(:pi0_reco, "pi0->gamma gamma: invariant mass 115-150 MeV; at least one photon in barrel; refit with pi0 mass constraint")
  .note(:photon_selection, "E>25 MeV barrel(|cos_theta|<0.8), E>50 MeV endcap(0.84<|cos_theta|<0.92); EMC TDC [0,700]ns; no track association")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID: TOF+dE/dx, highest probability hypothesis")
  .note(:dt_mbc_fit, "2D MBC fit: MC signal shape + 4-component bkg model; 36 DT pairs combined with weighted mean")
  .note(:peaking_bkg, "MC-determined peaking background subtraction; max ~2.5% for D+->Kspi3pi")
  .apply

alg_Dp.execute_on([psi3770_data, psi3770_incMC, exMC_Dp])