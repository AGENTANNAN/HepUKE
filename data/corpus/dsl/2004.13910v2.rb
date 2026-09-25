# Paper 2004.13910v2: D meson hadronic decays to eta at psi(3770) — DT method (TAG-BASED)
# e+e- -> psi(3770) -> D Dbar
# ST tag: D- (via 6 modes) or anti-D0 (via 3 modes)
# Signal D decays from remaining tracks: 14 modes with eta -> gamma gamma
# Two TagAnalysis objects needed (tag species differ: D+ vs D0)

decay_card = <<~DECAYCARD
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay psi(3770)
1.0 D+ D- PHSP;
Enddecay
End
DECAYCARD

psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ===== TagAnalysis 1: ST D- (tag D+), signal D0 -> eta + hadrons =====
alg_Dplus_tag = TagAnalysis.new("D0_to_eta_DT_DplusTag")
alg_Dplus_tag.set_header(["D0_to_eta_DT_DplusTagAlg/D0_to_eta_DT_DplusTag.h"])
alg_Dplus_tag.set_constant({ "ECMS" => [:double, 3.773] })
alg_Dplus_tag.with_decay_card(decay_card)

# ST tag D+ (D- decays): 6 modes
alg_Dplus_tag.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm 1
end

# Signal anti-D0 from remaining tracks: needs eta -> 2 photons + charged tracks
alg_Dplus_tag.signal_side do |s|
  s.photons 2
  s.charged(at_least: true)
  s.min_photon_angle 10.0
end

alg_Dplus_tag.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_Dplus_tag.dtag_reconstruction do |d|
  d.beam_energy [:constant, 1.8865]
end

# 6 D0 signal modes available in DTagAlg:
# D0toKPiEta, D0toKsPi0Eta, D0toKKEta, D0toKPiPi0Eta, D0toKsPiPiEta, D0toPiPiPi0Eta
# 2 D0 modes NOT in DTagAlg (D0toKsKsEta, D0toKsPi0Pi0Eta) — dropped
alg_Dplus_tag.note(:tag_mode_unavailable,
  "2 of 8 D0 signal modes not in DTagAlg: D0->K_S0 K_S0 eta, D0->K_S0 pi0 pi0 eta. " \
  "Available: D0toKPiEta, D0toKsPi0Eta, D0toKKEta, D0toKPiPi0Eta, D0toKsPiPiEta, D0toPiPiPi0Eta."
)

alg_Dplus_tag.note(:selection_details,
  "eta -> gamma gamma: |M(gamma gamma)-m(eta)| < 27.5 MeV, 1C Kalman fit. " \
  "K_S0 -> pi+ pi-: L/sigma_L > 2, mass window (0.486, 0.510) GeV. " \
  "Mass window rejections: K_S0 veto (0.468,0.528)/(0.438,0.538), " \
  "eta veto (0.498,0.578), omega veto (0.732,0.832), " \
  "eta' veto (0.908,1.008), phi veto (0.990,1.390) GeV. " \
  "DeltaE_sig requirements per mode (Table 1). " \
  "Opening angle between signal D and tag D > 160 deg. " \
  "DTag stores but does not cut mBC and deltaE; 2D fit in ROOT. " \
  "Best tag candidate: min |DeltaE_tag|. " \
  "Peaking background veto: D->K pi pi0 pi0 removed for K pi pi0 eta modes."
)

alg_Dplus_tag.apply
alg_Dplus_tag.execute_on([psip3770_data, psip3770_incMC])

# ===== TagAnalysis 2: ST anti-D0 (tag D0), signal D+ -> eta + hadrons =====
alg_D0_tag = TagAnalysis.new("Dplus_to_eta_DT_D0Tag")
alg_D0_tag.set_header(["Dplus_to_eta_DT_D0TagAlg/Dplus_to_eta_DT_D0Tag.h"])
alg_D0_tag.set_constant({ "ECMS" => [:double, 3.773] })
alg_D0_tag.with_decay_card(decay_card)

# ST tag anti-D0 (D0 decays): 3 modes
alg_D0_tag.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

# Signal D+ from remaining tracks
alg_D0_tag.signal_side do |s|
  s.photons 2
  s.charged(at_least: true)
  s.min_photon_angle 10.0
end

alg_D0_tag.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_D0_tag.dtag_reconstruction do |d|
  d.beam_energy [:constant, 1.8865]
end

# 5 D+ signal modes available in DTagAlg:
# DptoKsPiEta, DptoKsKEta, DptoKPiPiEta, DptoKsPiPi0Eta, DptoPiPiPiEta
# 1 D+ mode NOT in DTagAlg (DptoPiPi0Pi0Eta) — dropped
alg_D0_tag.note(:tag_mode_unavailable,
  "1 of 6 D+ signal modes not in DTagAlg: D+ -> pi+ pi0 pi0 eta. " \
  "Available: DptoKsPiEta, DptoKsKEta, DptoKPiPiEta, DptoKsPiPi0Eta, DptoPiPiPiEta."
)

alg_D0_tag.note(:selection_details,
  "Same selection details as D0 signal analysis above. " \
  "ST D- yields: N_ST_D- = 1558159+/-2113. " \
  "ST anti-D0 yields: N_ST_D0bar = 2327839+/-1860."
)

alg_D0_tag.apply
alg_D0_tag.execute_on([psip3770_data, psip3770_incMC])