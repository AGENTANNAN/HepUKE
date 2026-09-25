# =====================================================================
# psi(3770) D-tag analysis
#   signal D+ -> anti-K0 e+ nu_e   (anti-K0 -> K_S0 -> pi+ pi-)
#   signal D+ -> pi0 e+ nu_e       (pi0 -> gamma gamma)
#   tag   : the recoil D- (charm -1) in nine hadronic modes
# =====================================================================

# ------------------------------ datasets ------------------------------
psi3770_data  = DatasetManager.real_data.find("712_3773")     # 2.93 fb^-1 psi(3770) data
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC

# --------------------------- decay cards ------------------------------
# D+ -> anti-K0 e+ nu_e , anti-K0 -> K_S0 -> pi+ pi-  (semileptonic: PHOTOS + ISGW2)
decay_card_KsEnu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.000 anti-K0 e+ nu_e PHOTOS ISGW2;
  Enddecay

  Decay anti-K0
  1.000 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

# D+ -> pi0 e+ nu_e , pi0 -> gamma gamma  (semileptonic: PHOTOS + ISGW2)
decay_card_pi0Enu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.000 pi0 e+ nu_e PHOTOS ISGW2;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

# ------------------------ exclusive signal MC -------------------------
exMC_KsEnu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKsEnu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_KsEnu
  config.cross_section   = :default
end

exMC_pi0Enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToPi0Enu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_pi0Enu
  config.cross_section   = :default
end

# =====================================================================
# Channel 1 : D+ -> anti-K0 e+ nu_e , anti-K0 -> K_S0 -> pi+ pi-
# =====================================================================
alg_KsEnu = TagAnalysis.new("DpToKsEnu")
alg_KsEnu.set_header(["DpToKsEnuAlg/DpToKsEnu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:fsr_recovery, "FSR photons within 5 deg of the positron are added to the positron four-momentum before the kinematic fit")
         .note(:background_veto, "events with an unused photon of maximum energy E_gamma,max > 300 MeV are rejected")
         .note(:track_quality, "no explicit charged-track |cos(theta)|, Vz, Vr cuts and no PID method are encoded in the selection")
         .with_decay_card(decay_card_KsEnu)

# tag side: single tag on the D- (charm -1) in the nine hadronic modes
alg_KsEnu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsK, :DptoKKPi, :DptoKPiPiPi0,
          :DptoPiPiPi, :DptoKsPiPi0, :DptoKPiPiPiPi, :DptoKsPiPiPi
  t.charm(-1)
end

# signal side: what the tag did not use (K_S0 daughters + positron)
alg_KsEnu.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge(1)               # net charge +1
  s.min_photon_angle(10.0)
  s.missing(:nu_e)                  # massless missing neutrino
end

# kinematic fit: 4-momentum conservation + K_S0 mass constraint
alg_KsEnu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut(200)
end

alg_KsEnu.apply
alg_KsEnu.execute_on([psi3770_data, psi3770_incMC, exMC_KsEnu])

# =====================================================================
# Channel 2 : D+ -> pi0 e+ nu_e , pi0 -> gamma gamma
# =====================================================================
alg_pi0Enu = TagAnalysis.new("DpToPi0Enu")
alg_pi0Enu.set_header(["DpToPi0EnuAlg/DpToPi0Enu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:fsr_recovery, "FSR photons within 5 deg of the positron are added to the positron four-momentum before the kinematic fit")
          .note(:background_veto, "events with an unused photon of maximum energy E_gamma,max > 300 MeV are rejected")
          .note(:gamma_gamma_preselection, "the gamma gamma invariant mass is required to lie in (0.110, 0.150) GeV/c2 before the pi0 mass-constrained 1-C fit; the best pi0 candidate is the one with the smallest 1-C chi2")
          .note(:photon_timing, "photon showers are required within 700 ns of the event start")
          .note(:photon_energy_threshold, "signal photons must have E > 25 MeV (barrel) / 50 MeV (endcap)")
          .note(:track_quality, "no explicit charged-track |cos(theta)|, Vz, Vr cuts and no PID method are encoded in the selection")
          .with_decay_card(decay_card_pi0Enu)

# tag side: identical nine hadronic modes, D- (charm -1)
alg_pi0Enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsK, :DptoKKPi, :DptoKPiPiPi0,
          :DptoPiPiPi, :DptoKsPiPi0, :DptoKPiPiPiPi, :DptoKsPiPiPi
  t.charm(-1)
end

# signal side: positron + two photons of the pi0
alg_pi0Enu.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge(1)               # net charge +1
  s.photons(2)                      # pi0 -> gamma gamma
  s.min_photon_angle(10.0)
  s.min_photon_energy(0.025)        # barrel floor; endcap 50 MeV (see note)
  s.missing(:nu_e)                  # massless missing neutrino
end

# kinematic fit: 4-momentum conservation + pi0 mass constraint
alg_pi0Enu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.110, 0.150)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut(200)
end

alg_pi0Enu.apply
alg_pi0Enu.execute_on([psi3770_data, psi3770_incMC, exMC_pi0Enu])