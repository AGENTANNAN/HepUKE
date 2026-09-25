# ============================================================================
# Semileptonic decays D0 -> K- eta e+ nu_e, D+ -> K_S0 eta e+ nu_e,
# D+ -> eta eta e+ nu_e at psi(3770) with the DT method
# 2409.15044v3  —  BESIII
# ============================================================================

psip_data = DatasetManager.real_data.find("712_3773")
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================================
# Mode I: D0 -> K- eta e+ nu_e  (tag-based, D0 tag)
# ============================================================================

decay_card_d0_ke = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- eta e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_ke = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_K_eta_e_nu"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_d0_ke
  config.cross_section = :default
end

alg_d0_ke = TagAnalysis.new("D0ToKEtaENu")
alg_d0_ke.set_header(["D0ToKEtaENuAlg/D0ToKEtaENu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_d0_ke)

alg_d0_ke.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_d0_ke.signal_side do |s|
  s.photons 2
  s.charged(km: 1, ep: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e
  s.require_charge 0
end

alg_d0_ke.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.505, 0.575)
  f.chi2_cut 200
end

alg_d0_ke
  .note(:signal_mode, "D0 -> K- eta e+ nu_e with eta -> gamma gamma")
  .note(:pid_details, "Kaon/pion PID: L(K)>L(pi), L(pi)>L(K); e+: L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8, E/p in [0.8,1.1]")
  .note(:background_veto, "M_P_eta_e < 1.80 GeV veto; E_extra_gamma < 0.5 GeV, N_extra_pi0==0; M_K_eta > 1.30 GeV for D0->K-eta e nu")
  .note(:ks0_reconstruction, "K_S0: two oppositely charged tracks |Vz|<20, vertex fit chi2<100, mass [0.487,0.511], decay length>2sigma")
  .note(:pi0_eta_selection, "pi0: M_gg in [0.115,0.150], 1C fit chi2<50; eta: M_gg in [0.505,0.575], 1C fit chi2<50")
  .note(:photon_selection, "Barrel: E>25MeV |cos(theta)|<0.80; Endcap: E>50MeV 0.86<|cos(theta)|<0.92; TDC [0,700]ns; angle to nearest track>10deg")
  .apply

alg_d0_ke.execute_on([psip_data, psip_incMC, exMC_d0_ke])

# ============================================================================
# Mode II: D+ -> K_S0 eta e+ nu_e  (tag-based, D+ tag)
# ============================================================================

decay_card_dp_kse = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 eta e+ nu_e PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_kse = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_Ks_eta_e_nu"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_dp_kse
  config.cross_section = :default
end

alg_dp_kse = TagAnalysis.new("DpToKsEtaENu")
alg_dp_kse.set_header(["DpToKsEtaENuAlg/DpToKsEtaENu.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .with_decay_card(decay_card_dp_kse)

alg_dp_kse.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_dp_kse.signal_side do |s|
  s.photons 2
  s.charged(pip: 1, pim: 1, ep: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e
  s.require_charge 0
end

alg_dp_kse.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.505, 0.575)
  f.chi2_cut 200
end

alg_dp_kse
  .note(:signal_mode, "D+ -> K_S0 eta e+ nu_e with K_S0 -> pi+ pi-, eta -> gamma gamma")
  .note(:pid_details, "Kaon/pion PID: L(K)>L(pi), L(pi)>L(K); e+: L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8, E/p in [0.8,1.1]")
  .note(:background_veto, "M_P_eta_e < 1.80 GeV veto; E_extra_gamma < 0.5 GeV, N_extra_pi0==0")
  .note(:ks0_reconstruction, "K_S0: two oppositely charged tracks |Vz|<20, vertex fit chi2<100, mass [0.487,0.511], decay length>2sigma")
  .note(:pi0_eta_selection, "pi0: M_gg in [0.115,0.150], 1C fit chi2<50; eta: M_gg in [0.505,0.575], 1C fit chi2<50")
  .note(:photon_selection, "Barrel: E>25MeV |cos(theta)|<0.80; Endcap: E>50MeV 0.86<|cos(theta)|<0.92; TDC [0,700]ns; angle to nearest track>10deg")
  .apply

alg_dp_kse.execute_on([psip_data, psip_incMC, exMC_dp_kse])

# ============================================================================
# Mode III: D+ -> eta eta e+ nu_e  (tag-based, D+ tag)
# ============================================================================

decay_card_dp_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 eta eta e+ nu_e PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Alias eta1 eta
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    Decay eta1
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_eta_eta_e_nu"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_dp_ee
  config.cross_section = :default
end

alg_dp_ee = TagAnalysis.new("DpToEtaEtaENu")
alg_dp_ee.set_header(["DpToEtaEtaENuAlg/DpToEtaEtaENu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_dp_ee)

alg_dp_ee.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_dp_ee.signal_side do |s|
  s.photons 4
  s.charged(ep: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e
end

alg_dp_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_dp_ee
  .note(:signal_mode, "D+ -> eta eta e+ nu_e with both eta -> gamma gamma (use alias for identical daughters)")
  .note(:pid_details, "Kaon/pion PID: L(K)>L(pi), L(pi)>L(K); e+: L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8, E/p in [0.8,1.1]")
  .note(:background_veto, "M_P_eta_e < 1.80 GeV veto; E_extra_gamma < 0.5 GeV, N_extra_pi0==0")
  .note(:ks0_reconstruction, "K_S0: two oppositely charged tracks |Vz|<20, vertex fit chi2<100, mass [0.487,0.511], decay length>2sigma")
  .note(:pi0_eta_selection, "pi0: M_gg in [0.115,0.150], 1C fit chi2<50; eta: M_gg in [0.505,0.575], 1C fit chi2<50")
  .note(:photon_selection, "Barrel: E>25MeV |cos(theta)|<0.80; Endcap: E>50MeV 0.86<|cos(theta)|<0.92; TDC [0,700]ns; angle to nearest track>10deg")
  .apply

alg_dp_ee.execute_on([psip_data, psip_incMC, exMC_dp_ee])