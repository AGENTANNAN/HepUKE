# ============================================================================
# Semileptonic D0 -> pi- pi0 e+ nu_e at psi(3770) with DT method
# 2409.04276v1  —  BESIII
# ============================================================================

psip_data = DatasetManager.real_data.find("712_3773")
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi- pi0 e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_pi_pi0_e_nu"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("D0ToPiPi0ENu")
alg.set_header(["D0ToPiPi0ENuAlg/D0ToPiPi0ENu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card)

alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,
          :D0toKsPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg.signal_side do |s|
  s.photons 2
  s.charged(pim: 1, ep: 1)
  s.min_photon_angle 10.0
  s.missing :nu_e
  s.require_charge 0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)
  f.chi2_cut 200
end

alg
  .note(:pi0_1c_kinfit, "1-C kinematic fit constraining gamma gamma to nominal pi0 mass required chi2 < 50 before the main fit; tracks not originating from K_S0: |Vz| < 10 cm, |Vxy| < 1 cm; K_S0 candidates: vertex fit chi2 < 100, decay length > 2 sigma, |M_pi+pi- - m_KS0| < 12 MeV")
  .note(:background_veto, "E_max_extra_gamma < 0.25 GeV; no extra charged track; no extra pi0 from unused photons; M_pi-pi0e+ < 1.70 GeV; K veto: |M_pi-pi0 - m_K| < 70 MeV (BF measurement only); |U_miss| < 0.03 GeV for amplitude analysis; pi0 candidates with both photons from EMC endcap rejected; EMC time in [0, 700] ns; best pi0 selected by gamma-gamma invariant mass closest to PDG pi0 mass")
  .note(:photon_selection, "barrel: E > 25 MeV, |cos(theta)| < 0.80; endcap: E > 50 MeV, 0.86 < |cos(theta)| < 0.92; angle to nearest charged track > 10 deg")
  .note(:e_pid, "positron PID: L_e/(L_e+L_pi+L_K) > 0.8, L_e > 0.001, E/(c*p_e+) > 0.8; dE/dx and TOF used for pi/K separation: L_pi > L_K for pion, L_K > L_pi for kaon")
  .apply

alg.execute_on([psip_data, psip_incMC, exMC])