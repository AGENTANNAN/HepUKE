# BESIII analysis: Observation of D0 -> b1(1235)- e+ nu_e and evidence for D+ -> b1(1235)0 e+ nu_e
# arXiv: 2407.20551v2
# Double-tag method at sqrt(s)=3.773 GeV, 7.9 fb^-1
# b1 -> omega pi, omega -> pi+ pi- pi0

### Datasets ###
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

### Decay cards ###
# D0 -> b1(1235)- e+ nu_e, b1- -> omega pi-, omega -> pi+ pi- pi0
decay_card_d0_b1_enu = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 b1- e+ nu_e SLN;
    Enddecay

    Decay b1-
    1.000 omega pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# D+ -> b1(1235)0 e+ nu_e, b10 -> omega pi0, omega -> pi+ pi- pi0
decay_card_dp_b1_enu = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 b10 e+ nu_e SLN;
    Enddecay

    Decay b10
    1.000 omega pi0 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC ###
exMC_d0_b1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "d0_b1_enu_signal"
  config.events        = 500_000
  config.decay_card    = decay_card_d0_b1_enu
  config.related_dataset = data_3773
  config.cross_section = :default
end

exMC_dp_b1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "dp_b1_enu_signal"
  config.events        = 500_000
  config.decay_card    = decay_card_dp_b1_enu
  config.related_dataset = data_3773
  config.cross_section = :default
end

### Tag-based Analysis 1: D0 -> b1(1235)- e+ nu_e ###
alg_d0_b1_enu = TagAnalysis.new("D0B1ENu")
alg_d0_b1_enu.set_header(["D0B1ENuAlg/D0B1ENu.h"])
  .note(:dt_method, "Double-tag method at psi(3770). ST D0bar mesons reconstructed from hadronic decays; signal D0 -> b1- e+ nu_e selected from residual tracks/showers.")
  .note(:st_tag_modes, "6 ST D0bar modes: D0->K+pi-, D0->K+pi-pi0, D0->K+pi-pi-pi+, D0->KS0pi+pi-, D0->K+pi-pi0pi0, D0->K+pi+pi-pi-pi0.")
  .note(:tag_deltae, "DeltaE window per tag mode per Table 1. Keep candidate with minimum |DeltaE| per mode per charge per event.")
  .note(:tag_mbc_fit, "ST yields from unbinned ML fit to M_BC: MC signal shape convolved with double-Gaussian, ARGUS function for background. ST window: (1.859, 1.873) GeV/c^2 for D0bar.")
  .note(:ks_veto, "K_S0 veto: |M(pi+pi-) - m(K_S0)| > 0.008 GeV/c^2 to suppress K1(1270)->KS0 pi+ pi- background.")
  .note(:omega_window, "omega -> pi+ pi- pi0 signal region: M(pi+pi-pi0) in (0.757, 0.807) GeV/c^2. Sidebands: (0.697,0.742) and (0.822,0.867).")
  .note(:a0_veto, "M(pi+pi-pi0) > 0.6 GeV/c^2 to veto D->a0(980) e+ nu_e background.")
  .note(:positron_pid, "e+ PID: CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8. E/p > 0.8 for EMC energy deposit.")
  .note(:b1_mass_cut, "M(b1 e+) < 1.82 GeV/c^2 to suppress D->b1- pi+ peaking background.")
  .note(:extra_shower_veto, "E_extra_gamma_max < 0.30 GeV, N_extra_pi0 = 0, cos_theta(gamma,miss) < 0.3 (D0)/0.4 (D+).")
  .note(:signal_extraction, "2D unbinned ML fit to M(omega pi) vs U_miss. U_miss = E_miss - |p_miss|*c. Signal peaks at U_miss ~ 0 and M(omega pi) ~ b1 mass.")
  .note(:umiss, "U_miss computed using beam energy and known D mass for tag-side momentum: p_D = -hat{p}_Dbar * sqrt(E_beam^2/c^2 - m_Dbar^2*c^2).")

alg_d0_b1_enu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toK3Pi, :D0toKsPiPi,
          :D0toKPiPi0Pi0, :D0toK3PiPi0
  t.charm 1
end

alg_d0_b1_enu.signal_side do |s|
  s.photons 2
  s.charged(ep: 1, pip: 1, pim: 2)
  s.require_charge 0
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_d0_b1_enu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_d0_b1_enu.with_decay_card(decay_card_d0_b1_enu).apply
alg_d0_b1_enu.execute_on([data_3773] + [incMC_3773] + [exMC_d0_b1])

### Tag-based Analysis 2: D+ -> b1(1235)0 e+ nu_e ###
alg_dp_b1_enu = TagAnalysis.new("DpB1ENu")
alg_dp_b1_enu.set_header(["DpB1ENuAlg/DpB1ENu.h"])
  .note(:dt_method, "Double-tag method at psi(3770). ST D- mesons from hadronic decays; signal D+ -> b10 e+ nu_e.")
  .note(:st_tag_modes, "6 ST D- modes: D-->K+pi-pi-, D-->KS0pi-, D-->K+pi-pi-pi0, D-->KS0pi-pi0, D-->KS0pi+pi-pi-, D-->K+K-pi-.")
  .note(:tag_deltae, "DeltaE window per tag mode. Best candidate by minimum |DeltaE| per mode per charge per event.")
  .note(:tag_mbc_fit, "ST yields from M_BC fit. ST window: (1.863, 1.877) GeV/c^2 for D-.")
  .note(:ks_veto, "K_S0 veto: |M(pi+pi-) - m(K_S0)| > 0.008 GeV/c^2.")
  .note(:omega_window, "omega signal region: M(pi+pi-pi0) in (0.757, 0.807) GeV/c^2. Sidebands: (0.697,0.742) and (0.822,0.867).")
  .note(:a0_veto, "M(pi+pi-pi0) > 0.6 GeV/c^2 for a0 veto.")
  .note(:positron_pid, "e+ PID with CL_e > 0.001, CL_e/(CL_e+CL_pi+CL_K) > 0.8, E/p > 0.8.")
  .note(:b1_mass_cut, "M(b10 e+) < 1.82 GeV/c^2.")
  .note(:extra_shower_veto, "E_extra_gamma_max < 0.30 GeV, N_extra_pi0 = 0, cos_theta(gamma,miss) < 0.4 for D+.")
  .note(:signal_extraction, "2D unbinned ML fit to M(omega pi0) vs U_miss.")
  .note(:umiss, "U_miss computed with tag-side momentum from beam constraint.")

alg_dp_b1_enu.tag_side(:Dp) do |t|
  t.modes :DmtoKPiPi, :DmtoKsPi, :DmtoKPiPiPi0,
          :DmtoKsPiPi0, :DmtoKsPiPiPi, :DmtoKKPi
  t.charm -1
end

alg_dp_b1_enu.signal_side do |s|
  s.photons 4
  s.charged(ep: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_dp_b1_enu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_b1_enu.with_decay_card(decay_card_dp_b1_enu).apply
alg_dp_b1_enu.execute_on([data_3773] + [incMC_3773] + [exMC_dp_b1])