# 2202.13601v2: D0 -> K_L0 phi, K_L0 eta, K_L0 omega, K_L0 eta' BF measurement at psi(3770)
# Double-tag (DT) TagAnalysis with K_L0 as missing particle
# Multi-channel: separate TagAnalysis per sub-decay mode

# --- Common setup ---
HEADERS = ['kkmpw_evtgen/EvtGenMy/KKMCParticle.h']

DECAY_CARD = <<~DECAY
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 phi PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay phi
1.0 K+ K- VSS;
Enddecay
DECAY

# --- Real data and inclusive MC ---
real_3773 = DatasetManager.load_real_data.find('712_3773')
inc_3773 = DatasetManager.load_inclusive_mc.find('712_3773')

# --- Exclusive MC (per sub-decay) ---
# phi channel
mc_phi = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~PHI
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 phi PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay phi
1.0 K+ K- VSS;
Enddecay
  PHI
  c.sample_name '712_3773'
end

# eta -> gamma gamma channel
mc_eta_gg = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~ETA_GG
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 eta PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay eta
1.0 gamma gamma PHSP;
Enddecay
  ETA_GG
  c.sample_name '712_3773'
end

# eta -> pi+ pi- pi0 channel
mc_eta_3pi = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~ETA_3PI
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 eta PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay eta
1.0 pi+ pi- pi0 PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
  ETA_3PI
  c.sample_name '712_3773'
end

# omega -> pi+ pi- pi0 channel
mc_omega = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~OMEGA
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 omega PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay omega
1.0 pi+ pi- pi0 PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
  OMEGA
  c.sample_name '712_3773'
end

# eta' -> pi+ pi- eta channel
mc_etap = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~ETAP
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 eta' PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay eta'
1.0 pi+ pi- eta PHSP;
Enddecay
Decay eta
1.0 gamma gamma PHSP;
Enddecay
  ETAP
  c.sample_name '712_3773'
end

# eta' -> gamma rho0 channel
mc_etap_grho = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~ETAP_GRHO
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K_L0 eta' PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay eta'
1.0 gamma rho0 PHSP;
Enddecay
Decay rho0
1.0 pi+ pi- PHSP;
Enddecay
  ETAP_GRHO
  c.sample_name '712_3773'
end

# =====================================================================
# Channel 1: D0 -> K_L0 phi, phi -> K+ K-
# =====================================================================
alg_phi = TagAnalysis.new("D0_to_KL0_phi_KK")
alg_phi.set_header(HEADERS)
alg_phi.set_constant({ ECMS: 3.773 })

alg_phi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_phi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_phi.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.missing :K_L0
end

alg_phi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_phi.apply
alg_phi.execute_on([real_3773, inc_3773, mc_phi])

# =====================================================================
# Channel 2a: D0 -> K_L0 eta, eta -> gamma gamma
# =====================================================================
alg_eta_gg = TagAnalysis.new("D0_to_KL0_eta_gg")
alg_eta_gg.set_header(HEADERS)
alg_eta_gg.set_constant({ ECMS: 3.773 })

alg_eta_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_eta_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_eta_gg.signal_side do |s|
  s.photons 2
  s.missing :K_L0
end

alg_eta_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta_gg.apply
alg_eta_gg.execute_on([real_3773, inc_3773, mc_eta_gg])

# =====================================================================
# Channel 2b: D0 -> K_L0 eta, eta -> pi+ pi- pi0
# =====================================================================
alg_eta_3pi = TagAnalysis.new("D0_to_KL0_eta_3pi")
alg_eta_3pi.set_header(HEADERS)
alg_eta_3pi.set_constant({ ECMS: 3.773 })

alg_eta_3pi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_eta_3pi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_eta_3pi.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.missing :K_L0
end

alg_eta_3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_eta_3pi.apply
alg_eta_3pi.execute_on([real_3773, inc_3773, mc_eta_3pi])

# =====================================================================
# Channel 3: D0 -> K_L0 omega, omega -> pi+ pi- pi0
# =====================================================================
alg_omega = TagAnalysis.new("D0_to_KL0_omega")
alg_omega.set_header(HEADERS)
alg_omega.set_constant({ ECMS: 3.773 })

alg_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_omega.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.missing :K_L0
end

alg_omega.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_omega.apply
alg_omega.execute_on([real_3773, inc_3773, mc_omega])

# =====================================================================
# Channel 4a: D0 -> K_L0 eta', eta' -> pi+ pi- eta (eta -> gamma gamma)
# =====================================================================
alg_etap = TagAnalysis.new("D0_to_KL0_etap_pipi_eta")
alg_etap.set_header(HEADERS)
alg_etap.set_constant({ ECMS: 3.773 })

alg_etap.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_etap.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_etap.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.missing :K_L0
end

alg_etap.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_etap.apply
alg_etap.execute_on([real_3773, inc_3773, mc_etap])

# =====================================================================
# Channel 4b: D0 -> K_L0 eta', eta' -> gamma rho0 (rho0 -> pi+ pi-)
# =====================================================================
alg_etap_grho = TagAnalysis.new("D0_to_KL0_etap_gamma_rho")
alg_etap_grho.set_header(HEADERS)
alg_etap_grho.set_constant({ ECMS: 3.773 })

alg_etap_grho.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_etap_grho.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv
end

alg_etap_grho.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 1
  s.missing :K_L0
end

alg_etap_grho.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap_grho.apply
alg_etap_grho.execute_on([real_3773, inc_3773, mc_etap_grho])