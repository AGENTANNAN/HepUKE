# 2203.01555v2: D0 -> K+ pi- pi0 (DCS) BF measurement and search for D0 -> K+ pi- pi0 pi0
# Semileptonic tag: D0bar -> K+ e- nu_e_bar (signal side has pi- K+ missing nu_e)
# NB: The semileptonic tag means the tag side is the signal D0 -> K+ pi- pi0 (DCS),
# and the "tag" side actually tags through the partner's semileptonic decay.
# In the semileptonic-tag method: one D0 is tagged via D0bar -> K+ e- nu_e_bar,
# the signal D0 -> K+ pi- pi0 is found on the recoiling side.
#
# Since the DSL TagAnalysis has the tag side reconstructing D0 through DTagAlg modes
# and the signal side declaring what the tag did NOT use, we model this as:
#   - tag_side: D0 tag through hadronic modes (D0bar -> K+pi-, K+pi-pi0, K+pi-pi-pi+)
#   - signal_side: K+ pi- pi0 / K+ pi- pi0 pi0 with the D0 missing
#
# Actually, reading the paper more carefully: this uses the SEMILEPTONIC tag method.
# Tag: anti-D0 -> K+ e- nu_e_bar (semileptonic decay)
# Signal: D0 -> K+ pi- pi0 or D0 -> K+ pi- pi0 pi0 (DCS decays)
#
# In the DSL TagAnalysis, the signal_side declares what remains after tag consumption.
# With semileptonic tag, the tag consumes one D0 through hadronic modes,
# and the signal side has charged tracks + photons.
#
# Alternative interpretation: The paper tags through D0bar -> K+ e- nu_e_bar
# and finds D0 -> K+ pi- pi0 on recoil. In DSL terms:
#   tag_side: the DTagAlg hadronic tag partner
#   signal_side: charged tracks + pi0 from DCS decay
# Since the tag is hadronic (not semileptonic in DTagAlg), we use:
#   - tag_side(:D0) with hadronic modes
#   - signal_side with the DCS decay products

# --- Common setup ---
HEADERS = ['kkmpw_evtgen/EvtGenMy/KKMCParticle.h']

# --- Real data and inclusive MC ---
real_3773 = DatasetManager.load_real_data.find('712_3773')
inc_3773 = DatasetManager.load_inclusive_mc.find('712_3773')

# --- Exclusive MC for D0 -> K+ pi- pi0 (DCS) ---
mc_kpi_pi0 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~DECAY
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K+ pi- pi0 PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
  DECAY
  c.sample_name '712_3773'
end

# --- Exclusive MC for D0 -> K+ pi- pi0 pi0 (DCS) ---
mc_kpi_pi0pi0 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~DECAY
Decay psi(3770)
1.0 D0 anti-D0 PHSP;
Enddecay
Decay D0
1.0 K+ pi- pi0 pi0 PHSP;
Enddecay
Decay anti-D0
1.0 K+ pi- PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
  DECAY
  c.sample_name '712_3773'
end

# =====================================================================
# Channel 1: D0 -> K+ pi- pi0 (DCS decay, hadronic tag)
# =====================================================================
alg_kpi_pi0 = TagAnalysis.new("D0_to_Kpipi0_DCS")
alg_kpi_pi0.set_header(HEADERS)
alg_kpi_pi0.set_constant({ ECMS: 3.773 })

alg_kpi_pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_kpi_pi0.signal_side do |s|
  s.charged(kp: 1, pim: 1)
  s.photons 2
end

alg_kpi_pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_kpi_pi0.apply
alg_kpi_pi0.execute_on([real_3773, inc_3773, mc_kpi_pi0])

# =====================================================================
# Channel 2: D0 -> K+ pi- pi0 pi0 (DCS decay search, hadronic tag)
# =====================================================================
alg_kpi_pi0pi0 = TagAnalysis.new("D0_to_Kpipi0pi0_DCS")
alg_kpi_pi0pi0.set_header(HEADERS)
alg_kpi_pi0pi0.set_constant({ ECMS: 3.773 })

alg_kpi_pi0pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

alg_kpi_pi0pi0.signal_side do |s|
  s.charged(kp: 1, pim: 1)
  s.photons 4
end

alg_kpi_pi0pi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kpi_pi0pi0.apply
alg_kpi_pi0pi0.execute_on([real_3773, inc_3773, mc_kpi_pi0pi0])

# NOTE: The paper uses semileptonic tag (D0bar -> K+ e- nu_e_bar) rather than
# hadronic tag. However, the DSL TagAnalysis currently models the tag through
# DTagAlg hadronic modes. The semileptonic tag pattern (tag consumes K+ e- nu_e_bar,
# signal side has K+ pi- pi0) would need the signal_side to express lepton keys
# (ep/em) in the tag rather than signal. This is a known limitation.
# The above code uses hadronic tag as the best available approximation.
#
# For the actual semileptonic-tag method:
#   Tag side would be: D0bar -> K+ e- nu_e_bar
#   Signal side: D0 -> K+ pi- pi0 (recoil)
# This could potentially be expressed as:
#   tag_side(:D0) with semileptonic mode (not in DTagAlg hadronic table)
#   OR: Use ordinary Algorithm + Selection with partial_rec
# Neither approach maps cleanly to the current DSL surface for semileptonic tagging.