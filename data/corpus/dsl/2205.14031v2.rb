# 2205.14031v2: Absolute BFs of hadronic D-meson decays involving kaons and pions
#
# Tag-based (D double-tag at psi(3770)). Uses 2.93 fb-1 at sqrt(s)=3.773 GeV.
# ST: D0bar -> K+pi-, K+pi-pi0, K+pi-pi-pi+ and D- -> K+pi-pi-, KS0pi-,
#      K+pi-pi-pi0, KS0pi-pi0, KS0pi+pi-pi-, K+K-pi-
# DT: 7 signal modes with D0 and D+ decays into multi-body final states
#      containing KS0, pi0, pi+, pi-, K-.
# 2D fits to MBC(tag) vs MBC(sig). QC corrections for neutral D decays.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# --- Datasets (psi(3770), BOSS 712) ---
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================
# Algorithm 1: D0 tag — D0 decays
# Signal modes:
#   D0 -> KS0 pi0 pi0 pi0
#   D0 -> K- pi+ pi0 pi0 pi0
#   D0 -> KS0 pi+ pi- pi0 pi0
# ============================================

decay_card_psipp = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  End
DECAYCARD

# --- D0 tag analysis ---
alg_d0_signal = TagAnalysis.new("D0TagMultiPi0")
alg_d0_signal.set_header(["D0TagMultiPi0Alg/D0TagMultiPi0.h"])
               .set_constant({ "ECMS" => [:double, 3.773] })
               .with_decay_card(decay_card_psipp)

# Tag side 1: anti-D0 in tag modes
alg_d0_signal.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Tag side 2: D0 signal side. For the signal modes, they share
# common reconstruction: KS0, pi+, pi-, pi0.
# We use KPiPi0Pi0Pi0 as a representative; the individual modes
# (KS0pi0pi0pi0 etc.) differ only in the final-state assignment
# handled at the NTuple level.
alg_d0_signal.tag_side(:D0) do |t|
  t.modes :D0toKPiPi0Pi0Pi0
  t.charm 1
  t.rank_by :mbc
end

# Signal side: photons from pi0 decays + charged pions/kaons
alg_d0_signal.signal_side do |s|
  s.photons 6       # up to 6 photons for up to 3 pi0
  s.charged(km: 1, pip: 1)
  s.require_charge 0
  s.min_photon_angle 10.0
end

alg_d0_signal.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_d0_signal.note(:d0_signal_modes,
  "Signal modes: D0 -> KS0 pi0 pi0 pi0, D0 -> K- pi+ pi0 pi0 pi0, " \
  "D0 -> KS0 pi+ pi- pi0 pi0. Individual mode yields extracted via " \
  "2D fits to MBC(tag) vs MBC(sig) in ROOT. " \
  "KS0 -> pi+ pi- reconstructed via secondary vertex fit with " \
  "mass window [0.486, 0.510] GeV/c2. " \
  "pi0 -> gamma gamma with mass window [0.115, 0.150] GeV/c2.")

alg_d0_signal.note(:ks0_veto_sig,
  "For D0 -> KS0 pi+ pi- pi0 pi0: veto pi+ pi- pairs in [0.468, 0.528] " \
  "GeV/c2 to suppress D0 -> KS0 KS0 pi0 pi0. " \
  "For D0 -> KS0 pi0 pi0 pi0: peaking background from " \
  "D0 -> KS0 KS0(->pi0pi0) pi0 estimated from data (KS0->pi+pi-).")

alg_d0_signal.note(:qc_correction,
  "Quantum correlation correction factors: f_QC = 1.081 +/- 0.007 " \
  "for D0 -> KS0 pi0 pi0 pi0 (CP-even), f_QC = 0.956 +/- 0.006 " \
  "for D0 -> KS0 pi+ pi- pi0 pi0 (CP-odd). Measured using " \
  "CP-even tag D0 -> K+K- and CP-odd tag D0 -> KS0 pi0.")

alg_d0_signal.dtag_reconstruction do |d|
  d.beam_energy [:constant, 1.8865]
end

alg_d0_signal.apply
alg_d0_signal.execute_on([data_3773, incMC_3773])

# ============================================
# Algorithm 2: D+ tag — D+ decays
# Signal modes:
#   D+ -> KS0 pi+ pi0 pi0
#   D+ -> KS0 pi+ pi+ pi- pi0
#   D+ -> KS0 pi+ pi0 pi0 pi0
#   D+ -> K- pi+ pi+ pi0 pi0
# ============================================

decay_card_psipp_dp = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  End
DECAYCARD

alg_dp_signal = TagAnalysis.new("DpTagMultiBody")
alg_dp_signal.set_header(["DpTagMultiBodyAlg/DpTagMultiBody.h"])
               .set_constant({ "ECMS" => [:double, 3.773] })
               .with_decay_card(decay_card_psipp_dp)

# Tag side 1: D-
alg_dp_signal.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Tag side 2: D+ signal side
alg_dp_signal.tag_side(:Dplus) do |t|
  t.modes :DptoKsPiPi0Pi0
  t.charm 1
  t.rank_by :mbc
end

alg_dp_signal.signal_side do |s|
  s.photons 6       # up to 6 photons for up to 3 pi0
  s.charged(km: 1, pip: 2, pim: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
end

alg_dp_signal.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:"D+")
  f.chi2_cut 200
end

alg_dp_signal.note(:dp_signal_modes,
  "Signal modes: D+ -> KS0 pi+ pi0 pi0, D+ -> KS0 pi+ pi+ pi- pi0, " \
  "D+ -> KS0 pi+ pi0 pi0 pi0, D+ -> K- pi+ pi+ pi0 pi0. " \
  "Yields extracted via 2D fits to MBC(tag) vs MBC(sig) in ROOT. " \
  "For D+ -> KS0 pi+ pi+ pi- pi0: KS0 veto on pi+ pi- pairs " \
  "outside [0.468, 0.528] GeV/c2. " \
  "For D+ -> KS0 pi+ pi0 pi0: peaking background from " \
  "D+ -> KS0 KS0(->pi0pi0) pi+ estimated from data.")

alg_dp_signal.note(:dt_deltaE_cuts,
  "DeltaE(sig) mode-dependent windows applied at ROOT stage: " \
  "D0 modes: (-85,+65) MeV for pi0-containing, (-35,+35) MeV else. " \
  "D+ modes: (-55,+40) MeV for pi0-containing, (-25,+25) MeV else. " \
  "Best candidate per mode selected by minimum |DeltaE(sig)|.")

alg_dp_signal.dtag_reconstruction do |d|
  d.beam_energy [:constant, 1.8865]
end

alg_dp_signal.apply
alg_dp_signal.execute_on([data_3773, incMC_3773])