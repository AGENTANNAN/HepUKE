# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data (BOSS 712)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample at 3.773 GeV

# Decay card: ψ(3770) → D0 D0bar with the signal mode D0 → K−π+ (charge conjugate included), phase space.
# The CP-tag D is taken from the pre-stored DTagAlg candidates, so it is not listed in this card.
decay_card_d0d0bar = <<~DECAYCARD
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

### Event selection (BOSS) — tag-based (TagAnalysis) ###

# ===== Tag selection 1: CP-even tags (K+K−, π+π−, π0π0, K_S0π0π0) =====
alg_cp_even = TagAnalysis.new("D0CPEvenTag")
alg_cp_even.set_header(["D0CPEvenTagAlg/D0CPEvenTag.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .with_decay_card(decay_card_d0d0bar)

# Tag side: pre-stored DTagAlg candidates for the CP-even modes
# (both D0/D0bar charges are scanned, so charm is left unset; candidates ranked by invariant mass).
alg_cp_even.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toPiPi, :D0toPi0Pi0, :D0toKsPi0Pi0
end

# Signal side: the other D reconstructed as K−π+ from the tag's unused tracks.
alg_cp_even.signal_side do |s|
  s.photons 0                # require zero photons on the signal side
  s.charged(km: 1, pip: 1)   # exactly one K− and one π+
  s.require_charge(0)        # net charge zero
end

# Fit: four-momentum conservation + tag invariant mass constrained to the nominal D0 mass, χ² < 200.
alg_cp_even.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# CP tag modes that are not part of the DTagAlg tag vocabulary — handled separately.
alg_cp_even.note(:tag_modes_outside_vocabulary,
                 "CP tag modes K_L0 X, K_S0 omega, K_S0 phi, pi+pi-pi0 and K_S,L0 pi+pi- are not in the DTagAlg tag vocabulary; they are handled separately.")

alg_cp_even.apply
alg_cp_even.execute_on([data_3773, incMC_3773])

# ===== Tag selection 2: CP-odd tags (K_S0π0, K_S0η, K_S0η′) =====
alg_cp_odd = TagAnalysis.new("D0CPOddTag")
alg_cp_odd.set_header(["D0CPOddTagAlg/D0CPOddTag.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_d0d0bar)

# Tag side: pre-stored DTagAlg candidates for the CP-odd modes (K_S0η′ is via π+π−η).
alg_cp_odd.tag_side(:D0) do |t|
  t.modes :D0toKsPi0, :D0toKsEta, :D0toKsEtaPrime
end

# Signal side: the other D reconstructed as K−π+ from the tag's unused tracks.
alg_cp_odd.signal_side do |s|
  s.photons 0                # require zero photons on the signal side
  s.charged(km: 1, pip: 1)   # exactly one K− and one π+
  s.require_charge(0)        # net charge zero
end

# Fit: four-momentum conservation + tag invariant mass constrained to the nominal D0 mass, χ² < 200.
alg_cp_odd.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# CP tag modes that are not part of the DTagAlg tag vocabulary — handled separately.
alg_cp_odd.note(:tag_modes_outside_vocabulary,
                "CP tag modes K_L0 X, K_S0 omega, K_S0 phi, pi+pi-pi0 and K_S,L0 pi+pi- are not in the DTagAlg tag vocabulary; they are handled separately.")

alg_cp_odd.apply
alg_cp_odd.execute_on([data_3773, incMC_3773])