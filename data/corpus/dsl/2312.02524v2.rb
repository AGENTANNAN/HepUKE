# Paper 2312.02524v2: Amplitude analysis of D0→π+π-π+π- and D0→π+π-π0π0
# Tag-based DT method at ψ(3770), BOSS 712

data_3773 = DatasetManager.load_real_data.find("712_3773")
inc_mc_3773 = DatasetManager.load_inclusive_mc.find("712_3773")

# === Signal Mode I: D0→π+π-π+π- (charged 4π) ===
alg_4pi_charged = TagAnalysis.new("D0to4PiCharged")
alg_4pi_charged.set_header(["D0to4PiChargedAlg/D0to4PiCharged.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })

# ST: reconstruct D0bar via 3 hadronic modes
alg_4pi_charged.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0→π+π-π+π- (4 charged pions, no photons used by signal)
alg_4pi_charged.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 2)
end

# 1C kinematic fit constraining 4π invariant mass to D0 mass
alg_4pi_charged.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Background vetoes applied in ROOT analysis (post-BOSS)
alg_4pi_charged.note(:background_veto, "K_S0 veto: reject if any π+π- pair satisfies |M(ππ)-m_Ks0|<12 MeV/c2")
alg_4pi_charged.note(:background_veto, "D-→K+π-π- veto: reject if any K+π-π- or K-π+π+ combination with |M-M_D-|<10 MeV/c2")
alg_4pi_charged.note(:background_veto, "D0→π+π-π0 veto: reject if M(π+π-π0) near D0 mass")
alg_4pi_charged.note(:background_veto, "η veto: reject if any γγ combination with |M(γγ)-m_η|<20 MeV/c2")
alg_4pi_charged.note(:extra_pi0_veto, "Reject events with additional π0 candidates beyond those forming the signal")

alg_4pi_charged.apply
alg_4pi_charged.execute_on([data_3773])

# === Signal Mode II: D0→π+π-π0π0 (neutral 4π) ===
alg_4pi_neutral = TagAnalysis.new("D0to4PiNeutral")
alg_4pi_neutral.set_header(["D0to4PiNeutralAlg/D0to4PiNeutral.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })

# ST: same 3 D0bar modes
alg_4pi_neutral.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0→π+π-π0π0 (2 charged pions + 4 photons → 2π0)
alg_4pi_neutral.signal_side do |s|
  s.photons 4
  s.charged(pip: 1, pim: 1)
end

# 1C kinematic fit constraining 4π invariant mass to D0 mass
alg_4pi_neutral.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim, :pi0, :pi0).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Background vetoes
alg_4pi_neutral.note(:background_veto, "K_S0 veto: reject if any π+π- pair satisfies |M(ππ)-m_Ks0|<12 MeV/c2")
alg_4pi_neutral.note(:background_veto, "D-→K+π-π- veto: reject if any K+π-π- or K-π+π+ combination with |M-M_D-|<10 MeV/c2")
alg_4pi_neutral.note(:background_veto, "η veto: reject if any γγ combination with |M(γγ)-m_η|<20 MeV/c2")
alg_4pi_neutral.note(:extra_pi0_veto, "Reject events with additional π0 candidates beyond the two forming the signal")
alg_4pi_neutral.note(:pi0_reconstruction, "π0 candidates reconstructed via 1C Kalman kinematic fit constraining γγ mass to π0 mass; χ2_1C<20 required")

alg_4pi_neutral.apply
alg_4pi_neutral.execute_on([data_3773])