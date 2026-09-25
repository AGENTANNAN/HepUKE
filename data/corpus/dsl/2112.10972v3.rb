# 2112.10972v3: Search for BNV decays D0 → pbar e+ and D0 → p e-
# √s = 3.773 GeV, ψ(3770) → D0 D0bar, double-tag method
# Tag: anti-D0 reconstructed via hadronic modes.
# Signal: D0 → anti-p e+ (Alg I) and D0 → p e- (Alg II).
# Charge-conjugate channels implied throughout.

psip3770_data  = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process ψ(3770) → D0 anti-D0
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_bnv"
  config.related_dataset = psip3770_data
  config.events          = 600_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ── Algorithm I: D0 → anti-p e+ ──────────────────────────────────
alg_I = TagAnalysis.new("D0toPbarEp")
alg_I.set_header(["D0toPbarEpAlg/D0toPbarEp.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card)

# Tag side: anti-D0 tagged via three hadronic modes
alg_I.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0 → anti-p e+  (anti-proton + positron)
alg_I.signal_side do |s|
  s.charged(prm: 1, ep: 1)
  s.require_charge 0
end

alg_I.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_I
  .note(:electron_pid, "e+ PID: combine dE/dx, TOF, EMC → confidence levels CLe, CLπ, CLK. Require CLe > 0.001 and CLe/(CLe+CLπ+CLK) > 0.8. Also require E/p > 0.85c.")
  .note(:proton_pid, "anti-proton PID: combine dE/dx, TOF → CLp, CLK, CLπ. Require CLp > 0.001, CLp > CLK, CLp > CLπ.")
  .note(:fsr_recovery, "FSR recovery: add 4-momenta of EMC clusters within 10° of e+ direction to the e+ 4-momentum.")
  .note(:background_veto, "U_miss veto to suppress D0→K-e+ν_e: U_miss ≡ E_miss − |p_miss|·c outside (−0.15, 0.15) GeV.")
  .note(:signal_region, "Signal region in M_BC^sig vs ΔE^sig: 1.860 < M_BC^sig < 1.872 GeV/c² and −0.028 < ΔE^sig < 0.018 GeV. Signal yield via event counting, background from sideband method.")
  .apply

alg_I.execute_on([psip3770_data, psip3770_incMC, exMC_signal])

# ── Algorithm II: D0 → p e- ─────────────────────────────────────
alg_II = TagAnalysis.new("D0toPEm")
alg_II.set_header(["D0toPEmAlg/D0toPEm.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card)

# Tag side: anti-D0 (same as Alg I)
alg_II.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0 → p e-  (proton + electron)
alg_II.signal_side do |s|
  s.charged(prp: 1, em: 1)
  s.require_charge 0
end

alg_II.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_II
  .note(:electron_pid, "e+ PID: combine dE/dx, TOF, EMC → confidence levels CLe, CLπ, CLK. Require CLe > 0.001 and CLe/(CLe+CLπ+CLK) > 0.8. Also require E/p > 0.85c.")
  .note(:proton_pid, "proton PID: combine dE/dx, TOF → CLp, CLK, CLπ. Require CLp > 0.001, CLp > CLK, CLp > CLπ.")
  .note(:fsr_recovery, "FSR recovery for positron: add 4-momenta of EMC clusters within 10° of e+ direction.")
  .note(:background_veto, "U_miss veto to suppress D0→K-e+ν_e.")
  .note(:signal_region, "Signal region in M_BC^sig vs ΔE^sig: 1.860 < M_BC^sig < 1.872 GeV/c² and −0.028 < ΔE^sig < 0.018 GeV. Signal yield via event counting.")
  .apply

alg_II.execute_on([psip3770_data, psip3770_incMC, exMC_signal])