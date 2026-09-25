### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data at 3.773 GeV (BOSS 712)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Signal process for the exclusive MC: psi(3770) -> D0 anti-D0 (the D-tag machinery selects the tag modes)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    End
DECAYCARD

# 600k-event exclusive MC sample of psi(3770) -> D0 D0bar
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0D0bar_bnv_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events          = 600000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Alg I: D0 -> anti-p e+  (tag the opposite D0bar in the hadronic modes) ###
alg_I = TagAnalysis.new("BNV_D0ToPbarep")
alg_I.set_header(["BNV_D0ToPbarepAlg/BNV_D0ToPbarep.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .with_decay_card(decay_card_signal)
     .note(:signal_pid, "Signal-side particle ID combines dE/dx, TOF and EMC into confidence "
                        "levels: the positron requires CLe > 0.001, CLe/(CLe+CLpi+CLK) > 0.8 and "
                        "E/p > 0.85c; the antiproton requires CLp > 0.001 with CLp > CLK and "
                        "CLp > CLpi. These thresholds are not DSL-tunable (fixed v1 signal-side PID defaults).")
     .note(:fsr_recovery, "Final-state-radiation recovery: EMC clusters within 10 deg of the "
                          "positron are added to its four-momentum before the 4C kinematic fit.")
     .note(:background_veto, "Missing-energy veto on U_miss: reject events with U_miss outside "
                             "(-0.15, 0.15) GeV to suppress the D0 -> K- e+ nu_e background.")

# Tag side: anti-D0 reconstructed in K+pi-, K+pi-pi0, K+pi-pi+pi- (both tag charges scanned for charge conjugation)
alg_I.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: one antiproton + one positron, net charge zero
alg_I.signal_side do |s|
  s.charged(prm: 1, ep: 1)
  s.require_charge 0
end

# Four-momentum-constrained kinematic fit with chi2 < 200
alg_I.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_I.apply
alg_I.execute_on([psi3770_data, psi3770_incMC, exMC_signal])

### Alg II: D0bar -> p e-  (charge-conjugate signal channel) ###
alg_II = TagAnalysis.new("BNV_D0ToPem")
alg_II.set_header(["BNV_D0ToPemAlg/BNV_D0ToPem.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_signal)
      .note(:signal_pid, "Signal-side particle ID combines dE/dx, TOF and EMC into confidence "
                         "levels: the electron requires CLe > 0.001, CLe/(CLe+CLpi+CLK) > 0.8 and "
                         "E/p > 0.85c; the proton requires CLp > 0.001 with CLp > CLK and "
                         "CLp > CLpi. These thresholds are not DSL-tunable (fixed v1 signal-side PID defaults).")
      .note(:fsr_recovery, "Final-state-radiation recovery: EMC clusters within 10 deg of the "
                           "electron are added to its four-momentum before the 4C kinematic fit.")
      .note(:background_veto, "Missing-energy veto on U_miss: reject events with U_miss outside "
                              "(-0.15, 0.15) GeV to suppress the charge-conjugate semileptonic background.")

# Tag side: the opposite D0 reconstructed in the same three hadronic modes
alg_II.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: one proton + one electron, net charge zero
alg_II.signal_side do |s|
  s.charged(prp: 1, em: 1)
  s.require_charge 0
end

# Four-momentum-constrained kinematic fit with chi2 < 200
alg_II.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_II.apply
alg_II.execute_on([psi3770_data, psi3770_incMC, exMC_signal])