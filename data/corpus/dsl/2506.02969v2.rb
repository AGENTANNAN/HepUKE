# Dataset preparation -- multi-energy scan points for Lambda_c+ DT analysis
# 7 energy points from 4599.53 to 4698.82 MeV
scan_datasets = [
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
]

incMCs = scan_datasets.map { |ds|
  boss = ds.boss.gsub(".", "")
  DatasetManager.inclusive_mc.find("#{boss}_#{ds.sample_name}")
}.compact

# Decay card for Lambda_c+ -> (Lambda K_S0 K+, Sigma0 K_S0 K+, Xi0 K_S0 pi+)
# and Lambda_c- charge conjugate (tag side uses anti-Lambda_c-)
decay_card_signal = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.0000 anti-Lambda K_S0 K-              PHSP;
  1.0000 anti-Sigma0 K_S0 K-              PHSP;
  1.0000 anti-Xi0 K_S0 pi-               PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda K_S0 K+                   PHSP;
  1.0000 Sigma0 K_S0 K+                   PHSP;
  1.0000 Xi0 K_S0 pi+                     PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi-                           PHSP;
  Enddecay

  Decay anti-Lambda
  1.0000 anti-p- pi+                      PHSP;
  Enddecay

  Decay Sigma0
  1.0000 Lambda gamma                     PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 anti-Lambda gamma                PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda pi0                       PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda pi0                  PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                          PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "Lambdac_DT_KS0K_Lambda_Sigma0_Xi0"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ============================================================
# Signal mode 1: Lambda_c+ -> Lambda K_S0 K+
# Tag side: anti-Lambda_c- with hadronic tag modes
# Signal side: K_S0(pi+pi-) + K+ + missing Lambda
# ============================================================
alg1 = TagAnalysis.new("LambdacDTLambdaKS0K")
alg1.set_header(["LambdacDTLambdaKS0KAlg/LambdacDTLambdaKS0K.h"])
    .set_constant({"ECMS" => [:double, 4.600]})
    .with_decay_card(decay_card_signal)

alg1.tag_side(:Lambdac) do |t|
  t.modes :all
  t.charm -1
end

alg1.signal_side do |s|
  s.charged(pim: 1, pip: 1, kp: 1)
  s.require_charge 1
  s.missing :Lambda, mass: 1.115683
end

alg1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg1
  .note(:ks0_reconstruction,
        "K_S0 reconstructed from pi+pi- via secondary vertex fit on the signal side. " \
        "L/sigma_L > 2 for decay length significance. " \
        "Applied at ROOT level; inexpressible in TagAnalysis DSL signal_side.")
  .note(:partial_reconstruction,
        "The analysis uses partial reconstruction via missing mass M_miss " \
        "to identify Lambda, Sigma0, Xi0. The TagAnalysis 'missing' declaration " \
        "approximates this via a 1C kinematic fit with AddMissTrack. " \
        "Full partial rec with mass window cuts applied at ROOT level.")
  .note(:tag_modes,
        "12 exclusive anti-Lambda_c- tag modes are used. modes :all declared; " \
        "exact mode list should be narrowed per the paper's tag mode table. " \
        "tag_modes.yaml not available for authoritative mode symbols.")
  .note(:ks0_mass_window,
        "|M(pi+pi-) - m_K_S0| selection applied at ROOT level.")
  .note(:miss_mass_window,
        "M_miss selection windows for Lambda/Sigma0/Xi0 applied at ROOT level.")

alg1.apply
alg1.execute_on(scan_datasets + incMCs + exMC_signal)

# ============================================================
# Signal mode 2: Lambda_c+ -> Sigma0 K_S0 K+
# ============================================================
alg2 = TagAnalysis.new("LambdacDTSigma0KS0K")
alg2.set_header(["LambdacDTSigma0KS0KAlg/LambdacDTSigma0KS0K.h"])
    .set_constant({"ECMS" => [:double, 4.600]})
    .with_decay_card(decay_card_signal)

alg2.tag_side(:Lambdac) do |t|
  t.modes :all
  t.charm -1
end

alg2.signal_side do |s|
  s.charged(pim: 1, pip: 1, kp: 1)
  s.require_charge 1
  s.missing :Sigma0, mass: 1.192642
end

alg2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg2
  .note(:ks0_reconstruction,
        "K_S0 reconstructed from pi+pi- via secondary vertex fit on the signal side. " \
        "L/sigma_L > 2 for decay length significance. " \
        "Applied at ROOT level; inexpressible in TagAnalysis DSL signal_side.")
  .note(:partial_reconstruction,
        "The analysis uses partial reconstruction via missing mass M_miss " \
        "to identify Sigma0. The TagAnalysis 'missing' declaration " \
        "approximates this via a 1C kinematic fit with AddMissTrack. " \
        "Full partial rec with mass window cuts applied at ROOT level.")
  .note(:tag_modes,
        "12 exclusive anti-Lambda_c- tag modes are used. modes :all declared; " \
        "exact mode list should be narrowed per the paper's tag mode table.")
  .note(:sigma0_decay,
        "Sigma0 -> Lambda gamma; the Lambda and photon are NOT reconstructed -- " \
        "Sigma0 is identified via missing mass. Inexpressible in DSL.")

alg2.apply
alg2.execute_on(scan_datasets + incMCs + exMC_signal)

# ============================================================
# Signal mode 3: Lambda_c+ -> Xi0 K_S0 pi+
# ============================================================
alg3 = TagAnalysis.new("LambdacDTXi0KS0Pi")
alg3.set_header(["LambdacDTXi0KS0PiAlg/LambdacDTXi0KS0Pi.h"])
    .set_constant({"ECMS" => [:double, 4.600]})
    .with_decay_card(decay_card_signal)

alg3.tag_side(:Lambdac) do |t|
  t.modes :all
  t.charm -1
end

alg3.signal_side do |s|
  s.charged(pim: 1, pip: 2)
  s.require_charge 1
  s.missing :Xi0, mass: 1.31486
end

alg3.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg3
  .note(:ks0_reconstruction,
        "K_S0 reconstructed from pi+pi- via secondary vertex fit on the signal side. " \
        "L/sigma_L > 2 for decay length significance. " \
        "Applied at ROOT level; inexpressible in TagAnalysis DSL signal_side.")
  .note(:partial_reconstruction,
        "The analysis uses partial reconstruction via missing mass M_miss " \
        "to identify Xi0. The TagAnalysis 'missing' declaration " \
        "approximates this via a 1C kinematic fit with AddMissTrack. " \
        "Full partial rec with mass window cuts applied at ROOT level.")
  .note(:tag_modes,
        "12 exclusive anti-Lambda_c- tag modes are used. modes :all declared; " \
        "exact mode list should be narrowed per the paper's tag mode table.")
  .note(:xi0_decay,
        "Xi0 -> Lambda pi0; the Lambda and pi0 are NOT reconstructed -- " \
        "Xi0 is identified via missing mass. Inexpressible in DSL.")
  .note(:pion_pid,
        "Bachelor pi+ PID: L(pi) > L(K) and L(pi) > L(p). " \
        "Applied via SimplePIDSvc at ROOT level.")

alg3.apply
alg3.execute_on(scan_datasets + incMCs + exMC_signal)