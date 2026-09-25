### Dataset preparation: seven c.m. energy points (4.59953 – 4.69882 GeV) ###
# Real data scan samples
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV
# Matched inclusive MC samples
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

scan_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
scan_incMCs = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card covering all three signal modes (and their charge-conjugate chains).
# Continuum production e+e- -> Lambda_c+ anti-Lambda_c- : psi(4260) used as KKMC top mother.
decay_card_lambdac = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-    PHSP;
    Enddecay

    Decay Lambda_c+
    0.3333 Lambda0 K_S0 K+    PHSP;
    0.3333 Sigma0  K_S0 K+    PHSP;
    0.3333 Xi0     K_S0 pi+   PHSP;
    Enddecay

    Decay anti-Lambda_c-
    0.3333 anti-Lambda0 K_S0 K-    PHSP;
    0.3333 anti-Sigma0  K_S0 K-    PHSP;
    0.3333 anti-Xi0     K_S0 pi-   PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-          HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+     HypWK;
    Enddecay

    Decay Sigma0
    1.0000 Lambda0 gamma   PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 anti-Lambda0 gamma   PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0     PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0     PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma     PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event signal MC, same card for every energy point of the scan
exMC_lambdac = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "lambdac_signal_three_modes"
  config.events        = 200000
  config.decay_card    = decay_card_lambdac
  config.cross_section = :default
end

### Event selection (tag-based: TagAnalysis) ###
# Double-tag Lambda_c: hadronic tag side anti-Lambda_c- (charm -1, all tag modes),
# signal side reconstructed from the tracks/showers the tag did not use.
# One TagAnalysis per signal mode (different missing particle / signal content).
# The unreconstructed Lambda / Sigma0 / Xi0 are handled as a missing-particle
# constraint inside the tag fit (massive form).

# ---------- Mode I: Lambda_c+ -> Lambda K_S0 K+ (Lambda missing) ----------
alg_modeI = TagAnalysis.new("LambdacTagModeI")
alg_modeI.set_header(["LambdacTagModeIAlg/LambdacTagModeI.h"])
         .set_constant({"ECMS" => [:double, 4.600]})
         .with_decay_card(decay_card_lambdac)

alg_modeI.tag_side(:Lambdac) do |t|
  t.modes :all          # all tag modes
  t.charm -1            # tag side anti-Lambda_c-
end

alg_modeI.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)   # K+ plus the K_S0 daughters pi+ pi-
  s.require_charge 1                 # net charge +1
  s.missing :Lambda                  # Lambda unreconstructed (nominal mass 1.115683 GeV)
end

alg_modeI.fit do |f|
  f.constrain_four_momentum          # 4-momentum conservation against measured CMS
  f.chi2_cut 200                     # keep chi2 < 200
end

alg_modeI
  .note(:k_s0_selection, "one K_S0 reconstructed from pi+ pi- via a secondary-vertex fit; require decay-length significance L/sigma_L > 2 before the tag fit (inexpressible in the tag signal-side surface)")
  .note(:bachelor_pion_pid, "bachelor pi+ identified by likelihood selection requiring L(pi) > L(K) and L(pi) > L(p)")
  .note(:charge_conjugate_modes, "charge-conjugate signal modes (anti-Lambda_c- -> ...) are included in addition to the Lambda_c+ modes")

alg_modeI.apply
alg_modeI.execute_on(scan_points + scan_incMCs + exMC_lambdac)

# ---------- Mode II: Lambda_c+ -> Sigma0 K_S0 K+ (Sigma0 missing) ----------
alg_modeII = TagAnalysis.new("LambdacTagModeII")
alg_modeII.set_header(["LambdacTagModeIIAlg/LambdacTagModeII.h"])
          .set_constant({"ECMS" => [:double, 4.600]})
          .with_decay_card(decay_card_lambdac)

alg_modeII.tag_side(:Lambdac) do |t|
  t.modes :all
  t.charm -1
end

alg_modeII.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :Sigma0                  # Sigma0 unreconstructed (nominal mass 1.192642 GeV)
end

alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_modeII
  .note(:k_s0_selection, "one K_S0 reconstructed from pi+ pi- via a secondary-vertex fit; require decay-length significance L/sigma_L > 2 before the tag fit")
  .note(:bachelor_pion_pid, "bachelor pi+ identified by likelihood selection requiring L(pi) > L(K) and L(pi) > L(p)")
  .note(:charge_conjugate_modes, "charge-conjugate signal modes are included")

alg_modeII.apply
alg_modeII.execute_on(scan_points + scan_incMCs + exMC_lambdac)

# ---------- Mode III: Lambda_c+ -> Xi0 K_S0 pi+ (Xi0 missing) ----------
alg_modeIII = TagAnalysis.new("LambdacTagModeIII")
alg_modeIII.set_header(["LambdacTagModeIIIAlg/LambdacTagModeIII.h"])
           .set_constant({"ECMS" => [:double, 4.600]})
           .with_decay_card(decay_card_lambdac)

alg_modeIII.tag_side(:Lambdac) do |t|
  t.modes :all
  t.charm -1
end

alg_modeIII.signal_side do |s|
  s.charged(pip: 2, pim: 1)          # bachelor pi+ plus the K_S0 daughters pi+ pi-
  s.require_charge 1                 # net charge +1
  s.missing :Xi0                     # Xi0 unreconstructed (nominal mass 1.31486 GeV)
end

alg_modeIII.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_modeIII
  .note(:k_s0_selection, "one K_S0 reconstructed from pi+ pi- via a secondary-vertex fit; require decay-length significance L/sigma_L > 2 before the tag fit")
  .note(:bachelor_pion_pid, "bachelor pi+ identified by likelihood selection requiring L(pi) > L(K) and L(pi) > L(p)")
  .note(:charge_conjugate_modes, "charge-conjugate signal modes are included")

alg_modeIII.apply
alg_modeIII.execute_on(scan_points + scan_incMCs + exMC_lambdac)