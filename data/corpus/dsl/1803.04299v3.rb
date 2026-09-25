# BESIII: Measurements of Absolute Branching Fractions for
# Lambda_c+ -> Xi0 K+ and Xi(1530)0 K+
# Data: 567 pb^-1 at sqrt(s)=4.6 GeV; double-tag technique.

### Dataset description ###
data_4600   = DatasetManager.real_data.find("703_4600")
incMC_4600  = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for signal Lambda_c+ -> Xi0 K+
decay_card_Xi0K = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-          PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Xi0 K+                            PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0                       PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-                            HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                       PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0                      PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                           PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for signal Lambda_c+ -> Xi(1530)0 K+
decay_card_Xis0K = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-          PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Xi*0 K+                           PHSP;
  Enddecay

  Decay Xi*0
  1.0000 Xi0 pi0                           PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0                       PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-                            HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                       PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0                      PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                           PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for Lambda_c+ -> Xi0 K+ channel
exMC_Xi0K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lambdac_Xi0K_exclusive_mc"
  config.related_dataset = data_4600
  config.events = 200000
  config.decay_card = decay_card_Xi0K
  config.cross_section = :default
end

# Exclusive MC for Lambda_c+ -> Xi(1530)0 K+ channel
exMC_Xis0K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lambdac_Xis0K_exclusive_mc"
  config.related_dataset = data_4600
  config.events = 200000
  config.decay_card = decay_card_Xis0K
  config.cross_section = :default
end

# =============================================================================
# TagAnalysis: Lambda_c+ -> Xi0 K+  (ST + missing)
# Double-tag technique: ST anti-Lambda_c- tag + signal side K+ + missing Xi0
# =============================================================================
tag_Xi0K = TagAnalysis.new("LambdacXi0K", '00-00-01')
tag_Xi0K.set_header(["LambdacXi0KAlg/LambdacXi0K.h"])
          .set_constant({"ECMS" => [:double, 4.5995]})

# Tag side: reconstruct anti-Lambda_c- via 10 hadronic modes (12 in paper,
# 2 modes unavailable: Sigma_bar- pi0 and Sigma_bar- pi- pi+ have no
# corresponding DTagAlg channel symbols)
tag_Xi0K.tag_side(:Lambdac) do |t|
  t.modes(
    :LambdacPtoKsP,              # pbar K_S0
    :LambdacPtoKPiP,             # pbar K+ pi-
    :LambdacPtoKsPi0P,           # pbar K_S0 pi0
    :LambdacPtoKsPiPiP,          # pbar K_S0 pi- pi+
    :LambdacPtoKPiPi0P,          # pbar K+ pi- pi0
    :LambdacPtoPiPiP,            # pbar pi- pi+
    :LambdacPtoLambdaPi,         # Lambda_bar pi-
    :LambdacPtoLambdaPiPi0,      # Lambda_bar pi- pi0
    :LambdacPtoLambdaPiPiPi,     # Lambda_bar pi- pi+ pi-
    :LambdacPtoPiSIGMA0LambdaGam # Sigma_bar0 pi-
  )
  t.charm(-1)   # tag anti-Lambda_c-
  # DeltaE window: mode-dependent, ~3 sigma resolution (paper Table 2)
  # mBC is stored unconditionally, cut at ROOT level
end

# Signal side: one K+ track opposite the tag; missing Xi0 inferred from missing mass
tag_Xi0K.signal_side do |s|
  s.charged(kp: 1)
  s.missing :Xi0, mass: 1.31486   # nominal Xi0 mass (GeV/c^2)
  s.require_charge 1              # total signal-side charge = +1 (K+)
end

# Kinematic fit: constrain four-momentum of tag + signal + missing
tag_Xi0K.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_Xi0K.with_decay_card(decay_card_Xi0K).apply
tag_Xi0K
  .note(:tag_mode_unavailable, "Paper ST modes Sigma_bar- pi0 and Sigma_bar- pi- pi+ have no corresponding DTagAlg channel symbols in the current BOSS version; only 10 of 12 ST modes included.")
  .note(:st_selection, "ST anti-Lambda_c- reconstructed in 10 hadronic modes. Charged tracks: |cos(theta)|<0.93, Vz<10cm, Vr<1cm. PID: L(p)>L(K) and L(p)>L(pi) for protons; L(K)>L(pi) for kaons; L(pi)>L(K) for pions. Photons: E>25MeV barrel, E>50MeV endcap, EMC time 0-700ns. pi0 mass window [115,150] MeV/c^2 with 1C kinematic mass constraint fit. K_S0/Lambda: vertex fit chi2<100, decay vertex > 2*resolution from IP. K_S0 mass [487,511] MeV/c^2; Lambda_bar mass [1111,1121] MeV/c^2.")
  .note(:deltaE_cuts, "DeltaE cuts are mode-dependent (paper Table 2): approx 3*resolution, ranges span [-50,30] to [-20,20] MeV.")
  .note(:mBC_signal_region, "ST yield extracted from M_BC fit in signal region [2.282, 2.291] GeV/c^2.")
  .note(:dt_sel_signal, "DT: one K+ selected from remaining tracks opposite the ST anti-Lambda_c-. No multiple DT candidates. Missing mass M_miss = sqrt((E_beam - E_K+)^2 - |p_Lambda_c+ - p_K+|^2) used to infer Xi0/Xi*0. DT yields from unbinned maximum likelihood fit: N_DT_XiK = 68.2 +/- 9.9 (10.3 sigma), N_DT_XisK = 59.5 +/- 11.7 (6.4 sigma).")
  .note(:branching_fraction, "B(Lambda_c+ -> Xi0 K+) = (5.90 +/- 0.86 +/- 0.39)*10^-3. Systematic 6.7% total (MC model 3.2%, tracking 1.0%, PID 1.0%, fitting 5.2%, ST peaking bkg 0.8%, M_BC resolution 2.2%).")
tag_Xi0K.execute_on([data_4600, incMC_4600, exMC_Xi0K])

# =============================================================================
# TagAnalysis: Lambda_c+ -> Xi(1530)0 K+  (ST + missing)
# Same selection, different signal MC
# =============================================================================
tag_Xis0K = TagAnalysis.new("LambdacXis0K", '00-00-01')
tag_Xis0K.set_header(["LambdacXis0KAlg/LambdacXis0K.h"])
           .set_constant({"ECMS" => [:double, 4.5995]})

tag_Xis0K.tag_side(:Lambdac) do |t|
  t.modes(
    :LambdacPtoKsP,
    :LambdacPtoKPiP,
    :LambdacPtoKsPi0P,
    :LambdacPtoKsPiPiP,
    :LambdacPtoKPiPi0P,
    :LambdacPtoPiPiP,
    :LambdacPtoLambdaPi,
    :LambdacPtoLambdaPiPi0,
    :LambdacPtoLambdaPiPiPi,
    :LambdacPtoPiSIGMA0LambdaGam
  )
  t.charm(-1)
end

tag_Xis0K.signal_side do |s|
  s.charged(kp: 1)
  s.missing :Xi_star0, mass: 1.5318   # nominal Xi(1530)0 mass (GeV/c^2)
  s.require_charge 1
end

tag_Xis0K.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_Xis0K.with_decay_card(decay_card_Xis0K).apply
tag_Xis0K
  .note(:tag_mode_unavailable, "Paper ST modes Sigma_bar- pi0 and Sigma_bar- pi- pi+ have no corresponding DTagAlg channel symbols; only 10 of 12 ST modes included.")
  .note(:st_selection, "Same ST selection and systematics as Lambda_c+ -> Xi0 K+ analysis. See LambdacXi0K for details.")
  .note(:branching_fraction, "B(Lambda_c+ -> Xi(1530)0 K+) = (5.02 +/- 0.99 +/- 0.31)*10^-3. Systematic 6.1% total (MC model 3.9%, tracking 1.0%, PID 1.0%, fitting 3.7%, ST peaking bkg 0.8%, M_BC resolution 2.4%).")
tag_Xis0K.execute_on([data_4600, incMC_4600, exMC_Xis0K])