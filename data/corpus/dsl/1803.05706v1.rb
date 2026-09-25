# BESIII: Measurement of absolute branching fraction for Lambda_c+ -> Lambda + X (inclusive)
# and search for direct CP violation
# Data: 567 pb^-1 at sqrt(s)=4.6 GeV; double-tag method
# Lambda_c+ Lambda_c- pair production near threshold

### Dataset description ###
data_4600   = DatasetManager.real_data.find("703_4600")
incMC_4600  = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for Lambda_c+ -> Lambda X (inclusive, models Lambda_c+ -> Lambda + anything)
# Use Lambda_c+ -> Lambda + pi+ + pi- as representative exclusive channel for MC
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-          PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda0 pi+ pi-                   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-                            HypWK;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0                      PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                           PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal
exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lambdac_LambdaX_inclusive_mc"
  config.related_dataset = data_4600
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# =============================================================================
# TagAnalysis: Lambda_c+ -> Lambda + X (inclusive)
# ST: Lambda_c- tagged in 2 hadronic modes
# DT: Lambda reconstructed from remaining tracks opposite the tag
# =============================================================================
tag = TagAnalysis.new("LambdacLambdaX", '00-00-01')
tag.set_header(["LambdacLambdaXAlg/LambdacLambdaX.h"])
   .set_constant({"ECMS" => [:double, 4.5995]})

# Tag side: reconstruct Lambda_c- via 2 hadronic modes
tag.tag_side(:Lambdac) do |t|
  t.modes(
    :LambdacPtoKsP,         # anti-p K_S0
    :LambdacPtoKPiP         # anti-p K+ pi-
  )
  t.charm(-1)   # tag anti-Lambda_c-
  # DeltaE windows (paper Table I): mode-dependent
  # M_BC extracted from fit in signal region
end

# Signal side: Lambda reconstructed from p+ pi- among remaining tracks
# This is inclusive — any Lambda counts, no constraining the rest of X
tag.signal_side do |s|
  s.charged(prp: 1, pim: 1)     # Lambda -> p+ pi-
  s.require_charge 0              # Lambda is neutral
end

# Kinematic fit: constrain the Lambda mass from p pi-
tag.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag.with_decay_card(decay_card).apply
tag
  .note(:tag_mode_unavailable, "Paper uses only 2 of the available Lambda_c tag modes: anti-p K_S0 and anti-p K+ pi-.")
  .note(:st_selection, "Charged tracks: |cos(theta)|<0.93, Vz<10cm, Vr<1cm. PID: L(p)>L(K) and L(p)>L(pi) for protons; L(K)>L(pi) for kaons. No PID on pions from Lambda/KS0 decays. DeltaE cuts: mode-dependent at +/-2.5 sigma resolution. M_BC signal extraction via ML fit with ARGUS background.")
  .note(:ks_lambda_reconstruction, "K_S0 -> pi+ pi- and Lambda -> p pi-: vertex fit chi2<100, secondary vertex fit with momentum pointing back to IP constraint. Decay vertex > 2*resolution from IP. K_S0 mass [487,511] MeV/c^2; Lambda mass [1111,1121] MeV/c^2. Only one pair of tracks kept per event.")
  .note(:inclusive_signal, "Inclusive Lambda_c+ -> Lambda + X: signal yield from 2D distribution of M_BC vs M(p pi-). 5 sideband regions (A,B,C,D,E) used for background subtraction. Net signal yield computed in 5x4 (p, |cos(theta)|) Lambda kinematic intervals.")
  .note(:lambda_efficiency, "Lambda detection efficiency determined from control samples J/psi -> Lambda anti-Lambda and J/psi -> anti-p K+ Lambda using tag-and-probe technique. Efficiencies weighted by Lambda momentum and polar angle distributions in DT signal. BF calculation: B_sig = sum_j(N_sig_j / eps_sig_j) / sum_i(N_tag_i).")
  .note(:branching_fraction, "B(Lambda_c+ -> Lambda + X) = (38.2 +2.8/-2.2 +/- 0.8)%. Systematic 2.3% total (control sample stats 0.6%, Lambda efficiency bias 1.1%, tag efficiency bias 1.6%, interval choice 0.5%, tag yields 0.9%, background subtraction 0.3%).")
  .note(:cp_violation, "A_CP = (2.1 +7.0/-6.6 +/- 1.4)%. No evidence of CP violation. Separate BFs: B(Lambda_c+ -> Lambda+X) = (39.4 +4.7/-3.4)%, B(anti-Lambda_c- -> anti-Lambda+X) = (37.8 +3.8/-2.9)%.")
tag.execute_on([data_4600, incMC_4600, exMC])