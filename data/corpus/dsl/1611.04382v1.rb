# BESIII: absolute branching fraction for Lambda_c+ -> Lambda mu+ nu_mu
# (567 pb^-1 at sqrt(s) = 4.6 GeV)
#
# Single-tag (ST) / double-tag (DT) technique at the Lambda_c+ Lambda_c- mass
# threshold: the anti-Lambda_c- is reconstructed in eleven hadronic modes (ST
# sample) and the semileptonic decay Lambda_c+ -> Lambda mu+ nu_mu is searched
# for in the system recoiling against it (DT sample).  The neutrino is not
# detected; its kinematics is inferred from four-momentum conservation through
# U_miss = E_miss - |p_miss| c, which peaks at zero for signal.  The absolute
# branching fraction follows from N_obs / (N_tot x epsilon x B(Lambda -> p pi-)),
# so the ST-side uncertainties cancel.
#
# The tag modes and their selection follow the companion measurement
# (Ref. [14] of the paper); the signal side contributes a Lambda -> p pi- pair
# plus a muon.

### Dataset description ###
data_4600  = DatasetManager.real_data.find("703_4600")    # 567 pb^-1 at sqrt(s) = 4.6 GeV
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600") # generic Lc+ Lc-, D(s)*, ISR psi, QED

### Decay card for the exclusive signal MC ###
# Lambda_c+ decays only to Lambda mu+ nu_mu (form-factor model of Ref. [10];
# phase space used here since the model shape is not expressible) while the
# anti-Lambda_c- tag decays inclusively.
decay_card_lmu = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-                        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Lambda0  mu+  nu_mu                              PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                                         PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_lmu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "LambdacToLambdaMuNu"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_lmu
  config.cross_section   = :default
end
exMC_lmu.save_to_config(format: :yaml, file_path: 'temp_for_test')

datasets = [data_4600, incMC_4600]

### Event selection (BOSS) — ST anti-Lambda_c- + Lambda_c+ -> Lambda mu+ nu_mu ###
alg_name = "LcTagLambdaMuNu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.600] })
   .set_alias({ "std::vector<double>" => "Vdouble" })
   .note(:tag_side_selection,
         "The anti-Lambda_c- is reconstructed in eleven hadronic modes with the same " \
         "charged-track, photon, pi0, K_S0 and Lambda selection criteria as the " \
         "companion single-tag measurement (Ref. [14] of the paper): |cos(theta)| < 0.93 " \
         "for charged tracks, good tracks (except K_S0 daughters) within " \
         "V_xy < 1.0 cm and V_z < 10.0 cm, pion/kaon separation from the dE/dx and TOF " \
         "confidence levels, EMC photon energy thresholds 25 MeV barrel / 50 MeV endcap, " \
         "0.110 < M(gamma gamma) < 0.155 GeV/c^2 for pi0 with a 1C mass constraint, " \
         "|M(pi+ pi-) - M(K_S0)^PDG| window with a positive decay length for K_S0 and " \
         "M(p pi-) within the Lambda window for Lambda.")
   .note(:tag_mbc_deltae,
         "Tag candidates are identified with M_BC = sqrt(E_beam^2/c^4 - |p_Lambdabar_c-|^2/c^2) " \
         "and Delta E = E_beam - E_Lambdabar_c-, each mode keeping the approximately " \
         "3 sigma_DeltaE window around the Delta E peak (mode-dependent bounds spanning " \
         "-0.025..0.062 GeV).  The per-mode yields are obtained from a fit to the M_BC " \
         "distributions in the signal region (2.280, 2.296) GeV/c^2, giving a total ST " \
         "yield N_tot = 14415 +- 159 summed over the eleven modes.  Stored and applied at " \
         "the ROOT stage.")
   .note(:signal_lambda_selection,
         "The signal Lambda is formed from a p pi- combination constrained by a common " \
         "vertex fit to have a positive decay length L; when several candidates are " \
         "formed the one with the largest L/sigma_L is retained, where sigma_L is the " \
         "resolution of the measured L.  The Lambda and mu+ are taken from the tracks " \
         "recoiling against the ST anti-Lambda_c-.")
   .note(:signal_muon_pid,
         "Muon PID uses probabilities from the MDC dE/dx, the TOF and the EMC: a mu+ " \
         "candidate must satisfy L'_mu > 0.001, L'_mu > L'_e and L'_mu > L'_K, where " \
         "L'_mu, L'_e and L'_K are the muon, electron and kaon probabilities.")
   .note(:background_veto,
         "Backgrounds are dominated by Lambda_c+ -> Lambda pi+, Sigma0 pi+ and " \
         "Lambda pi+ pi0.  Lambda_c+ -> Lambda pi+ and Lambda_c+ -> Sigma0 pi+ are " \
         "rejected by requiring M(Lambda mu+) < 2.12 GeV/c^2.  The " \
         "Lambda_c+ -> Lambda pi+ pi0 background is suppressed by requiring the largest " \
         "energy of any unused photon, E_gamma_max, below 0.25 GeV and the muon " \
         "candidate's EMC deposited energy below 0.30 GeV.  The surviving " \
         "Lambda_c+ -> Lambda pi+ pi0 peaking background is fixed in the fit to " \
         "N_bkg = 37.1 +- 2.3, estimated as " \
         "N_bkg = N_tot x B(Lambda_c+ -> Lambda pi+ pi0) x eta with " \
         "eta = (3.67 +- 0.05)% including B(Lambda -> p pi-) and B(pi0 -> gamma gamma).")
   .note(:missing_neutrino_observable,
         "The neutrino is undetected: U_miss = E_miss - |p_miss| c is formed from " \
         "E_miss = E_beam - E_Lambda - E_mu+ and p_miss = p_Lambda_c+ - p_Lambda - p_mu+, " \
         "with p_Lambda_c+ = -p_hat_tag sqrt(E_beam^2/c^2 - m^2_Lambdabar_c-) using the " \
         "ST Lambda_c- direction and the nominal anti-Lambda_c- mass.  The signal peaks " \
         "at zero.  The yield is extracted by fitting the U_miss distribution with a " \
         "Gaussian core plus two power-law tails for ISR/FSR effects (shape parameters " \
         "fixed from signal MC), a double Gaussian for the Lambda_c+ -> Lambda pi+ pi0 " \
         "peaking background and an MC-derived shape for the other combinatorial " \
         "backgrounds, giving N_obs = 78.7 +- 10.5.  Handled in the ROOT analysis.")
   .note(:signal_mc_model,
         "The Lambda_c+ -> Lambda mu+ nu_mu signal is generated with the form factor " \
         "obtained using Heavy Quark Effective Theory and QCD sum rules; varying the " \
         "form-factor parameterisation and accounting for the q^2 dependence observed in " \
         "data gives a 5.2% model uncertainty.")
   .note(:absolute_branching_fraction,
         "B(Lambda_c+ -> Lambda mu+ nu_mu) = N_obs / (N_tot x epsilon x B(Lambda -> p pi-)) " \
         "with epsilon = (24.5 +- 0.2)% (DT/ST efficiency ratio weighted by the ST yields, " \
         "excluding B(Lambda -> p pi-)); result (3.49 +- 0.46 +- 0.27)%, giving " \
         "B(Lambda_c+ -> Lambda mu+ nu_mu)/B(Lambda_c+ -> Lambda e+ nu_e) = 0.96 +- 0.16 " \
         "+- 0.04.")

# Tag side: the eleven hadronic anti-Lambda_c- modes, as in the companion
# single-tag measurement.  charm -1 pins the tagged anti-Lambda_c-.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,           # p K_S0
          :LambdacPtoKPiP,          # p K- pi+
          :LambdacPtoKsPi0P,        # p K_S0 pi0
          :LambdacPtoKsPiPiP,       # p K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,       # p K- pi+ pi0
          :LambdacPtoPiPiP,         # p pi+ pi-
          :LambdacPtoLambdaPi,      # Lambda pi+
          :LambdacPtoLambdaPiPi0,   # Lambda pi+ pi0
          :LambdacPtoLambdaPiPiPi,  # Lambda pi+ pi+ pi-
          :LambdacPtoLambdaPiOmega, # Lambda pi+ omega
          :LambdacPtoPiSIGMA0LambdaGam # Sigma0 pi+ (Sigma0 -> gamma Lambda)
  t.charm -1
end

# Signal side: the Lambda -> p pi- pair (recoiling against the tag) plus the
# muon, and the undetected neutrino.
alg.signal_side do |s|
  s.photons 0
  s.charged(prp: 1, pim: 1, mup: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.missing :nu_mu
end

# Kinematic fit: 4-momentum conservation including the missing neutrino, with
# the Lambda mass constraint over its two daughters.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.chi2_cut 200
end

alg.with_decay_card(decay_card_lmu).apply
alg.execute_on(datasets + [exMC_lmu])
