### Dataset description ###
# sqrt(s) = 4.600 GeV real data (~567/pb) and the matching inclusive MC (BOSS 703)
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for the signal process (EvtGen format).
# Signal: Lambda_c+ -> Lambda mu+ nu_mu in phase space, Lambda -> p pi-.
# The anti-Lambda_c- is left undecayed so EvtGen treats its decay inclusively.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 Lambda mu+ nu_mu PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Lc2LmuNu_signal_exclusive_mc"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — Lambda_c tag + semileptonic signal ###
alg_name = "Lc2LmuNuTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]})
   .with_decay_card(decay_card_signal)

# Tag side: fully reconstruct the anti-Lambda_c- in the eleven hadronic modes.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,          # Lambda_c+ -> p K_S0
          :LambdacPtoKPiP,         # Lambda_c+ -> p K- pi+
          :LambdacPtoKsPi0P,       # Lambda_c+ -> p K_S0 pi0
          :LambdacPtoKsPiPiP,      # Lambda_c+ -> p K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,      # Lambda_c+ -> p K- pi+ pi0
          :LambdacPtoPiPiP,        # Lambda_c+ -> p pi+ pi-
          :LambdacPtoLambdaPi,     # Lambda_c+ -> Lambda pi+
          :LambdacPtoLambdaPiPi0,  # Lambda_c+ -> Lambda pi+ pi0
          :LambdacPtoLambdaPiPiPi, # Lambda_c+ -> Lambda pi+ pi+ pi-
          :LambdacPtoLambdaPiOmega,# Lambda_c+ -> Lambda pi+ omega
          :LambdacPtoSigma0Pi      # Lambda_c+ -> Sigma0 pi+  (Sigma0 -> gamma Lambda)
  t.charm -1                       # pin the tagged side to the anti-Lambda_c-
  # Explicit tag-side windows requested by the analysis (store-not-cut otherwise)
  t.window :mBC,    min: 2.280, max: 2.296   # M_BC signal region (GeV/c^2)
  t.window :deltaE, min: -0.025, max: 0.062  # mode-dependent DeltaE envelope (GeV)
end

# Signal side: Lambda formed from p pi-, one mu+, no photons, net charge +1,
# and one undetected (massless) nu_mu.
alg.signal_side do |s|
  s.photons 0                              # no photons on the signal side
  s.charged(prp: 1, pim: 1, mup: 1)        # exactly one p, one pi-, one mu+
  s.require_charge 1                       # net charge +1
  s.min_photon_angle 10.0                  # shower-track isolation angle (deg)
  s.missing :nu_mu                         # undetected neutrino (massless form)
end

# Kinematic fit: 4-momentum conservation including the missing neutrino,
# Lambda mass constraint on the p pi- pair, and M(Lambda mu+) < 2.12 GeV/c^2
# veto against Lambda_c+ -> Lambda pi+ / Sigma0 pi+.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:prp, :pim, :mup).between(0.0, 2.12)
  f.chi2_cut 200
end

# BOSS-side procedures that the tag DSL cannot express directly.
alg.note(:tag_side_track_selection,
         "tag-side charged tracks required |cos(theta)| < 0.93, Vxy < 1 cm, |Vz| < 10 cm (applied inside DTagAlg reconstruction)")
   .note(:tag_side_pid,
         "tag-side pi/K separation from combined dE/dx and TOF confidence levels (DTagAlg internal tight-PID, type()==1)")
   .note(:tag_side_photon,
         "tag-side EMC photon energy thresholds 25 MeV (barrel) / 50 MeV (endcap)")
   .note(:tag_side_pi0,
         "tag-side pi0: mass window 0.110-0.155 GeV/c^2 with a 1C mass constraint (handled inside DTagAlg mode reconstruction)")
   .note(:tag_side_ks,
         "tag-side K_S0: mass window and positive decay-length requirement")
   .note(:tag_side_lambda,
         "tag-side Lambda: mass window (reconstructed from p pi-)")
   .note(:signal_lambda_vertex,
         "signal-side Lambda formed from p pi- with a common vertex fit requiring positive decay length; among passing candidates keep the one with the largest L/sigma_L")
   .note(:background_veto,
         "Lambda_c+ -> Lambda pi+ pi0 suppressed by requiring the largest unused-photon energy < 0.25 GeV and the muon EMC energy < 0.30 GeV; not expressible in the tag DSL")
   .note(:semileptonic_observables,
         "neutrino kinematics inferred from U_miss = E_miss - |p_miss|c (auto-stored as m_Umiss / m_Umiss2 alongside m_P4_miss_fit and m_q2)")

alg.apply
root_files = alg.execute_on([data_4600, incMC_4600, exMC_signal])