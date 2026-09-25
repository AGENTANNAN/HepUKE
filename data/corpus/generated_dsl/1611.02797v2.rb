# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# 4.599 GeV data (BOSS 703); 567 pb^-1 -> sample name follows [BOSS]_[CMS_MeV]
data_4600  = DatasetManager.real_data.find("703_4600")      # Real data at 4.599 GeV
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")   # Corresponding inclusive MC

# Decay card for the signal process: Lambda_c+ -> n K_S0 pi+, anti-Lambda_c- generic, K_S0 -> pi+ pi-
# (no intermediate resonance at the ee -> Lambda_c anti-Lambda_c threshold: use psi(4260) as KKMC top mother)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 n0 K_S0 pi+ PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC (500k events), 4.599 GeV
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_Lambdac_to_nKsPi"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — tag-based analysis ###
alg_name = "LambdacToNKsPi"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.599]})          # sqrt(s) = 4.599 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})
   .with_decay_card(decay_card_signal)

# --- Tag side: anti-Lambda_c- reconstructed in eleven hadronic modes (single tag) ---
alg.tag_side(:Lambdac) do |t|
  # full DTagAlg channel-name symbols (verify against config/tag_modes.yaml)
  t.modes :LambdacPtoKsP,
          :LambdacPtoKsPPi0,
          :LambdacPtoKsPPiPi,
          :LambdacPtoKPiP,
          :LambdacPtoKPiPPi0,
          :LambdacPtoKPiPPiPi,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigma0Pi,
          :LambdacPtoSigmaPPiPi
  t.charm -1                                           # pin the tagged anti-Lambda_c- (charm -1)
  # explicit tag-side windows requested by the description (store-not-cut bypass)
  t.window :mBC, min: 2.280, max: 2.296               # M_BC in (2.280, 2.296) GeV/c^2
  t.window :deltaE, abs: 0.02                          # ~3 sigma DeltaE window
end

# --- Signal side: Lambda_c+ -> n K_S0 pi+ (K_S0 -> pi+ pi-, missing neutron) ---
alg.signal_side do |s|
  s.photons 0                                          # no photons on the signal side
  s.charged(pip: 2, pim: 1)                            # K_S0 -> pi+ pi- (two tracks) + the direct signal pi+
  s.require_charge 1                                   # net signal-side charge +1
  s.missing :n                                         # undetected neutron (massive missing particle)
end

# --- Kinematic fit: 4-momentum conservation incl. the missing neutron + K_S0 mass constraint ---
alg.fit do |f|
  f.constrain_four_momentum                            # tag + signal + missing neutron = measured CMS 4-vector
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)  # K_S0 nominal-mass constraint
  f.chi2_cut 200                                       # loose chi2 cut (tight cut applied in ROOT)
end

# --- Inexpressible tag-side procedures (handled by DTagAlg) captured as notes ---
alg.note(:tag_side_selection,
         "DTagAlg tag-side reconstruction: charged tracks |cos(theta)|<0.93, |Vz|<10 cm, " \
         "|Vr|<1 cm (V0 daughters exempt); PID from dE/dx, TOF and EMC requiring pi (L_pi>L_K), " \
         "K (L_K>L_pi) and proton (L'_p>L'_pi and L'_p>L'_K); photons E>25 MeV barrel " \
         "(|cos(theta)|<=0.80) or E>50 MeV endcap (0.86<=|cos(theta)|<=0.92), >10 deg from the " \
         "nearest charged track, EMC time in (0,700) ns; pi0 from gamma gamma with invariant mass " \
         "in (0.110,0.155) GeV/c^2 using a 1C mass constraint with chi2<20; K_S0 -> pi+ pi- and " \
         "Lambda -> p-bar pi+ built from vertex fits with flight length L>0 and masses " \
         "M(pi+pi-) in (0.485,0.510) GeV/c^2, M(p-bar pi+) in (1.110,1.121) GeV/c^2.")
   .note(:background_veto,
         "Tag-side background vetoes: Sigma0 suppressed via M(gamma Lambda-bar) in (1.179,1.205) " \
         "GeV/c^2; Sigma- suppressed via M(p-bar pi0) in (1.173,1.200) GeV/c^2; additional " \
         "Lambda/Sigma- and K_S0/Lambda vetoes applied for specific tag modes.")

# Render the tag spec (no Selection argument) and run on data + MC
alg.apply
root_files = alg.execute_on([data_4600, incMC_4600, exMC_signal])