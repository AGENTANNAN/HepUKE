# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Signal process: e+e- -> Ds*+ Ds-, with the electromagnetic decay Ds*+ -> e+e- Ds+.
# Data taken at the 4.130-4.230 GeV energy-scan points (eight points, ~7.33 fb^-1 in total).
data_4130 = DatasetManager.real_data.find("705_4130")
data_4160 = DatasetManager.real_data.find("705_4160")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
scan_data_points = [data_4130, data_4160, data_4180, data_4190,
                    data_4200, data_4210, data_4220, data_4230]

# Matching inclusive MC for each energy point
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
scan_incMC = [incMC_4130, incMC_4160, incMC_4180, incMC_4190,
              incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card for the signal process (EvtGen format, EvtGen particle names).
# Ds*+ -> e+e- Ds+ is an electromagnetic transition; PHSP is used for the 3-body final state.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s-    PHSP;
    Enddecay

    Decay D_s*+
    1.0000 e+ e- D_s+    PHSP;
    Enddecay

    End
DECAYCARD

# 200k exclusive-MC events for the signal mode, one sample per scan energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name   = "sig_dsstar_to_ee_ds"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS part, tag-based) ###
alg_name = "DsTagDsStarEE"
alg = TagAnalysis.new(alg_name)                       # tag layer: TagAnalysis < Algorithm
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.180]})        # representative CM energy of the scan (per-run beam energy comes from MeasuredEcmsSvc)
   .set_alias({"std::vector<double>" => "Vdouble"})
   .with_decay_card(decay_card_signal)

# --- Tag side: one Ds reconstructed in eleven tag modes ---
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,                     # K_S K
          :DstoKKPi,                    # K K pi
          :DstoKsKPi0,                  # K_S K pi0
          :DstoKKPiPi0,                 # K K pi pi0
          :DstoKsKPiPi,                 # K_S K pi pi
          :DstoPiPiPi,                  # pi pi pi
          :DstoPiEta,                   # pi eta
          :DstoPiPi0Eta,                # pi pi0 eta
          :DstoPiEtaPrimePiPiEta,       # pi eta'(-> pi+ pi- eta)
          :DstoPiEtaPrimeGammaRho,      # pi eta'(-> gamma rho0)
          :DstoKPiPi                    # K pi pi
  # Tag Ds mass window (store-not-cut opt-in; Ds mass handled as the beam-constrained mass)
  t.window :mBC, min: 1.85, max: 2.06
end

# --- Signal side: e+e- from Ds*+ -> e+e- Ds+, with the Ds+ left missing ---
alg.signal_side do |s|
  s.charged(ep: 1, em: 1)     # exactly one e+ and one e- on the signal side
  s.require_charge 0          # net charge zero
  s.min_photon_angle 10.0     # minimum photon angle to charged tracks (degrees)
  s.missing :Ds               # one missing Ds (massive form)
end

# --- 4C kinematic fit: tag + e+e- pair + missing Ds constrained to the measured CMS ---
alg.fit do |f|
  f.constrain_four_momentum   # tag + signal + missing = measured CMS four-momentum
  f.chi2_cut 200              # loose cut; tighter cuts applied downstream
  f.store_fitted_momenta      # store the fitted four-momenta per participant
end

# --- BOSS-side procedures that cannot be expressed in the tag DSL ---
alg.note(:tag_side_reconstruction,
    "tag tracks require |Vz| < 10 cm, Vxy < 1 cm and |cos(theta)| < 0.93, with K/pi " \
    "separation from MDC dE/dx plus TOF; tag-internal intermediate states are reconstructed " \
    "as K_S -> pi+pi- ([0.487, 0.511] GeV/c^2, decay length > 2 sigma, |Vz| < 20 cm), " \
    "pi0 -> gamma gamma ([0.115, 0.150] GeV/c^2, 1C mass-constrained fit), " \
    "eta -> gamma gamma ([0.490, 0.580] GeV/c^2, 1C mass-constrained fit), " \
    "rho0 -> pi+pi- ([0.570, 0.970] GeV/c^2), eta' -> pi+pi-eta ([0.943, 0.973] GeV/c^2) or " \
    "eta' -> gamma rho0 ([0.946, 0.970] GeV/c^2); extra tag pi/pi0 are required to have " \
    "momentum above 100 MeV/c, and a pi+pi- K_S veto in [0.468, 0.528] GeV/c^2 is applied " \
    "for the K pi pi and pi pi pi modes; all handled inside DTagAlg tag reconstruction")
   .note(:recoil_mass_windows,
    "the tag Ds recoil mass is required to lie in the energy-dependent window ~[2.04, 2.22] GeV/c^2 " \
    "and the e+e-Ds recoil mass in (1.93, 2.03) GeV/c^2; these are stored by the tag diagnostic " \
    "branches and windowed in the ROOT stage")
   .note(:signal_side_selection,
    "signal electrons/positrons are required to have momentum < 200 MeV/c (small " \
    "Ds*+ -> e+e- Ds+ phase space) and are identified by dE/dx-only likelihoods " \
    "(L(e) > 0, L(e) > L(K), L(e) > L(pi)); a gamma-conversion veto rejects e+e- pairs " \
    "whose vertex satisfies Rxy < 2.0 cm; the v1 tag DSL lepton PID thresholds are fixed " \
    "and not tunable, so this is applied at the ROOT stage")
   .note(:best_candidate_selection,
    "when several candidates survive per tag mode, the candidate whose e+e-Ds recoil mass is " \
    "closest to the Ds nominal mass is retained")

# Validate and render the tag specification (apply takes no Selection argument)
alg.apply

# Execute on real data, inclusive MC and the signal exclusive-MC samples
root_files = alg.execute_on(scan_data_points + scan_incMC + exMC_signal)