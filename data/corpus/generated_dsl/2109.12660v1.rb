# ================================================================
# Dataset preparation — six energy points spanning 4.178–4.226 GeV
# (sample name convention: [BOSS version]_[CMS energy in MeV])
# ================================================================
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

data_points   = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
inc_mc_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card for the signal process e+e- -> D_s*+ D_s-, with the signal-side chain
# D_s*+ -> gamma D_s+, D_s+ -> pi+ pi0 pi0, pi0 -> gamma gamma.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0  D_s*+  D_s-    PHSP;
    Enddecay

    Decay D_s*+
    1.0  gamma  D_s+    PHSP;
    Enddecay

    Decay D_s+
    1.0  pi+  pi0  pi0  PHSP;
    Enddecay

    Decay pi0
    1.0  gamma  gamma   PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC generated once per energy point (√s scan)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsst_ds_gamma_pip_pi0pi0"  # auto-suffixed per energy point
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ================================================================
# Event selection — tag-based analysis (TagAnalysis)
# ================================================================
alg_name = "DsstDsPiPi0Pi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.226] })   # nominal value; per-run beam energy taken from MeasuredEcmsSvc across the scan
   .with_decay_card(decay_card_signal)

# ---- Tag side: hadronic D_s- in seven modes ----
alg.tag_side(:Ds) do |t|
  t.modes :DstoKSK,          # K_S0 K-
          :DstoKKPi,         # K+ K- pi-
          :DstoKSKPi0,       # K_S0 K+ pi0
          :DstoKSKPiPi,      # K_S0 K- pi- pi+
          :DstoKSKPimPim,    # K_S0 K+ pi- pi-
          :DstoPiEtaPrime,   # pi- eta' (eta' -> pi+ pi- eta)
          :DstoKPiPi         # K- pi+ pi-
  t.charm -1                 # pin the tagged D_s-
end

# ---- Signal side: five photons + one pi+ ----
alg.signal_side do |s|
  s.photons 5                # four from the two pi0 + the D_s* transition photon
  s.charged(pip: 1)          # the single pi+ from D_s+
  s.min_photon_angle 10.0    # photon angle > 10 degrees
  s.min_photon_energy 0.025  # barrel energy floor 25 MeV (endcap 50 MeV inside EMC)
end

# ---- Kinematic fit: 4C + pi0 mass constraint ----
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 -> gamma gamma mass constraint
  f.chi2_cut 200            # loose nominal chi2 (final optimisation in ROOT)
end

# ---- BOSS-side criteria that the current DSL surface cannot express formally ----
alg.note(:tag_mass_window,
         "per-tag-mode M(D_s-) windows ~[1.94, 1.996] GeV/c^2; tag mass is stored "
         "unconditionally and windowed in the ROOT analysis (store-not-cut)")
alg.note(:recoil_mass_window,
         "recoil-mass window against the tagged D_s- (U_miss / M_recoil) applied in ROOT")
alg.note(:background_veto,
         "candidate ranking by minimum chi2 of the 8C fit (4C + tag D_s and D_s* mass "
         "constraints) and the resonance windows pi0 [0.115,0.150], eta [0.490,0.580], "
         "eta' [0.946,0.970] GeV/c^2 are not fully expressible with the v1 photon arity cap")

alg.apply
root_files = alg.execute_on(data_points + inc_mc_points + exMCs_signal)