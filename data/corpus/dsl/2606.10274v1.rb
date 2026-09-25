### Dataset description ###
# Energy-scan samples 4.128 - 4.226 GeV (total 7.33 fb^-1), used for Ds*+ Ds- production
data_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230")
]

inc_mc_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Signal decay card: e+e- -> Ds*+ Ds- ; Ds*+ -> e+ e- Ds+  (EM Dalitz decay)
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0000  D_s*+ D_s-                                 PHSP;
    Enddecay

    Decay D_s*+
    1.0000  e+ e- D_s+                                 PHSP;
    Enddecay

    Decay D_s-
    1.0000  K+ K- pi-                                  PHSP;
    Enddecay

    Decay D_s+
    1.0000  K+ K- pi+                                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsStar_ee_Ds_signal"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

########################################################################
# Tag-based analysis: Ds tag on one side, signal side = e+e- + missing Ds
# Eleven Ds tag modes:
#   K_S0 K+/-, K+K- pi+/-, K_S0 K+/- pi0, K+K- pi+/- pi0,
#   K_S0 K-/+ pi+/- pi+/-, pi- pi+ pi+/-, pi+/- eta, pi+/- pi0 eta,
#   pi+/- eta'(pi+pi-eta), pi+/- eta'(gamma rho0), K+/- pi+ pi-
########################################################################
alg = TagAnalysis.new("DsStarEEDsDalitz")
alg.set_header(["DsStarEEDsDalitzAlg/DsStarEEDsDalitz.h"])
   .set_constant({"ECMS" => [:double, 4.178]})     # nominal; per-run measured Ecms used at runtime
   .with_decay_card(decay_card_signal)

# ---- Ds single tag ----
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,
          :DstoKKPi,
          :DstoKsKPi0,
          :DstoKKPiPi0,
          :DstoKsKPiPi,
          :DstoPiPiPi,
          :DstoPiEta,
          :DstoPiPi0Eta,
          :DstoPiEPPiPiEta,
          :DstoPiEtaPrimeGammaRho0,
          :DstoKPiPi
end

# ---- Signal side: an extra e+e- pair and the missing bachelor/daughter Ds ----
alg.signal_side do |s|
  s.charged(ep: 1, em: 1)         # reconstruct the additional e+ e- pair (electron ID)
  s.require_charge 0
  s.min_photon_angle 10.0
  s.missing :"D_s+"               # missing Ds recoiling against tag + e+e-
end

# ---- Kinematic fit ----
alg.fit do |f|
  f.constrain_four_momentum       # tag + e+ e- + missing Ds = ecms_lab (measured boost)
  f.chi2_cut 200                  # loose cut; tight cuts applied in ROOT
  f.store_fitted_momenta
end

# Capture BOSS-side procedures / cuts not directly expressible in DSL primitives.
alg.note(:signal_electron_momentum,
         "additional e+ and e- candidates are required to have momentum < 200 MeV/c due " \
         "to the small phase space of Ds*+ -> e+ e- Ds+ decay.")
    .note(:signal_electron_pid,
         "electron PID uses dE/dx-only likelihoods: L(e) > 0, L(e) > L(K) and L(e) > L(pi).")
    .note(:gamma_conversion_veto,
         "gamma-conversion background suppressed by requiring the reconstructed conversion " \
         "vertex Rxy of the e+e- pair to be < 2.0 cm.")
    .note(:tag_side_selection,
         "tag Ds+/-: track |Vz| < 10 cm and Vxy < 1 cm (for non-KS tracks); |cos(theta)| < 0.93; " \
         "K/pi PID by L(K) vs L(pi) from MDC dE/dx + TOF; KS -> pi+pi- via secondary vertex with " \
         "|Vz| < 20 cm, M(pi+pi-) in [0.487, 0.511] GeV/c^2 and decay length > 2 sigma; " \
         "pi0 -> gamma gamma with M(gg) in [0.115, 0.150] GeV/c^2 and 1C mass-constrained fit; " \
         "eta -> gamma gamma with M(gg) in [0.490, 0.580] GeV/c^2 and 1C fit; " \
         "rho0 -> pi+pi- with M(pi+pi-) in [0.570, 0.970] GeV/c^2; " \
         "eta'-> pi+pi-eta in [0.943, 0.973] GeV/c^2 and eta'-> gamma rho0 in [0.946, 0.970] GeV/c^2.")
    .note(:tag_ds_mass_window,
         "reconstructed Ds+/- invariant mass required within [1.85, 2.06] GeV/c^2; " \
         "additional pi+/- and pi0 in the tag are required to have momentum > 100 MeV/c to suppress " \
         "wrong combinations from Ds*+/-, D*0 and D*+/- decays.")
    .note(:tag_ks_veto,
         "for Ds -> K pi+ pi- and Ds -> pi+ pi- pi+/- modes, any pi+pi- combination with mass " \
         "in [0.468, 0.528] GeV/c^2 is vetoed to remove KS-associated background.")
    .note(:m_recoil_window,
         "recoil mass Mrec of tag Ds is required in the energy-dependent window given by " \
         "Table I of the paper (roughly [2.04, 2.22] GeV/c^2) and the joint e+e-Ds recoil " \
         "M_rec_{e+e-Ds+/-} is required in (1.93, 2.03) GeV/c^2.")
    .note(:multi_candidate_arbitration,
         "if multiple tag + e+e- combinations survive for a mode, the one with M_rec_{e+e-Ds+/-} " \
         "closest to the Ds+ nominal mass is retained.")

alg.apply
alg.execute_on(data_points + inc_mc_points + exMC_signal)
