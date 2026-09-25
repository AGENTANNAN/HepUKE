### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data at 3.773 GeV (20.3 fb^-1)
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")   # Inclusive MC at psi(3770)

# Decay card for the signal process: D+ -> pi+ pi0 pi0 (with tag D- decays)
signal_decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.000 pi+ pi0 pi0                        PHSP;
  Enddecay

  Decay D-
  0.1667 K+ pi- pi-                        PHSP;
  0.1667 K+ pi- pi- pi0                    PHSP;
  0.1667 K_S0 pi-                          PHSP;
  0.1667 K_S0 pi- pi0                      PHSP;
  0.1667 K_S0 pi- pi- pi+                  PHSP;
  0.1667 K+ K- pi-                         PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi-                            PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Signal MC sample (used to compute detection efficiency; generated per amplitude-analysis results)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Dp_pippi0pi0_signal"
  config.related_dataset = psipp_data
  config.events         = 2_000_000
  config.decay_card     = signal_decay_card
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS, tag-based) ###
alg_name = "DpToPipPi0Pi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

# ---- Tag side: D- reconstructed through six hadronic tag modes ----
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKPiPiPi0,
          :DptoKsPi,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

# ---- Signal side: D+ -> pi+ pi0 pi0 ----
alg.signal_side do |s|
  s.charged(pip: 1)
  s.photons 4..48
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.require_charge 1
end

# ---- Kinematic fit ----
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:ks0_pi0pi0_veto,
         "K_S0 veto on the signal side of the amplitude analysis: reject events with " \
         "0.428 < M(pi0 pi0) < 0.548 GeV/c^2 to suppress D+ -> K_S0 pi+ (K_S0 -> pi0 pi0).")
   .note(:mispartition_veto_DDbar,
         "Mispartition veto for D+ -> pi+ pi0 pi0 vs D- -> K+ K- pi-: reject events that " \
         "simultaneously satisfy |M(K- pi+ pi0) - 1.865| < 0.05 GeV/c^2 and " \
         "|M(K+ pi- pi0) - 1.865| < 0.05 GeV/c^2.")
   .note(:pi0_1c_kalman_fit,
         "Each pi0 is reconstructed from a gamma-gamma pair with 0.115 < M(gg) < 0.150 GeV/c^2, " \
         "at least one photon in the barrel EMC, and a 1-C kinematic fit constraining M(gg) to the " \
         "nominal pi0 mass with chi2 < 50.")
   .note(:seven_c_kinematic_fit,
         "Amplitude analysis uses a seven-constraint kinematic fit: four-momenta of the " \
         "final-state particles constrained to the initial e+e- 4-momentum, plus D+ mass constraint " \
         "and two pi0 mass constraints.")
   .note(:best_candidate_selection,
         "When multiple signal-side candidates coexist in an event, keep the candidate that " \
         "minimises dE_tag^2 + dE_sig^2.")
   .note(:mbc_signal_region,
         "Signal region 1.865 < M_BC^tag < 1.875 GeV/c^2 and 1.865 < M_BC^sig < 1.875 GeV/c^2.")
   .note(:delta_e_windows,
         "Per-mode dE windows: D+->pi+pi0pi0 (-0.100,0.045); D-->K+pi-pi- (-0.025,0.024); " \
         "D-->K+pi-pi-pi0 (-0.057,0.046); D-->KS0pi- (-0.025,0.026); D-->KS0pi-pi0 (-0.062,0.049); " \
         "D-->KS0pi-pi-pi+ (-0.028,0.027); D-->K+K-pi- (-0.024,0.023) GeV.")

alg.apply
root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal])
