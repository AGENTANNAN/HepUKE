# D0 -> K+ K- pi0 pi0 amplitude analysis and branching fraction (BESIII, psi(3770), 20.3 fb^-1)
# Double-tag technique: signal D0 -> K+ K- pi0 pi0 vs Dbar0 -> K+pi-, K+pi-pi0, K+pi-pi-pi+

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal D0 -> K+ K- pi0 pi0 with a K*(892)+ K*(892)- dominant intermediate
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D0    anti-D0                 VSS_MIX;
    Enddecay

    Decay D0
    1.0000  K+  K-  pi0  pi0              PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                       PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0_KKpi0pi0_signal"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — Double-tag D-tag analysis ###
alg_name = "D0toKKpi0pi0DT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

# Tag side: anti-D0 reconstructed in K+pi-, K+pi-pi0, K+pi-pi-pi+
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0 -> K+ K- pi0 pi0 (2 kaons of opposite sign + 2 pi0 -> 4 photons)
alg.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.require_charge 0
  s.photons 4
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: full 4-momentum constraint + two pi0 mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures captured as notes
alg.note(:pi0_mass_window,
         "pi0 candidates require |M(gg) - m_pi0| in [0.115, 0.150] GeV/c^2 with at least one photon in the barrel EMC; a 1-C Kalman fit to the nominal pi0 mass with chi2 < 50 is applied.")
   .note(:ks0_veto,
         "Ks0 veto on the signal side: M(pi0 pi0) not in [0.460, 0.525] GeV/c^2 to suppress D0 -> Ks0 K+ K-.")
   .note(:opening_angle_DDbar,
         "Opening angle between the reconstructed D0 and anti-D0 candidates required to be > 167 degrees to suppress non-DDbar backgrounds.")
   .note(:tag_deltaE_windows,
         "Tag-side DeltaE windows: Dbar0 -> K+pi-: (-0.025,0.025) GeV; Dbar0 -> K+pi-pi0: (-0.055,0.040) GeV; Dbar0 -> K+pi-pi-pi+: (-0.025,0.025) GeV. Applied in ROOT after storage.")
   .note(:signal_deltaE_window,
         "Signal-side DeltaE window (-0.020, 0.020) GeV applied in ROOT.")
   .note(:best_candidate_selection,
         "Best candidate per tag mode chosen by minimum |DeltaE_tag|; best signal candidate by minimum |DeltaE_sig|.")
   .note(:amplitude_analysis_region,
         "Amplitude analysis signal region 1.861 < MBC < 1.870 GeV/c^2 for both signal and tag; applied in ROOT.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])
