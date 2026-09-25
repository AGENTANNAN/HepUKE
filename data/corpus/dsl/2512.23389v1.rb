### Dataset description ###
psipp_data   = DatasetManager.real_data.find("712_3773")     # 20.3 fb^-1 at sqrt(s) = 3.773 GeV (psi(3770))
psipp_incMC  = DatasetManager.inclusive_mc.find("712_3773")  # Corresponding inclusive MC

# Decay card for the signal double-tag process:
# psi(3770) -> D0 Dbar0, with D0 -> Ks0 pi0 eta on signal side and Dbar0 -> K+ pi- pi0 on tag side.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0                       PHSP;
    Enddecay

    Decay D0
    1.000  K_S0  pi0  eta                    PHSP;
    Enddecay

    Decay anti-D0
    1.000  K+   pi-  pi0                     PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                          PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                      PHSP;
    Enddecay

    Decay eta
    1.000  gamma  gamma                      PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0_to_KsPi0Eta_signal_MC"
  config.related_dataset = psipp_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — TagAnalysis (double tag) ###
alg_name = "D0KsPi0EtaAmp"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })     # sqrt(s) = 3.773 GeV
   .note(:pi0_mispairing_veto,
         "Reject events where any pair of the four photons (two from pi0 and two from eta) " \
         "falls in the pi0 mass window 0.115 < M_gg < 0.150 GeV/c^2 to suppress " \
         "mis-pairing of photons between the pi0 and the eta.")
   .note(:kinematic_fit_details,
         "Four-constraint (4C) kinematic fit constrains the invariant masses of (gamma,gamma)_eta, " \
         "(gamma,gamma)_pi0, (pi+,pi-)_Ks0 and the D0 candidate to their nominal PDG masses.")
   .note(:best_candidate,
         "In case of multiple candidates per event, the one with the minimum |DeltaE| is chosen.")

# --- Tag side: Dbar0 reconstructed in three hadronic tag modes ---
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                                          # anti-D0 tag
  t.window :mBC, min: 1.859, max: 1.873               # 1.859 < M_BC(tag) < 1.873 GeV/c^2
end

# --- Signal side: D0 -> Ks0 pi0 eta reconstructed as second tag ---
alg.tag_side(:D0) do |t|
  t.modes :D0toKsPi0Eta
  t.charm  1                                          # D0
  t.rank_by :inv
  t.window :mBC,    min: 1.858, max: 1.871            # 1.858 < M_BC(sig) < 1.871 GeV/c^2
  t.window :deltaE, min: -0.050, max: 0.045           # -50 < DeltaE(sig) < 45 MeV
end

alg.signal_side do |s|
  s.photons 0                                         # all showers used by tag2 (pi0 + eta)
end

# --- Kinematic fit: 4C total four-momentum plus intermediate mass constraints ---
alg.fit do |f|
  f.constrain_four_momentum                                                            # constrain to CMS 4-vector
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)                         # M(Dbar0 tag) = M_D0
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)                         # M(signal D0) = M_D0
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg.apply
root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal])
