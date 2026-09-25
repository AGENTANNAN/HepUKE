### dataset description ###
# psi(3770) data at sqrt(s) = 3.773 GeV, 20.3 fb^-1
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay card: psi(3770) -> D+ D-, with signal side D+ -> pi+ eta eta
# (nominal amplitude model: D+ -> a0(980)+ eta, a0(980)+ -> pi+ eta) and
# the D- tag reconstructed in one of six hadronic modes.
decay_card_signal = <<~DECAYCARD
    Alias a0plus a_0+

    Decay psi(3770)
    1.000  D+  D-                          PHSP;
    Enddecay

    Decay D+
    1.000  a0plus  eta                     PHSP;
    Enddecay

    Decay a0plus
    1.000  pi+  eta                        PHSP;
    Enddecay

    Decay D-
    0.1667  K+  pi-  pi-                   PHSP;
    0.1667  K+  pi-  pi-  pi0              PHSP;
    0.1667  K_S0  pi-                      PHSP;
    0.1667  K_S0  pi-  pi0                 PHSP;
    0.1667  K_S0  pi-  pi-  pi+            PHSP;
    0.1667  K+   K-   pi-                  PHSP;
    Enddecay

    Decay eta
    1.000  gamma  gamma                    PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                    PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                        PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_Dp_pi_eta_eta_Dm_tag"
  config.related_dataset = data_3773
  config.events = 200000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### event selection (BOSS) --- Double-tag D+ -> pi+ eta eta ###
alg = TagAnalysis.new("DpToPiEtaEta")
alg.set_header(["DpToPiEtaEtaAlg/DpToPiEtaEta.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# Tag side: D- reconstructed in six hadronic modes
alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi,
          :DptoKPiPiPi0,
          :DptoKsPi,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm(-1)
end

# Signal side: D+ -> pi+ eta eta, eta -> gamma gamma (four photons + one pi+)
alg.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge(+1)
  s.photons 4
  s.min_photon_angle 10.0
end

# Kinematic fit: 4-momentum conservation plus mass constraints on the two eta
# candidates (gamma gamma pairs) and on the signal D+.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg
  .note(:tag_selection_reference, "tracking, PID, K_S0, pi0 and eta reconstruction follow the standard BESIII D-tag prescription; the eta signal window is widened to 0.45 < M(gamma gamma)_eta < 0.65 GeV/c^2 to increase efficiency")
  .note(:tag_mBC_deltaE_windows, "for both the tag and signal sides, candidates with M_BC < 1.83 GeV/c^2 or |Delta E| > 0.1 GeV are rejected before the fit; per-mode Delta E windows use 3.5 sigma of the corresponding resolution")
  .note(:signal_mBC_window, "signal-side M_BC window [1.860, 1.880] GeV/c^2 and |Delta E| < 0.040 GeV applied before the kinematic fit")
  .note(:best_candidate_choice, "one candidate per tag mode is kept with |Delta E| closest to zero; the signal candidate with M_BC closest to the nominal D+ mass is retained")
  .note(:bdtg_background_suppression, "a Gradient Boosted Decision Tree (BDTG) trained on inclusive MC is applied to suppress no-eta backgrounds. Inputs: ln(chi2) of the double-eta mass-constraint fit, and the invariant mass and helicity-angle cosine of the higher-energy photon in each eta candidate. Working point retains 78.3% signal and rejects 89.7% background; Dalitz-plot distortion is negligible")
  .apply

root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])
