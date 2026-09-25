# D^0 -> K_S0 pi- l+ nu_l (l = e, mu) semileptonic form-factor analysis
# 20.3 fb^-1 psi(3770) data
# Tag-based analysis: anti-D0 tag reconstructed via 3 hadronic modes:
#   K+ pi-, K+ pi- pi- pi+, K+ pi- pi0
# Signal side: K_S0 (-> pi+ pi-) + pi- + l+ + missing nu_l
# Two separate SL channels -> two TagAnalysis objects (Rule T1).

### Datasets ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# --- Signal MC decay cards (e mode / mu mode) ---
decay_card_e = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0                             VSS_BMIX dm;
    Enddecay

    Decay D0
    1.000  anti-K0  pi-  e+  nu_e                  PHSP;
    Enddecay

    Decay anti-K0
    1.000  K_S0                                    PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                                PHSP;
    Enddecay

    Decay anti-D0
    0.5    K+  pi-                                 PHSP;
    0.5    K+  pi-  pi0                            PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mu = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0                             VSS_BMIX dm;
    Enddecay

    Decay D0
    1.000  anti-K0  pi-  mu+  nu_mu                PHSP;
    Enddecay

    Decay anti-K0
    1.000  K_S0                                    PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                                PHSP;
    Enddecay

    Decay anti-D0
    0.5    K+  pi-                                 PHSP;
    0.5    K+  pi-  pi0                            PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay

    End
DECAYCARD

exMC_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "D0toKsPiENu_signal_mc"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_e
  config.cross_section   = :default
end
exMC_e.save_to_config(format: :yaml, file_path: 'exMC_D0KsPiENu_config')

exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "D0toKsPiMuNu_signal_mc"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_mu
  config.cross_section   = :default
end
exMC_mu.save_to_config(format: :yaml, file_path: 'exMC_D0KsPiMuNu_config')

######################################################################
### Electron mode : D0 -> K_S0 pi- e+ nu_e                          ###
######################################################################
alg_e = TagAnalysis.new("D0KsPiENu")
alg_e.set_header(["D0KsPiENuAlg/D0KsPiENu.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .with_decay_card(decay_card_e)

# ST anti-D0 tag with three hadronic modes; charm = +1 for anti-D0
alg_e.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm 1
end

# Signal side: K_S0 -> pi+ pi- (2 tracks) + pi- + e+, plus missing nu_e (massless)
alg_e.signal_side do |s|
  s.charged(pip: 1, pim: 2, ep: 1)   # K_S0 daughters pi+ pi-, plus signal pi-, plus e+
  s.require_charge 0                 # +1 (pi+) -1 (pi-) -1 (pi-) +1 (e+) = 0
  s.missing :nu_e                    # massless neutrino
end

alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_e
  .note(:tag_deltaE_windows,
        "Mode-dependent tag deltaE windows: K+pi- [-0.027, 0.027] GeV, K+pi-pi+pi- [-0.026, 0.024] GeV, K+pi-pi0 [-0.062, 0.049] GeV; applied in ROOT fit stage.")
  .note(:tag_mBC_signal_region,
        "Tag M_BC signal region 1.859 < M_BC < 1.873 GeV/c^2 for ST yield extraction.")
  .note(:Ks_reconstruction,
        "K_S0 candidate reconstructed from pi+ pi- pair with secondary-vertex-fit and mass window ~M(K_S0). Handled inside DTagTool for K_S0-containing tag modes; signal-side K_S0 handled analogously.")
  .note(:Umiss_signal_region,
        "Signal region -0.05 < U_miss < 0.05 GeV for the electron mode.")
  .note(:signal_yield_fit,
        "Signal yield from unbinned max-likelihood fit to U_miss with MC signal shape convolved with a Gaussian, background from inclusive MC.")
  .note(:E_gamma_max_veto,
        "E_gamma_max cut (energy of most energetic extra photon) applied to suppress residual photon backgrounds; adopted from Ref. [12,13].")
  .note(:M_Kspi_le_cut,
        "M(K_S0 pi- l+) cut applied to suppress D0 -> K_S0 pi- pi+ pi0 peaking background.")
  .apply
alg_e.execute_on([psi3770_data, psi3770_incMC, exMC_e])

######################################################################
### Muon mode : D0 -> K_S0 pi- mu+ nu_mu                            ###
######################################################################
alg_mu = TagAnalysis.new("D0KsPiMuNu")
alg_mu.set_header(["D0KsPiMuNuAlg/D0KsPiMuNu.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_mu)

alg_mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm 1
end

alg_mu.signal_side do |s|
  s.charged(pip: 1, pim: 2, mup: 1)  # K_S0 pi+ pi- + signal pi- + mu+
  s.require_charge 0
  s.missing :nu_mu                   # massless neutrino
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu
  .note(:tag_deltaE_windows,
        "Mode-dependent tag deltaE windows: K+pi- [-0.027, 0.027] GeV, K+pi-pi+pi- [-0.026, 0.024] GeV, K+pi-pi0 [-0.062, 0.049] GeV.")
  .note(:tag_mBC_signal_region,
        "Tag M_BC signal region 1.859 < M_BC < 1.873 GeV/c^2.")
  .note(:Umiss_signal_region,
        "Signal region -0.02 < U_miss < 0.02 GeV for the muon mode (tighter than electron mode).")
  .note(:M_Kspi_mu_cut,
        "M(K_S0 pi- mu+ (pi0)) cut applied to suppress D0 -> K_S0 pi- pi+ pi0 peaking background in the muon channel.")
  .note(:E_gamma_max_veto,
        "E_gamma_max cut applied to suppress extra-photon backgrounds.")
  .apply
alg_mu.execute_on([psi3770_data, psi3770_incMC, exMC_mu])
