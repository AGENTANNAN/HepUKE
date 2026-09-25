# Paper: 2109.12812v1
# Title: J/psi -> gamma eta and absolute eta decay branching fractions
# Dataset: J/psi (708_3097), 1.0087e10 events
#
# Analysis structure:
#   1. Inclusive channel: J/psi -> gamma eta, gamma -> e+e- (PCF) — mostly inexpressible
#   2. Exclusive channel A: J/psi -> gamma eta, eta -> gamma gamma (4C fit, chi2<80)
#   3. Exclusive channel B: J/psi -> gamma eta, eta -> pi0 pi0 pi0 (7C fit, chi2<100)
#   4. Exclusive channel C: J/psi -> gamma eta, eta -> pi+ pi- pi0 (5C fit, chi2<100)
#   5. Exclusive channel D: J/psi -> gamma eta, eta -> pi+ pi- gamma (4C fit, chi2<60)
#      + competing hypothesis: 5C fit under pi+ pi- pi0 hypothesis

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi  = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Algorithm 0: Inclusive channel — J/psi -> gamma eta (PCF)
# The Photon Conversion Finder (PCF) package is specialized BOSS
# tooling for reconstructing converted photons (gamma -> e+e-).
# This channel is largely inexpressible in the current DSL.
# ============================================================

alg_incl = Algorithm.new("JPsiGammaEtaIncl")
alg_incl.set_header(["JPsiGammaEtaInclAlg/JPsiGammaEtaIncl.h"])
        .set_constant({ "ECMS" => [:double, 3.097] })
alg_incl.note(:pcf, "Photon Conversion Finder (PCF) used to reconstruct radiative photon " \
                      "from e+e- conversion pair. PCF-specific cuts: R_xy>2cm, |Delta_xy|<0.2cm, " \
                      "Delta_z<1.5cm, Psi_pair in [-0.5,0.5], cos_theta_eg>0.8. " \
                      "At least one extra EMC photon required. Photon energy < 1.4 GeV. " \
                      "cos_theta_gammagamma in [-0.998,0] for N_gamma<5. " \
                      "|cos_theta_miss|<0.98 for 2-track <4-photon events.")
alg_incl.note(:inclusive_method, "Inclusive eta tagging via recoil mass of e+e- pair. " \
                                  "Signal yield from unbinned ML fit to M_recoil(e+e-). " \
                                  "Backgrounds: e+e-->gammagamma, e+e-->e+e-, J/psi->e+e-eta, " \
                                  "J/psi->pi+pi-pi0, J/psi->omega eta, J/psi->omega pi0, J/psi->gamma eta'.")
alg_incl.apply(nil)  # No Selection — PCF is inexpressible
alg_incl.execute_on([data_jpsi, incMC_jpsi])

# ============================================================
# Decay cards for exclusive channels
# ============================================================

# Channel A: eta -> gamma gamma
decay_eta_gg = <<~DECAYCARD
  Decay J/psi
  1 gamma eta HELAMP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Channel B: eta -> pi0 pi0 pi0
decay_eta_3pi0 = <<~DECAYCARD
  Decay J/psi
  1 gamma eta HELAMP;
  Enddecay
  Decay eta
  1 pi0 pi0 pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Channel C: eta -> pi+ pi- pi0
decay_eta_pipipim_pi0 = <<~DECAYCARD
  Decay J/psi
  1 gamma eta HELAMP;
  Enddecay
  Decay eta
  1 pi+ pi- pi0 D_DALITZ;
  Enddecay
  End
DECAYCARD

# Channel D: eta -> pi+ pi- gamma
decay_eta_pipi_gamma = <<~DECAYCARD
  Decay J/psi
  1 gamma eta HELAMP;
  Enddecay
  Decay eta
  1 pi+ pi- gamma BOX_ANOMALY;
  Enddecay
  End
DECAYCARD

# Signal MC for all exclusive channels (same J/psi -> gamma eta topology)
sig_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_gg"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_eta_gg
  config.cross_section   = :default
end

sig_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_3pi0"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_eta_3pi0
  config.cross_section   = :default
end

sig_eta_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_pipipi0"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_eta_pipipim_pi0
  config.cross_section   = :default
end

sig_eta_pipi_g = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_pipig"
  config.related_dataset = data_jpsi
  config.events          = 1_000_000
  config.decay_card      = decay_eta_pipi_gamma
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: Exclusive channel A — J/psi -> gamma eta, eta -> gamma gamma
# 4C kinematic fit, chi2_4C < 80
# Photon energy > 0.07 GeV to suppress e+e- -> gamma gamma(gamma)
# Barrel EMC only: |cos_theta| < 0.80
# All gamma-gamma combinations kept (radiative photon not separable)
# ============================================================

alg_eta_gg = Algorithm.new("JPsiGammaEta_GG")
alg_eta_gg.set_header(["JPsiGammaEtaGGAlg/JPsiGammaEtaGG.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })

sel_eta_gg = Selection.new
  .select_track { nChrp "==0"; nChrn "==0" }
  .select_photon { nGam ">=3"; energyThreshold_b 0.07; energyThreshold_e 0.07 }
  .kinematic_fit([:gamma, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 80
    nominal
  end

alg_eta_gg.note(:photon_selection, "Only barrel EMC photons (|cos_theta|<0.80). " \
                                    "EMC time within [-500,500] ns of most energetic photon. " \
                                    "Radiative photon is most energetic photon (E_gamma >> E_eta-decay). " \
                                    "All gamma-gamma combinations kept; wrong combinations produce flat background.")
alg_eta_gg.with_decay_card(decay_eta_gg).apply(sel_eta_gg)
alg_eta_gg.execute_on([data_jpsi, incMC_jpsi, sig_eta_gg])

# ============================================================
# Algorithm 2: Exclusive channel B — J/psi -> gamma eta, eta -> pi0 pi0 pi0
# 7C kinematic fit (4C + 3 pi0 mass constraints), chi2_7C < 100
# 6 photons minimum for 3pi0
# ============================================================

alg_eta_3pi0 = Algorithm.new("JPsiGammaEta_3Pi0")
alg_eta_3pi0.set_header(["JPsiGammaEta3Pi0Alg/JPsiGammaEta3Pi0.h"])
            .set_constant({ "ECMS" => [:double, 3.097] })

sel_eta_3pi0 = Selection.new
  .select_track { nChrp "==0"; nChrn "==0" }
  .select_photon { nGam ">=6"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 100
    nominal
  end

alg_eta_3pi0.note(:photon_selection, "Barrel EMC only. EMC time within [-500,500] ns of most energetic photon. " \
                                     "Radiative photon is most energetic. Best 3pi0 combination chosen by min chi2_7C.")
alg_eta_3pi0.note(:peaking_bg, "Peaking backgrounds: J/psi -> omega eta (omega->gamma pi0), " \
                               "J/psi -> gamma f0(2100) (f0->eta eta).")
alg_eta_3pi0.with_decay_card(decay_eta_3pi0).apply(sel_eta_3pi0)
alg_eta_3pi0.execute_on([data_jpsi, incMC_jpsi, sig_eta_3pi0])

# ============================================================
# Algorithm 3: Exclusive channel C — J/psi -> gamma eta, eta -> pi+ pi- pi0
# 5C kinematic fit (4C + pi0 mass constraint), chi2_5C < 100
# 2 charged tracks (pi+ pi-) + 2 photons (pi0) + radiative photon
# ============================================================

alg_eta_pipipi0 = Algorithm.new("JPsiGammaEta_PiPiPi0")
alg_eta_pipipi0.set_header(["JPsiGammaEtaPiPiPi0Alg/JPsiGammaEtaPiPiPi0.h"])
               .set_constant({ "ECMS" => [:double, 3.097] })

sel_eta_pipipi0 = Selection.new
  .select_track { nChrp "==1"; nChrn "==1"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=3"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 100
    nominal
  end

alg_eta_pipipi0.note(:photon_selection, "Barrel EMC only. Radiative photon is most energetic photon. " \
                                        "Best combination chosen by min chi2_5C.")
alg_eta_pipipi0.with_decay_card(decay_eta_pipipim_pi0).apply(sel_eta_pipipi0)
alg_eta_pipipi0.execute_on([data_jpsi, incMC_jpsi, sig_eta_pipipi0])

# ============================================================
# Algorithm 4: Exclusive channel D — J/psi -> gamma eta, eta -> pi+ pi- gamma
# 4C kinematic fit, chi2_4C < 60
# Competing hypothesis: if >2 good photons, 5C fit under pi+pi-pi0 hypothesis
# Require P(gamma pi+ pi- gamma) > P(gamma pi+ pi- pi0)
# ============================================================

alg_eta_pipig = Algorithm.new("JPsiGammaEta_PiPiGamma")
alg_eta_pipig.set_header(["JPsiGammaEtaPiPiGammaAlg/JPsiGammaEtaPiPiGamma.h"])
             .set_constant({ "ECMS" => [:double, 3.097] })

# Decay card: same J/psi -> gamma eta, eta -> pi+ pi- gamma
# But with 5C competing hypothesis for pi+ pi- pi0 background veto

sel_eta_pipig = Selection.new
  .select_track { nChrp "==1"; nChrn "==1"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=2"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  # Nominal 4C fit under pi+ pi- gamma hypothesis
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 60
    nominal
  end
  # Competing hypothesis: 5C fit under pi+ pi- pi0 hypothesis (Rule T2)
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) do
    use_track_index_from_nominal_kmfit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    # No nominal, no chi2_cut — store chi2 for ROOT-level comparison
  end

alg_eta_pipig.note(:photon_selection, "Barrel EMC only. Radiative photon = most energetic photon. " \
                                      "If >2 good photons: perform 5C fit under pi+pi-pi0 hypothesis, " \
                                      "require P(gamma pi+ pi- gamma) > P(gamma pi+ pi- pi0).")
alg_eta_pipig.note(:bg, "Background from eta -> pi+ pi- pi0 forms broad bump on left side of eta peak.")
alg_eta_pipig.with_decay_card(decay_eta_pipi_gamma).apply(sel_eta_pipig)
alg_eta_pipig.execute_on([data_jpsi, incMC_jpsi, sig_eta_pipi_g])