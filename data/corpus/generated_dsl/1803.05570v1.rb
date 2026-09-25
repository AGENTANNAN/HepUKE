# frozen_string_literal: true

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at 3.773 GeV (2.93 fb^-1)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC at 3.773 GeV

# ---- Decay cards for the D+ semileptonic signal chains (tag D- decays hadronically) ----
# Chain I: D+ -> eta(->gamma gamma) e+ nu_e
decay_card_eta_2gamma = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 eta e+ nu_e  PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Chain II: D+ -> eta(->pi+ pi- pi0) e+ nu_e
decay_card_eta_3pi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 eta e+ nu_e  PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-  PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Chain III: D+ -> eta'(->pi+ pi- eta, eta->gamma gamma) e+ nu_e
decay_card_etap_2gamma = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 eta' e+ nu_e  PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC samples (200k events each) ----
exMC_eta_2gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToEtaEtaNu_eta2gam"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_eta_2gamma
  config.cross_section   = :default
end

exMC_eta_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToEtaEtaNu_eta3pi"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_eta_3pi
  config.cross_section   = :default
end

exMC_etap_2gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToEtapEtaNu_etap2pipiEta2gam"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_etap_2gamma
  config.cross_section   = :default
end

### ================= Chain I : D+ -> eta(->gamma gamma) e+ nu_e ================= ###
alg_eta_2gamma = TagAnalysis.new("Dtag_SL_Eta2Gam")
alg_eta_2gamma.set_header(["Dtag_SL_Eta2GamAlg/Dtag_SL_Eta2Gam.h"])
               .set_constant({ "ECMS" => [:double, 3.773] })
               .with_decay_card(decay_card_eta_2gamma)

# Tag side: fully reconstructed hadronic D- (single tag). charm -1 pins the D-.
alg_eta_2gamma.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :mBC, min: 1.86, max: 1.88     # explicit tag M_BC window request
end

# Signal side: eta -> gamma gamma, one e+ (charge +1), missing nu_e
alg_eta_2gamma.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: constrain gamma gamma to nominal eta mass + 4-momentum conservation
alg_eta_2gamma.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta_2gamma
  .note(:tag_reconstruction, "D- tag fully reconstructed by DTagAlg from six hadronic modes (K+pi-pi-, K+pi-pi-pi0, K_S0 pi-, K_S0 pi-pi0, K_S0 pi+pi-pi-, K+K-pi-). Internal charged-track cuts (|cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm), pion/kaon likelihood PID (L(pi)>L(K) for pi, L(K)>L(pi) for K), photon cuts (E>25 MeV barrel / >50 MeV endcap, EMC time<700 ns), pi0->gamma gamma mass (0.115,0.150) GeV with 1C mass constraint, K_S0->pi+pi- vertex chi2<100 with mass (0.487,0.511) GeV, and the mode-dependent |deltaE|<3sigma tag window are applied inside the tag reconstruction and are not expressible in the DSL.")
  .note(:signal_positron_pid, "positron must satisfy L(e)>0 and L(e)/(L(e)+L(pi)+L(K))>0.8, with E_EMC/p>0.8 (eta channels); the DSL signal-side ep key uses fixed SimplePIDSvc thresholds, so these tighter positron requirements must be re-applied in the reconstruction.")
  .note(:background_veto, "no extra unused EMC shower above 250 MeV allowed in the event.")
  .note(:efficiency_curve, "signal eta->gamma gamma mass window (0.50,0.58) GeV/c^2 with a 1C fit to the nominal eta mass and chi2<20, selecting the best-chi2 candidate, applied after the nominal kinematic fit (ROOT level).")

alg_eta_2gamma.apply
alg_eta_2gamma.execute_on([data_3773, incMC_3773, exMC_eta_2gamma])

### ================= Chain II : D+ -> eta(->pi+ pi- pi0) e+ nu_e ================= ###
alg_eta_3pi = TagAnalysis.new("Dtag_SL_Eta3Pi")
alg_eta_3pi.set_header(["Dtag_SL_Eta3PiAlg/Dtag_SL_Eta3Pi.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .with_decay_card(decay_card_eta_3pi)

alg_eta_3pi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :mBC, min: 1.86, max: 1.88
end

# Signal side: eta -> pi+ pi- pi0 (pi0 -> gamma gamma), one e+ (charge +1), missing nu_e
alg_eta_3pi.signal_side do |s|
  s.photons 2
  s.charged(ep: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: constrain gamma gamma to nominal pi0 mass + 4-momentum conservation
alg_eta_3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_eta_3pi
  .note(:tag_reconstruction, "same DTagAlg D- tag as Chain I: charged-track and photon quality cuts, pion/kaon likelihood PID, pi0 1C mass constraint, K_S0 vertex fit and mass window, and the mode-dependent |deltaE|<3sigma window are internal to the tag reconstruction.")
  .note(:signal_positron_pid, "positron L(e)>0, L(e)/(L(e)+L(pi)+L(K))>0.8 and E_EMC/p>0.8 are signal-side requirements not expressible by the fixed-threshold ep key.")
  .note(:background_veto, "no extra unused EMC shower above 250 MeV allowed in the event.")
  .note(:efficiency_curve, "signal eta->pi+pi-pi0 mass window (0.52,0.58) GeV/c^2, choosing the candidate closest to the nominal eta mass, applied after the nominal kinematic fit (ROOT level).")

alg_eta_3pi.apply
alg_eta_3pi.execute_on([data_3773, incMC_3773, exMC_eta_3pi])

### ============ Chain III : D+ -> eta'(->pi+ pi- eta, eta->gamma gamma) e+ nu_e ============ ###
alg_etap_2gamma = TagAnalysis.new("Dtag_SL_Etap2PiPiEta2Gam")
alg_etap_2gamma.set_header(["Dtag_SL_Etap2PiPiEta2GamAlg/Dtag_SL_Etap2PiPiEta2Gam.h"])
                .set_constant({ "ECMS" => [:double, 3.773] })
                .with_decay_card(decay_card_etap_2gamma)

alg_etap_2gamma.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :mBC, min: 1.86, max: 1.88
end

# Signal side: eta' -> pi+ pi- eta (eta -> gamma gamma), one e+ (charge +1), missing nu_e
alg_etap_2gamma.signal_side do |s|
  s.photons 2
  s.charged(ep: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: constrain gamma gamma to nominal eta mass + 4-momentum conservation
alg_etap_2gamma.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_etap_2gamma
  .note(:tag_reconstruction, "same DTagAlg D- tag as Chain I: charged-track and photon quality cuts, pion/kaon likelihood PID, pi0 1C mass constraint, K_S0 vertex fit and mass window, and the mode-dependent |deltaE|<3sigma window are internal to the tag reconstruction.")
  .note(:signal_positron_pid, "positron L(e)>0, L(e)/(L(e)+L(pi)+L(K))>0.8, E_EMC/p>0.6 and p>0.2 GeV/c (eta' channel) are signal-side requirements not expressible by the fixed-threshold ep key.")
  .note(:background_veto, "no extra unused EMC shower above 250 MeV allowed in the event.")
  .note(:signal_mass_window, "eta'(->pi+pi-eta_2gamma) mass window (0.935,0.980) GeV/c^2 and eta'(->pi+pi-eta_3pi) mass window (0.930,0.980) GeV/c^2; the additional eta'->gamma rho0 (rho0->pi+pi-) sub-mode with window (0.55,0.90) GeV/c^2 requires radiative photon E>0.1 GeV, |cos(theta_{pi,rho})|<0.85, angle to e+ >0.20 rad and angle to tag tracks >0.52 rad. These sub-mode selections are applied after the nominal kinematic fit (ROOT level).")

alg_etap_2gamma.apply
alg_etap_2gamma.execute_on([data_3773, incMC_3773, exMC_etap_2gamma])