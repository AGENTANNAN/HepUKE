# BESIII: Study of D+ -> eta(') e+ nu_e
# Data: 2.93 fb^-1 at sqrt(s)=3.773 GeV (psi(3770))
# Double-tag technique: ST D- hadronic tag + SL D+ signal

### Dataset description ###
psi3770_data   = DatasetManager.real_data.find("712_3773")
psi3770_incMC  = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for D+ -> eta e+ nu_e (eta -> gamma gamma)
decay_card_eta_gg = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 eta e+ nu_e                        PHOTOS ISGW2;
  Enddecay

  Decay eta
  1.0000 gamma gamma                        PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for D+ -> eta e+ nu_e (eta -> pi+ pi- pi0)
decay_card_eta_3pi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 eta e+ nu_e                        PHOTOS ISGW2;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0                        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for D+ -> eta' e+ nu_e (eta' -> pi+ pi- eta, eta -> gamma gamma)
decay_card_etap_pipieta_gg = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 eta' e+ nu_e                       PHOTOS ISGW2;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta                        PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                        PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC samples
exMC_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_eta_gg_enu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_eta_gg
  config.cross_section = :default
end

exMC_eta_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_eta_3pi_enu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_eta_3pi
  config.cross_section = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_etap_enu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_etap_pipieta_gg
  config.cross_section = :default
end

# =============================================================================
# TagAnalysis: D+ -> eta e+ nu_e  (eta -> gamma gamma)
# ST D- hadronic tag + SL signal: eta(->gamma gamma) + e+ + missing nu_e
# =============================================================================
tag_eta_gg = TagAnalysis.new("DpEtaGGEnu", '00-00-01')
tag_eta_gg.set_header(["DpEtaGGEnuAlg/DpEtaGGEnu.h"])
            .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: reconstruct D- via 6 hadronic modes (paper lists 6 tag modes)
tag_eta_gg.tag_side(:Dplus) do |t|
  t.modes(
    :DptoKPiPi,        # K+ pi- pi-
    :DptoKPiPiPi0,     # K+ pi- pi- pi0
    :DptoKsPi,          # K_S0 pi-
    :DptoKsPiPi0,       # K_S0 pi- pi0
    :DptoKsPiPiPi,     # K_S0 pi+ pi- pi-
    :DptoKKPi           # K+ K- pi-
  )
  t.charm(-1)   # tag D-
end

# Signal side: eta -> gamma gamma + e+ + missing nu_e
tag_eta_gg.signal_side do |s|
  s.photons 2                     # two photons for eta -> gamma gamma
  s.charged(ep: 1)                # one positron
  s.missing :nu_e                 # missing neutrino (massless)
  s.require_charge 1              # signal side charge = +1 (e+)
  s.min_photon_energy 0.025       # photon energy > 25 MeV barrel
end

# Kinematic fit: constrain eta mass and four-momentum
tag_eta_gg.fit do |f|
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_eta_gg.with_decay_card(decay_card_eta_gg).apply
tag_eta_gg
  .note(:tag_mode_unavailable, "Paper also uses D+->K+pi-pi-pi0 tag mode which includes pi0 reconstruction. The DSL mode :DptoKPiPiPi0 corresponds to this mode; if the mode is unavailable in the current BOSS version, tag mode count reduces accordingly.")
  .note(:st_selection, "Charged tracks: |cos(theta)|<0.93, Vz<10cm, Vr<1cm. PID: L(pi)>L(K) for pions, L(K)>L(pi) for kaons. Photons: E>25MeV barrel, E>50MeV endcap, EMC time<700ns. pi0->gamma gamma: mass in (0.115,0.150) GeV/c^2 with 1C kinematic mass constraint. K_S0->pi+pi-: vertex fit chi2<100, mass (0.487,0.511) GeV/c^2. DeltaE cuts: mode-dependent at +/-3sigma. M_BC signal region [1.86,1.88] GeV/c^2.")
  .note(:d_eta_gg_signal, "eta->gamma gamma: mass (0.50,0.58) GeV/c^2, 1C fit constraining to nominal eta mass, chi2<20. Multiple candidates: keep best chi2.")
  .note(:d_eta_3pi_signal, "eta->pi+pi-pi0: mass (0.52,0.58) GeV/c^2. Multiple candidates: keep closest to nominal eta mass. (Uses separate TagAnalysis or combined at ROOT level due to different track/photon requirements)")
  .note(:d_etap_signal, "eta' -> pi+pi-eta_2gamma: mass (0.935,0.980) GeV/c^2; eta' -> pi+pi-eta_3pi: mass (0.930,0.980) GeV/c^2; eta' -> gamma rho0, rho0->pi+pi-: mass (0.55,0.90) GeV/c^2, radiative photon E>0.1GeV, |cos(theta_pi_rho)|<0.85, angle to positron>0.20 rad, angle to tag tracks>0.52 rad.")
  .note(:positron_selection, "Positron: L(e)>0 and L(e)/(L(e)+L(pi)+L(K))>0.8. E_EMC/p > 0.8 for eta channel, > 0.6 for eta' channel. p>0.2 GeV/c for eta' channel. Extra shower veto: no unused EMC shower with E>250 MeV.")
  .note(:U_miss, "U_miss = E_miss - c|p_miss| where E_miss = E_beam - E_eta(') - E_e+ and p_miss = p_D+ - p_eta(') - p_e+. U_miss peaks at zero for signal. DT yield from simultaneous unbinned fit over decay modes.")
  .note(:branching_fraction, "B(D+->eta e+ nu_e) = (10.74 +/- 0.81 +/- 0.51)*10^-4. B(D+->eta' e+ nu_e) = (1.91 +/- 0.51 +/- 0.13)*10^-4. Total systematic: 4.7% (eta), 6.9% (eta').")
  .note(:form_factor, "Form factor f_+(0)|V_cd| measured in 3 q^2 bins for D+->eta e+ nu_e. Three parameterizations: simple pole, modified pole, series expansion. Results in paper Table IV.")
tag_eta_gg.execute_on([psi3770_data, psi3770_incMC, exMC_eta_gg])

# =============================================================================
# TagAnalysis: D+ -> eta e+ nu_e  (eta -> pi+ pi- pi0)
# Same tag side, signal side with 4 charged + 2 photons
# =============================================================================
tag_eta_3pi = TagAnalysis.new("DpEta3PiEnu", '00-00-01')
tag_eta_3pi.set_header(["DpEta3PiEnuAlg/DpEta3PiEnu.h"])
             .set_constant({"ECMS" => [:double, 3.773]})

tag_eta_3pi.tag_side(:Dplus) do |t|
  t.modes(:DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi)
  t.charm(-1)
end

tag_eta_3pi.signal_side do |s|
  s.photons 2                     # two photons for pi0 -> gamma gamma
  s.charged(pip: 1, pim: 1, ep: 1)  # pi+ pi- for eta + e+
  s.missing :nu_e
  s.require_charge 1
  s.min_photon_energy 0.025
end

# Fit: pi0 mass from gammas, eta mass from pi+pi-pi0, 4-momentum
tag_eta_3pi.fit do |f|
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_eta_3pi.with_decay_card(decay_card_eta_3pi).apply
tag_eta_3pi
  .note(:signal_selection, "Same ST selection and systematics as DpEtaGGEnu. See DpEtaGGEnu for shared details.")
  .note(:eta_3pi_reco_note, "eta -> pi+ pi- pi0 is reconstructed by applying a mass window (0.52, 0.58) GeV/c^2 after the pi0 constraint. The invariant_mass_of(:pip, :pim, :pi0) -> eta constraint is applied at ROOT level if not expressible in the TagAnalysis fit block.")
tag_eta_3pi.execute_on([psi3770_data, psi3770_incMC, exMC_eta_3pi])

# =============================================================================
# TagAnalysis: D+ -> eta' e+ nu_e  (eta' -> pi+ pi- eta, eta -> gamma gamma)
# =============================================================================
tag_etap = TagAnalysis.new("DpEtapEnu", '00-00-01')
tag_etap.set_header(["DpEtapEnuAlg/DpEtapEnu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})

tag_etap.tag_side(:Dplus) do |t|
  t.modes(:DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi)
  t.charm(-1)
end

tag_etap.signal_side do |s|
  s.photons 2                      # two photons for eta -> gamma gamma
  s.charged(pip: 2, pim: 1, ep: 1)  # pi+ pi+ pi- for eta' recoil + e+
  s.missing :nu_e
  s.require_charge 1
  s.min_photon_energy 0.025
end

tag_etap.fit do |f|
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_etap.with_decay_card(decay_card_etap_pipieta_gg).apply
tag_etap
  .note(:etap_decay_modes, "Three eta' decay modes studied: (1) eta'->pi+pi-eta_2gamma (mass 0.935-0.980 GeV/c^2), (2) eta'->pi+pi-eta_3pi (mass 0.930-0.980 GeV/c^2), (3) eta'->gamma rho0, rho0->pi+pi- (mass 0.55-0.90 GeV/c^2, radiative gamma E>0.1GeV). Only mode (1) is expressed in the DSL; modes (2) and (3) require separate TagAnalysis instances due to different particle requirements or custom reconstruction.")
  .note(:signal_selection, "Same ST selection and systematics as DpEtaGGEnu. See DpEtaGGEnu notes for shared details.")
  .note(:positron_selection, "For eta' channel: E_EMC/p > 0.6, p>0.2 GeV/c to reduce mis-PID. Angles between radiative photon and all charged tracks in tag final state > 0.52 rad.")
  .note(:branching_fraction_etap, "B(D+->eta' e+ nu_e) = (1.91 +/- 0.51 +/- 0.13)*10^-4.")
tag_etap.execute_on([psi3770_data, psi3770_incMC, exMC_etap])