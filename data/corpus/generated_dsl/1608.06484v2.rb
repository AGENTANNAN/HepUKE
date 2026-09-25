# =============================================================================
# Dataset preparation — sqrt(s) = 4.009 GeV (psi(4040)), L = 482 pb^-1
# =============================================================================
# BESIII sample-name convention: [BOSS version]_[CMS energy in MeV] -> 703_4009
data_4009  = DatasetManager.real_data.find("703_4009")        # 482 pb^-1 real data at 4.009 GeV
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")     # corresponding inclusive MC

# ---------------------------------------------------------------------------
# Signal decay cards (EvtGen syntax), one per signal mode
# ---------------------------------------------------------------------------
# Mode I: Ds+ -> eta e+ nu_e ,  eta -> gamma gamma
decay_card_eta_ev = <<~DECAYCARD
  Decay psi(4040)
  1.0000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.0000 eta e+ nu_e PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: Ds+ -> eta' e+ nu_e , eta' -> eta pi+ pi- , eta -> gamma gamma
decay_card_etap_etapipi_ev = <<~DECAYCARD
  Decay psi(4040)
  1.0000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.0000 eta' e+ nu_e PHSP;
  Enddecay

  Decay eta'
  1.0000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: Ds+ -> eta' e+ nu_e , eta' -> gamma rho0 , rho0 -> pi+ pi-
decay_card_etap_gammarho_ev = <<~DECAYCARD
  Decay psi(4040)
  1.0000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.0000 eta' e+ nu_e PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma rho0 PHSP;
  Enddecay

  Decay rho0
  1.0000 pi+ pi- VSS;
  Enddecay

  End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples: 500k events per signal mode, matched to the 4.009 data
# ---------------------------------------------------------------------------
exMC_eta_ev = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ds_semilep_eta_ev"
  config.related_dataset = data_4009
  config.events          = 500_000
  config.decay_card      = decay_card_eta_ev
  config.cross_section   = :default
end

exMC_etap_etapipi_ev = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ds_semilep_etap_etapipi_ev"
  config.related_dataset = data_4009
  config.events          = 500_000
  config.decay_card      = decay_card_etap_etapipi_ev
  config.cross_section   = :default
end

exMC_etap_gammarho_ev = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ds_semilep_etap_gammarho_ev"
  config.related_dataset = data_4009
  config.events          = 500_000
  config.decay_card      = decay_card_etap_gammarho_ev
  config.cross_section   = :default
end

# =============================================================================
# Event selection (BOSS) — double tag: Ds- tag (hadronic) + Ds+ semileptonic
# =============================================================================
# ---------------------------------------------------------------------------
# Mode I : Ds+ -> eta (-> gamma gamma) e+ nu_e
# ---------------------------------------------------------------------------
alg_eta_ev = TagAnalysis.new("DsTagEtaEv")
alg_eta_ev.set_header(["DsTagEtaEvAlg/DsTagEtaEv.h"])
          .set_constant({"ECMS" => [:double, 4.009]})
          .note(:pid_correction_method,
                "Tag-side reconstruction is delegated to DTagAlg (EvtRecDTag): charged tracks |cos(theta)|<0.93, |Vz|<10 cm, |Vxy|<1 cm, pion/kaon separation by likelihood (pion if CL_pi>0 and CL_pi>CL_K; kaon if CL_K>0 and CL_K>CL_pi); photons with TDC in [0,700] ns, >10 deg from any charged track and E>25 MeV (barrel, |cos(theta)|<0.80) or E>50 MeV (endcap, 0.86<|cos(theta)|<0.92). Tag-side resonance windows: 0.115<M(gamma gamma)<0.150 (pi0), 0.510<M(gamma gamma)<0.570 (eta), 1.005<M(K+ K-)<1.040 (phi), 0.570<M(pi0 pi-)<0.970 (rho-), 0.943<M(eta pi+ pi-)<0.973 (eta'->eta pi pi), 0.932<M(gamma rho0)<0.980 (eta'->gamma rho0), K_S0 -> pi+ pi- with |Vz|<20 cm, 0.487<M(pi+ pi-)<0.511 and a positive decay length; pi0 and eta mass-constrained (1C).")
          .note(:efficiency_curve,
                "Tag candidate selected as the one with minimum |deltaE|, |deltaE| required within 3 sigma of the fitted peak; mBC signal region (-4 sigma, +5 sigma) about the fitted peak. Tag mBC/deltaE are stored and windowed in the ROOT analysis.")
          .note(:background_veto,
                "Extra photon with energy > 300 MeV vetoed; no additional charged tracks allowed. Final signal window on U_miss = (-0.10, +0.12) GeV applied in ROOT.")
          .with_decay_card(decay_card_eta_ev)

# Tag side: Ds- fully reconstructed in hadronic modes (single tag)
alg_eta_ev.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,              # Ds- -> K+ K- pi-
          :DstoPhiRho,            # Ds- -> phi rho-
          :DstoKsKPiPi,           # Ds- -> K_S0 K+ pi- pi-
          :DstoKsKPimPip,         # Ds- -> K_S0 K- pi+ pi-
          :DstoKsK,               # Ds- -> K_S0 K-
          :DstoPiPiPi,            # Ds- -> pi+ pi- pi-
          :DstoEtaPi,             # Ds- -> eta pi-
          :DstoEtaPrimePi,        # Ds- -> eta' pi-, eta' -> eta pi pi
          :DstoEtaPrimePiGamRho,  # Ds- -> eta' pi-, eta' -> gamma rho0
          :DstoEtaRho             # Ds- -> eta rho-
  t.charm(-1)                     # tag the Ds- side
end

# Signal side: everything the tag did not use -> eta -> gamma gamma, e+, nu_e
alg_eta_ev.signal_side do |s|
  s.photons 2                 # two photons from eta -> gamma gamma
  s.min_photon_angle 10.0     # photon > 10 deg from any charged track
  s.min_photon_energy 0.025   # EMC shower energy floor
  s.charged(ep: 1)            # the signal electron (opposite charge to the tag)
  s.require_charge 1          # net signal-side charge +1
  s.missing :nu_e             # missing neutrino (massless)
end

# Kinematic fit: 4-momentum conservation + eta -> gamma gamma mass constraint
alg_eta_ev.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta_ev.apply
alg_eta_ev.execute_on([data_4009, incMC_4009, exMC_eta_ev])

# ---------------------------------------------------------------------------
# Mode II : Ds+ -> eta' (-> eta pi+ pi-, eta -> gamma gamma) e+ nu_e
# ---------------------------------------------------------------------------
alg_etap_etapipi_ev = TagAnalysis.new("DsTagEtaPrimeEtaPiPiEv")
alg_etap_etapipi_ev.set_header(["DsTagEtaPrimeEtaPiPiEvAlg/DsTagEtaPrimeEtaPiPiEv.h"])
                   .set_constant({"ECMS" => [:double, 4.009]})
                   .note(:pid_correction_method,
                         "Tag-side reconstruction is delegated to DTagAlg (EvtRecDTag): charged tracks |cos(theta)|<0.93, |Vz|<10 cm, |Vxy|<1 cm, pion/kaon separation by likelihood (pion if CL_pi>0 and CL_pi>CL_K; kaon if CL_K>0 and CL_K>CL_pi); photons with TDC in [0,700] ns, >10 deg from any charged track and E>25 MeV (barrel) or E>50 MeV (endcap). Tag-side resonance windows: 0.115<M(gamma gamma)<0.150 (pi0), 0.510<M(gamma gamma)<0.570 (eta), 1.005<M(K+ K-)<1.040 (phi), 0.570<M(pi0 pi-)<0.970 (rho-), 0.943<M(eta pi+ pi-)<0.973 (eta'->eta pi pi), 0.932<M(gamma rho0)<0.980 (eta'->gamma rho0), K_S0 -> pi+ pi- with |Vz|<20 cm, 0.487<M(pi+ pi-)<0.511 and a positive decay length; pi0 and eta mass-constrained (1C).")
                   .note(:efficiency_curve,
                         "Tag candidate selected as the one with minimum |deltaE|, |deltaE| required within 3 sigma of the fitted peak; mBC signal region (-4 sigma, +5 sigma) about the fitted peak. Tag mBC/deltaE are stored and windowed in the ROOT analysis.")
                   .note(:background_veto,
                         "Extra photon with energy > 300 MeV vetoed; no additional charged tracks and no extra pi0 allowed. Final signal window on U_miss = (-0.10, +0.12) GeV applied in ROOT.")
                   .with_decay_card(decay_card_etap_etapipi_ev)

alg_etap_etapipi_ev.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,              # Ds- -> K+ K- pi-
          :DstoPhiRho,            # Ds- -> phi rho-
          :DstoKsKPiPi,           # Ds- -> K_S0 K+ pi- pi-
          :DstoKsKPimPip,         # Ds- -> K_S0 K- pi+ pi-
          :DstoKsK,               # Ds- -> K_S0 K-
          :DstoPiPiPi,            # Ds- -> pi+ pi- pi-
          :DstoEtaPi,             # Ds- -> eta pi-
          :DstoEtaPrimePi,        # Ds- -> eta' pi-, eta' -> eta pi pi
          :DstoEtaPrimePiGamRho,  # Ds- -> eta' pi-, eta' -> gamma rho0
          :DstoEtaRho             # Ds- -> eta rho-
  t.charm(-1)                     # tag the Ds- side
end

# Signal side: eta' -> eta pi+ pi- (eta -> gamma gamma) + e+ + nu_e
alg_etap_etapipi_ev.signal_side do |s|
  s.photons 2                 # two photons from eta -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(ep: 1, pip: 1, pim: 1)  # e+ from the semileptonic decay, pi+ pi- from eta'
  s.require_charge 1                # net signal-side charge +1
  s.missing :nu_e                   # missing neutrino (massless)
end

alg_etap_etapipi_ev.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_etap_etapipi_ev.apply
alg_etap_etapipi_ev.execute_on([data_4009, incMC_4009, exMC_etap_etapipi_ev])

# ---------------------------------------------------------------------------
# Mode III : Ds+ -> eta' (-> gamma rho0, rho0 -> pi+ pi-) e+ nu_e
# ---------------------------------------------------------------------------
alg_etap_gammarho_ev = TagAnalysis.new("DsTagEtaPrimeGammaRhoEv")
alg_etap_gammarho_ev.set_header(["DsTagEtaPrimeGammaRhoEvAlg/DsTagEtaPrimeGammaRhoEv.h"])
                    .set_constant({"ECMS" => [:double, 4.009]})
                    .note(:pid_correction_method,
                          "Tag-side reconstruction is delegated to DTagAlg (EvtRecDTag): charged tracks |cos(theta)|<0.93, |Vz|<10 cm, |Vxy|<1 cm, pion/kaon separation by likelihood (pion if CL_pi>0 and CL_pi>CL_K; kaon if CL_K>0 and CL_K>CL_pi); photons with TDC in [0,700] ns, >10 deg from any charged track and E>25 MeV (barrel) or E>50 MeV (endcap). Tag-side resonance windows: 0.115<M(gamma gamma)<0.150 (pi0), 0.510<M(gamma gamma)<0.570 (eta), 1.005<M(K+ K-)<1.040 (phi), 0.570<M(pi0 pi-)<0.970 (rho-), 0.943<M(eta pi+ pi-)<0.973 (eta'->eta pi pi), 0.932<M(gamma rho0)<0.980 (eta'->gamma rho0), K_S0 -> pi+ pi- with |Vz|<20 cm, 0.487<M(pi+ pi-)<0.511 and a positive decay length; pi0 and eta mass-constrained (1C).")
                    .note(:efficiency_curve,
                          "Tag candidate selected as the one with minimum |deltaE|, |deltaE| required within 3 sigma of the fitted peak; mBC signal region (-4 sigma, +5 sigma) about the fitted peak. Tag mBC/deltaE are stored and windowed in the ROOT analysis.")
                    .note(:background_veto,
                          "Extra photon with energy > 300 MeV vetoed; no additional charged tracks and no extra pi0 allowed. Final signal window on U_miss = (-0.08, +0.10) GeV applied in ROOT.")
                    .with_decay_card(decay_card_etap_gammarho_ev)

alg_etap_gammarho_ev.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,              # Ds- -> K+ K- pi-
          :DstoPhiRho,            # Ds- -> phi rho-
          :DstoKsKPiPi,           # Ds- -> K_S0 K+ pi- pi-
          :DstoKsKPimPip,         # Ds- -> K_S0 K- pi+ pi-
          :DstoKsK,               # Ds- -> K_S0 K-
          :DstoPiPiPi,            # Ds- -> pi+ pi- pi-
          :DstoEtaPi,             # Ds- -> eta pi-
          :DstoEtaPrimePi,        # Ds- -> eta' pi-, eta' -> eta pi pi
          :DstoEtaPrimePiGamRho,  # Ds- -> eta' pi-, eta' -> gamma rho0
          :DstoEtaRho             # Ds- -> eta rho-
  t.charm(-1)                     # tag the Ds- side
end

# Signal side: eta' -> gamma rho0 (rho0 -> pi+ pi-) + e+ + nu_e
alg_etap_gammarho_ev.signal_side do |s|
  s.photons 1                       # single photon from eta' -> gamma rho0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(ep: 1, pip: 1, pim: 1)  # e+ from the semileptonic decay, pi+ pi- from rho0
  s.require_charge 1                # net signal-side charge +1
  s.missing :nu_e                   # missing neutrino (massless)
end

# Kinematic fit: 4-momentum conservation only (no eta -> gamma gamma here)
alg_etap_gammarho_ev.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap_gammarho_ev.apply
alg_etap_gammarho_ev.execute_on([data_4009, incMC_4009, exMC_etap_gammarho_ev])