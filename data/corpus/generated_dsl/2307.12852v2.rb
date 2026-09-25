# =============================================================================
#  Ds+ -> eta(') mu+ nu_mu   |   tag-based (Ds- tag) missing-mass analysis
#  Eight energy points spanning 4.128 - 4.226 GeV
# =============================================================================

### Dataset description ###
# BOSS 705: 4130 / 4160 MeV ; BOSS 703: 4180 / 4190 / 4200 / 4210 / 4220 / 4230 MeV
scan_points = %w[705_4130 705_4160 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230]

data_points  = scan_points.map { |name| DatasetManager.real_data.find(name) }
incMC_points = scan_points.map { |name| DatasetManager.inclusive_mc.find(name) }

### Decay cards (EvtGen syntax, exclusive signal processes) ###
# Ds+ -> eta mu+ nu_mu, ISGW2 form factors; eta -> gamma gamma treated as phase space
decay_card_eta = <<~DECAYCARD
  Decay D_s+
  1.0000 eta mu+ nu_mu ISGW2;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Ds+ -> eta' mu+ nu_mu, ISGW2 form factors; eta' -> gamma pi+ pi- treated as phase space
decay_card_etap = <<~DECAYCARD
  Decay D_s+
  1.0000 eta' mu+ nu_mu ISGW2;
  Enddecay

  Decay eta'
  1.0000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC: 200k events per mode at every energy point ###
exMC_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_eta_munu"      # suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

exMC_etap = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_etap_munu"     # suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_etap
  config.cross_section = :default
end

### =========================================================================
### Event selection - mode I: Ds+ -> eta(-> gamma gamma) mu+ nu_mu
### =========================================================================
alg_name_eta = "DsToEtaMuNu"
alg_eta = TagAnalysis.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 4.178]})   # nominal energy of the 8-point scan (4.128-4.226 GeV)
       .note(:soft_transition_photon_selection,
             "the soft transition photon from the Ds* is chosen by minimum |DeltaE| inside the " \
             "tagging (DTagAlg reconstruction detail, not tunable in the DSL)")
       .with_decay_card(decay_card_eta)

# --- Ds- tag side: 14 hadronic tag modes (full DTagAlg channel names) ---
alg_eta.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,              # K+K-pi-
          :DstoKKPiPi0,           # K+K-pi-pi0
          :DstoKsK,               # Ks K-
          :DstoKsKPi0,            # Ks K- pi0
          :DstoKsKsPi,            # Ks Ks pi-
          :DstoKsKPiPi,           # Ks K- pi+ pi-
          :DstoPiPiPi,            # pi+ pi- pi-
          :DstoEtaGGPi,           # eta(gamma gamma) pi-
          :DstoEtaPiPiPi0Pi,      # eta(pi+ pi- pi0) pi-
          :DstoEtaPrimeGGPiPiPi,  # eta'(gamma gamma) pi pi pi
          :DstoEtaPrimeGamPiPiPi, # eta'(gamma pi+ pi-) pi
          :DstoEtaGGRho,          # eta(gamma gamma) rho-
          :DstoEtaPiPiPi0Rho      # eta(pi+ pi- pi0) rho-
  t.charm -1   # tagged side pinned to Ds- (tag mBC / DeltaE stored, not cut here)
end

# --- Signal side: everything the tag did not use ---
alg_eta.signal_side do |s|
  s.photons 2                        # eta -> gamma gamma
  s.charged(mup: 1, at_least: true)  # at least one mu+ (extra-track veto is a ROOT-level cut)
  s.missing :nu_mu                   # missing neutrino
  s.min_photon_energy 0.025          # photon energy above 25 MeV
end

# --- 3C kinematic fit: four-momentum conservation + M(gamma gamma) = m(eta) ---
alg_eta.fit do |f|
  f.constrain_four_momentum                                                 # 4C
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)    # 1C -> 3C
  f.chi2_cut 200
end

alg_eta.apply
alg_eta.execute_on(data_points + incMC_points + exMC_eta)

### =========================================================================
### Event selection - mode II: Ds+ -> eta'(-> gamma pi+ pi-) mu+ nu_mu
### =========================================================================
alg_name_etap = "DsToEtaPrimeMuNu"
alg_etap = TagAnalysis.new(alg_name_etap)
alg_etap.set_header(["#{alg_name_etap}Alg/#{alg_name_etap}.h"])
        .set_constant({"ECMS" => [:double, 4.178]})   # nominal energy of the 8-point scan
        .note(:soft_transition_photon_selection,
              "the soft transition photon from the Ds* is chosen by minimum |DeltaE| inside the " \
              "tagging (DTagAlg reconstruction detail, not tunable in the DSL)")
        .with_decay_card(decay_card_etap)

# --- Ds- tag side: same 14 hadronic tag modes ---
alg_etap.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,              # K+K-pi-
          :DstoKKPiPi0,           # K+K-pi-pi0
          :DstoKsK,               # Ks K-
          :DstoKsKPi0,            # Ks K- pi0
          :DstoKsKsPi,            # Ks Ks pi-
          :DstoKsKPiPi,           # Ks K- pi+ pi-
          :DstoPiPiPi,            # pi+ pi- pi-
          :DstoEtaGGPi,           # eta(gamma gamma) pi-
          :DstoEtaPiPiPi0Pi,      # eta(pi+ pi- pi0) pi-
          :DstoEtaPrimeGGPiPiPi,  # eta'(gamma gamma) pi pi pi
          :DstoEtaPrimeGamPiPiPi, # eta'(gamma pi+ pi-) pi
          :DstoEtaGGRho,          # eta(gamma gamma) rho-
          :DstoEtaPiPiPi0Rho      # eta(pi+ pi- pi0) rho-
  t.charm -1   # tagged side pinned to Ds- (tag mBC / DeltaE stored, not cut here)
end

# --- Signal side: gamma pi+ pi- from eta', plus mu+ and a missing neutrino ---
alg_etap.signal_side do |s|
  s.photons 1                                        # eta' -> gamma pi+ pi-
  s.charged(pip: 1, pim: 1, mup: 1, at_least: true)  # at least one mu+ (extra-track veto is ROOT-level)
  s.missing :nu_mu                                   # missing neutrino
  s.min_photon_energy 0.025                          # photon energy above 25 MeV
end

# --- 3C kinematic fit: four-momentum conservation + M(gamma pi+ pi-) = m(eta');
#     the tighter chi2 < 30 suppresses the non-DsDs* background ---
alg_etap.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :pip, :pim).constrain_to_nominal_mass_of(:etap)  # 1C -> 3C
  f.chi2_cut 30
end

alg_etap.apply
alg_etap.execute_on(data_points + incMC_points + exMC_etap)