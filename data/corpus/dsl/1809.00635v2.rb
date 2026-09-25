# 1809.00635v2: J/psi → eta' e+ e- EM Dalitz decay and dark photon search at BESIII
# Two eta' decay modes: I (gamma pi+ pi-) + II (pi+ pi- eta(gammagamma))
# Separate Algorithm per mode per Rule T1

# === Datasets ===
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# =========================================================
# Mode I: J/psi → eta' e+ e-, eta' → gamma pi+ pi-
# Final state: gamma pi+ pi- e+ e-
# =========================================================
decay_card_I = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' e+ e- PHSP;
  Enddecay
  Decay eta'
  1.0000 gamma pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_etap_ee_gam_pipi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_I
  config.cross_section   = :default
end

alg_I = Algorithm.new("JpsietapEE_ModeI")
alg_I.set_header(["JpsietapEE_ModeIAlg/JpsietapEE_ModeI.h"])
     .set_constant({ "ECMS" => [:double, 3.097] })

sel_I = Selection.new

# Charged track selection: 4 tracks, net zero
sel_I.select_track do
  cos_theta   0.93
  Vz          100.0
  Vr          10.0
  nChrp       "==2"
  nChrn       "==2"
  nNet        "==0"
end

# Photon selection: ≥1 photon (eta' → gamma pi+ pi-)
sel_I.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=1"
end

# PID: high-momentum leptons for e+/e-; no PID requirement on pions
sel_I.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
end

# Assign remaining charged tracks as pions
sel_I.assign({ chrgp: :pip, chrgn: :pim })

# 4C kinematic fit with vertex constraint for the 4 charged tracks
# chi2 < 100; best combination by minimum chi2
sel_I.kinematic_fit([:gamma, :pip, :pim, :ep, :em]) do
  nominal
  constrain_four_momentum
  vertex_fit([1, 2, 3, 4])
  chi2_cut 100
end

alg_I.with_decay_card(decay_card_I)
     .note(:gamma_conversion_veto,
       "Veto γ→e+e- conversion: δxy = sqrt(Rx²+Ry²) < 2 cm; " \
       "rejects events with e+e- conversion vertex at beam pipe or MDC inner wall")
     .note(:pi0_dalitz_veto,
       "Mode I only: veto M(γ e+e-) ∈ [0.10, 0.16] GeV/c² to reject π0 Dalitz background; " \
       "applied in ROOT")
     .note(:bf_measurement,
       "Unbinned ML fit to M(γπ+π-) with signal MC shape convolved with Gaussian; " \
       "peaking backgrounds fixed from MC; non-peaking by 2nd-order Chebyshev")
     .note(:dark_photon_search,
       "Dark photon search in M(e+e-) with η' mass window [0.93, 0.98] GeV/c²; " \
       "excluding ω [0.74,0.84] and φ [1.00,1.04] mass regions; " \
       "γ conversion veto removed for dark photon search")
     .note(:tracking_correction,
       "Electron tracking efficiency corrected by factor δ=1.012 per e±; " \
       "2D correction in pT vs cosθ from radiative Bhabha control sample")
     .apply(sel_I)

alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_I])

# =========================================================
# Mode II: J/psi → eta' e+ e-, eta' → pi+ pi- eta, eta → gamma gamma
# Final state: gamma gamma pi+ pi- e+ e-
# =========================================================
decay_card_II = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' e+ e- PHSP;
  Enddecay
  Decay eta'
  1.0000 pi+ pi- eta PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_etap_ee_pipi_eta"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_II
  config.cross_section   = :default
end

alg_II = Algorithm.new("JpsietapEE_ModeII")
alg_II.set_header(["JpsietapEE_ModeIIAlg/JpsietapEE_ModeII.h"])
      .set_constant({ "ECMS" => [:double, 3.097] })

sel_II = Selection.new

# Charged track selection: 4 tracks, net zero
sel_II.select_track do
  cos_theta   0.93
  Vz          100.0
  Vr          10.0
  nChrp       "==2"
  nChrn       "==2"
  nNet        "==0"
end

# Photon selection: ≥2 photons (eta → gamma gamma)
sel_II.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

# PID: high-momentum leptons for e+/e-; no PID requirement on pions
sel_II.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
end

# Assign remaining charged tracks as pions
sel_II.assign({ chrgp: :pip, chrgn: :pim })

# Reconstruct eta → gamma gamma with Kalman mass constraint
# Pre-fit: M(γγ) ∈ [0.48, 0.60] GeV/c²
sel_II.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# 4C kinematic fit with vertex constraint for the 4 charged tracks
sel_II.kinematic_fit([:eta, :pip, :pim, :ep, :em]) do
  nominal
  constrain_four_momentum
  vertex_fit([1, 2, 3, 4])
  chi2_cut 100
end

alg_II.with_decay_card(decay_card_II)
      .note(:gamma_conversion_veto,
        "Veto γ→e+e- conversion: δxy = sqrt(Rx²+Ry²) < 2 cm; " \
        "rejects events with e+e- conversion vertex at beam pipe or MDC inner wall")
      .note(:eta_mass_window,
        "η → γγ: M(γγ) ∈ [0.48, 0.60] GeV/c² pre-fit constraint")
      .note(:bf_measurement,
        "Unbinned ML fit to M(γγπ+π-) with signal MC shape convolved with Gaussian; " \
        "peaking backgrounds fixed from MC; non-peaking by 1st-order Chebyshev")
      .note(:dark_photon_search,
        "Dark photon search in M(e+e-) with η' mass window [0.93, 0.98] GeV/c²; " \
        "excluding ω [0.74,0.84] and φ [1.00,1.04] mass regions")
      .apply(sel_II)

alg_II.execute_on([jpsi_data, jpsi_incMC, exMC_II])