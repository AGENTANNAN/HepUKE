# ============================================================
# Paper: Measurement of the Branching Fraction of and Search for
#        a CP-Violating Asymmetry in eta' -> pi+pi-e+e- at BESIII
# arXiv: 2011.07902v2
# ============================================================

###
### Dataset: J/psi data (BOSS 708, 3097 MeV)
###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Signal Channel: J/psi -> gamma eta', eta' -> pi+ pi- e+ e-
# ============================================================
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta'  HELAMP 0 0 1 0 1 0 0 0 1 0 1 0 0 0;
    Enddecay

    Decay eta'
    1.000 pi+ pi- e+ e-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_gamma_etap_pipi_ee"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

alg_signal = Algorithm.new("JpsiEtaPipPimEE")
alg_signal.set_header(["JpsiEtaPipPimEEAlg/JpsiEtaPipPimEE.h"])
alg_signal.set_constant({"ECMS" => [:double, 3.097]})

sel_signal = Selection.new

# Charged track selection: 4 tracks (e+, e-, pi+, pi-), net charge zero
sel_signal.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       "==2"
  nChrn       "==2"
  nNet        "==0"
end

# Photon selection: at least 1 photon (radiative photon from J/psi)
# Barrel: 25 MeV, Endcap: 50 MeV, angle to track > 15°, timing within 700 ns
sel_signal.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    15.0
  nGam              ">=1"
end

# PID: chi2-sum combinatorial PID — assign each track as pion or electron
# The combined chi2 = chi2_4C + sum chi2_PID(j) is minimized
sel_signal.pid(method: :chi2_sum) do
  chi_min_cut 4
  identify :pion, :electron
end

# Nominal 4C kinematic fit under gamma pi+ pi- e+ e- hypothesis
sel_signal.kinematic_fit([:gamma, :pip, :pim, :ep, :em]) do
  nominal
  constrain_four_momentum
  chi2_cut 62
end

# Competing hypothesis 1: gamma pi+ pi- mu+ mu-
# Reassign electron tracks to muon, refit with same track indices
sel_signal.assign({ep: :mup, em: :mum})
          .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) do
            use_track_index_from_nominal_kmfit
            constrain_four_momentum
          end

# Competing hypothesis 2: gamma pi+ pi- pi+ pi-
# Reassign electron tracks to pion, refit with same track indices
sel_signal.assign({ep: :pip, em: :pim})
          .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) do
            use_track_index_from_nominal_kmfit
            constrain_four_momentum
          end

alg_signal.note(:photon_conversion_veto, "Photon conversion veto using R_xy (ee-vertex distance from IP), M_ee^BP (ee invariant mass at beam pipe), and Phi_ee (z-projection of ee opening angle); veto windows defined by polygons in M_ee^BP vs R_xy and Phi_ee vs R_xy planes")
alg_signal.note(:helix_correction, "Helix parameter correction applied for photon conversion veto: ee track helices re-referenced to beam pipe intersection point for M_ee^BP calculation")
alg_signal.note(:best_hypothesis_selection, "Events kept where gamma pi+pi-e+e- hypothesis has smallest combined chi2_4C+PID among the three tested hypotheses (ee, mumu, 4pi)")
alg_signal.note(:etap_mass_window, "Signal region: |M(pi+pi-e+e-) - m_eta'| < 0.02 GeV/c2 applied in ROOT")
alg_signal.note(:bin_migration, "sin(2phi) bin migration correction (alpha = 21.3%) applied in ROOT-level asymmetry calculation")

alg_signal.with_decay_card(decay_card_signal)
          .apply(sel_signal)

# ============================================================
# Normalization Channel: J/psi -> gamma eta', eta' -> pi+ pi- gamma
# ============================================================
decay_card_norm = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta'  HELAMP 0 0 1 0 1 0 0 0 1 0 1 0 0 0;
    Enddecay

    Decay eta'
    1.000 pi+ pi- gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_norm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_gamma_etap_pipi_gamma"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_norm
  config.cross_section   = :default
end

alg_norm = Algorithm.new("JpsiEtaPipPimGam")
alg_norm.set_header(["JpsiEtaPipPimGamAlg/JpsiEtaPipPimGam.h"])
alg_norm.set_constant({"ECMS" => [:double, 3.097]})

sel_norm = Selection.new

# Charged track selection: 2 tracks (pi+, pi-), net charge zero
sel_norm.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       "==1"
  nChrn       "==1"
  nNet        "==0"
end

# Photon selection: at least 2 photons (radiative from J/psi + from eta' decay)
sel_norm.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    15.0
  nGam              ">=2"
end

# PID: chi2-sum combinatorial PID — assign each track as pion
sel_norm.pid(method: :chi2_sum) do
  chi_min_cut 4
  identify :pion
end

# 4C kinematic fit under gamma pi+ pi- gamma hypothesis
sel_norm.kinematic_fit([:gamma, :pip, :pim, :gamma]) do
  nominal
  constrain_four_momentum
  chi2_cut 140
end

alg_norm.note(:radiative_photon_energy, "Radiative photon energy: |E_gamma - 1.4 GeV| < 20 MeV (J/psi -> gamma eta' monoenergetic photon); candidate with closest energy selected")
alg_norm.note(:etap_decay_photon_energy, "Photon from eta' decay: E_gamma > 0.15 GeV to suppress pi0->gamma gamma and eta->gamma gamma backgrounds")
alg_norm.note(:pi0_eta_veto, "pi0->gamma gamma and eta->gamma gamma contamination suppressed by E_gamma > 0.15 GeV requirement")

alg_norm.with_decay_card(decay_card_norm)
         .apply(sel_norm)

# ============================================================
# Execute both algorithms
# ============================================================
alg_signal.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
alg_norm.execute_on([jpsi_data, jpsi_incMC, exMC_norm])