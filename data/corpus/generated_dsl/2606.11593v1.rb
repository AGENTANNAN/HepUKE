# frozen_string_literal: true
# =============================================================================
#  e+ e- -> Ds*+ Ds-     (4.130 - 4.230 GeV energy scan, total 7.33 fb^-1)
#    Ds*+ -> gamma Ds+ (93.5%) | pi0 Ds+ (6.5%)
#    Ds+  -> tau+ nu_tau
#    tau+ -> e+ nu_e nu_tau | mu+ nu_mu nu_tau | pi+ nu_tau | pi+ pi0 nu_tau
#  Tag side : Ds- reconstructed in eleven hadronic modes (DTag, pre-stored tags)
# =============================================================================

### ------------------------- real data (energy scan) -------------------------
data_4130 = DatasetManager.real_data.find("705_4130")   # 4128.5 MeV
data_4160 = DatasetManager.real_data.find("705_4160")   # 4157.4 MeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4178.0 MeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4188.8 MeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4198.9 MeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4209.2 MeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4218.7 MeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4226.3 MeV
scan_points = [data_4130, data_4160, data_4180, data_4190,
               data_4200, data_4210, data_4220, data_4230]

### ---------------------------- inclusive MC --------------------------------
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

### ------------------------- decay card (EvtGen) ----------------------------
# Signal process: e+e- -> Ds*+ Ds-, the Ds*+ transition photon/pi0, Ds+ -> tau+ nu,
# and the four tau+ decay sub-channels. The tag side (Ds-) is reconstructed by
# DTagAlg from the event and is therefore not forced here.
decay_card_signal = <<~DECAYCARD
    Decay psi(4160)
    1.0000  Ds*+  Ds-        PHSP;
    Enddecay

    Decay Ds*+
    0.935   gamma  Ds+       PHSP;
    0.065   pi0    Ds+       PHSP;
    Enddecay

    Decay Ds+
    1.000   tau+  nu_tau     VLL;
    Enddecay

    Decay tau+
    0.25    e+   nu_e   anti-nu_tau     PHOTOS VLL;
    0.25    mu+  nu_mu  anti-nu_tau     PHOTOS VLL;
    0.25    pi+  anti-nu_tau            VLL;
    0.25    pi+  pi0    anti-nu_tau     PHSP;
    Enddecay

    Decay pi0
    1.000   gamma  gamma     PHSP;
    Enddecay

    End
DECAYCARD

### --------------- exclusive signal MC (500k events per point) ---------------
# One signal MC sample per energy point; the generation energy is fixed to 4.178 GeV.
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "DsStarTauNu_signal"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### ------------------- common dataset list for execution --------------------
all_datasets = [data_4130, incMC_4130, data_4160, incMC_4160,
                data_4180, incMC_4180, data_4190, incMC_4190,
                data_4200, incMC_4200, data_4210, incMC_4210,
                data_4220, incMC_4220, data_4230, incMC_4230] + exMCs_signal

### ------------------- eleven hadronic Ds- tag modes ------------------------
ds_tag_modes = [:DstoKsK,                    # K_S0 K-
                :DstoKKPi,                   # K+ K- pi-
                :DstoKKPiPi0,                # K+ K- pi- pi0
                :DstoKsKPiPi,                # K_S0 K- pi+ pi-
                :DstoKsKPimPim,              # K_S0 K+ pi- pi-
                :DstoPiPiPi,                 # pi+ pi- pi-
                :DstoPiEta,                  # pi- eta
                :DstoPiPi0Eta,               # pi- pi0 eta
                :DstoPiEtaPrimeEtaPiPi,      # pi- eta'(-> pi+ pi- eta)
                :DstoPiEtaPrimeGammaRho,     # pi- eta'(-> gamma rho0)
                :DstoKPiPi]                  # K- pi+ pi-

# =============================================================================
#  Channel 1 :  tau+ -> e+ nu_e anti-nu_tau
# =============================================================================
alg_e = TagAnalysis.new("DsStarTauNuE")
alg_e.set_header(["DsStarTauNuEAlg/DsStarTauNuE.h"])
     .set_constant({ "ECMS" => [:double, 4.178] })
     .with_decay_card(decay_card_signal)

alg_e.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)   # eleven hadronic Ds- modes -> single tag
  t.charm(-1)              # pin the tagged side to Ds-
end

alg_e.signal_side do |s|
  s.photons 1..48          # Ds*+ transition photon / pi0 gammas + soft extra showers
  s.min_photon_energy 0.025
  s.min_photon_angle  10.0
  s.charged(ep: 1)         # exactly one signal charged track, identified as e+
  s.missing :nu            # massless missing system (nu_e + nu_tau)
end

alg_e.fit do |f|
  f.constrain_four_momentum  # tag + signal + missing = measured CMS four-momentum
  f.chi2_cut 200
end

alg_e.note(:efficiency_curve,
      "tag-side track quality |Vz|<10 cm, |Vxy|<1 cm, |cos(theta)|<0.93 and all " \
      "K_S0 / pi0 / eta / rho0 / eta' reconstruction windows are handled inside " \
      "DTagAlg; tag Ds- recoil mass [2.050,2.195] GeV/c^2 and the mode-dependent " \
      "+-3 sigma invariant-mass window are stored and cut in ROOT.")
alg_e.note(:pid_correction_method,
      "signal-side e+ identified from MDC+TOF likelihoods: L(e)>0.001, " \
      "L(e)/(L(e)+L(pi)+L(K))>0.8, p>0.2 GeV/c and EMC/p>0.8.")
alg_e.note(:background_veto,
      "the K- pi+ pi- tag mode vetoes M(pi+ pi-) in [0.480,0.515] GeV/c^2; " \
      "DeltaE = Ecm - E_ST - E'_miss - E_gamma(pi0) required in [-0.2,0.2] GeV " \
      "(ROOT-level cut on stored quantities).")

alg_e.apply                     # takes NO Selection argument
alg_e.execute_on(all_datasets)

# =============================================================================
#  Channel 2 :  tau+ -> mu+ nu_mu anti-nu_tau
# =============================================================================
alg_mu = TagAnalysis.new("DsStarTauNuMu")
alg_mu.set_header(["DsStarTauNuMuAlg/DsStarTauNuMu.h"])
      .set_constant({ "ECMS" => [:double, 4.178] })
      .with_decay_card(decay_card_signal)

alg_mu.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm(-1)
end

alg_mu.signal_side do |s|
  s.photons 1..48
  s.min_photon_energy 0.025
  s.min_photon_angle  10.0
  s.charged(mup: 1)        # exactly one signal charged track, identified as mu+
  s.missing :nu
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu.note(:efficiency_curve,
      "tag-side K_S0 / pi0 / eta / rho0 / eta' windows and the tag recoil-mass / " \
      "+-3 sigma Ds- mass windows are applied as stored-variable ROOT cuts.")
alg_mu.note(:pid_correction_method,
      "signal-side mu+ identified from MDC+TOF likelihoods: L(mu)>0.001 and greater " \
      "than L(K) and L(e), EMC energy 0.1-0.3 GeV, p>0.5 GeV/c and MUC hit-depth cuts.")
alg_mu.note(:background_veto,
      "K- pi+ pi- tag mode vetoes M(pi+ pi-) in [0.480,0.515] GeV/c^2; " \
      "DeltaE in [-0.2,0.2] GeV applied in ROOT.")

alg_mu.apply
alg_mu.execute_on(all_datasets)

# =============================================================================
#  Channel 3 :  tau+ -> pi+ anti-nu_tau
# =============================================================================
alg_pi = TagAnalysis.new("DsStarTauNuPi")
alg_pi.set_header(["DsStarTauNuPiAlg/DsStarTauNuPi.h"])
      .set_constant({ "ECMS" => [:double, 4.178] })
      .with_decay_card(decay_card_signal)

alg_pi.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm(-1)
end

alg_pi.signal_side do |s|
  s.photons 1..48
  s.min_photon_energy 0.025
  s.min_photon_angle  10.0
  s.charged(pip: 1)        # exactly one signal charged track, identified as pi+
  s.missing :nu
end

alg_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_pi.note(:efficiency_curve,
      "tag-side hadron reconstruction windows plus the tag recoil-mass window " \
      "[2.050,2.195] GeV/c^2 and the +-3 sigma Ds- mass window are stored and cut in ROOT.")
alg_pi.note(:pid_correction_method,
      "signal-side pi+ identified from MDC+TOF likelihoods: L(pi)>L(K), EMC/p<0.9, " \
      "no pi0, residual photon energy <0.3 GeV and |cos(theta_miss)|<0.9.")
alg_pi.note(:background_veto,
      "K- pi+ pi- tag mode vetoes M(pi+ pi-) in [0.480,0.515] GeV/c^2; " \
      "DeltaE in [-0.2,0.2] GeV applied in ROOT.")

alg_pi.apply
alg_pi.execute_on(all_datasets)

# =============================================================================
#  Channel 4 :  tau+ -> pi+ pi0 anti-nu_tau
# =============================================================================
alg_pipi0 = TagAnalysis.new("DsStarTauNuPiPi0")
alg_pipi0.set_header(["DsStarTauNuPiPi0Alg/DsStarTauNuPiPi0.h"])
         .set_constant({ "ECMS" => [:double, 4.178] })
         .with_decay_card(decay_card_signal)

alg_pipi0.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm(-1)
end

alg_pipi0.signal_side do |s|
  s.photons 1..48          # transition gamma + the pi0 -> gamma gamma pair
  s.min_photon_energy 0.025
  s.min_photon_angle  10.0
  s.charged(pip: 1)        # the pi+ from tau+ -> pi+ pi0 nu
  s.missing :nu
end

alg_pipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # gamma gamma -> pi0 (1C)
  f.chi2_cut 200
end

alg_pipi0.note(:efficiency_curve,
      "tag-side reconstruction windows, the tag recoil-mass window [2.050,2.195] GeV/c^2 " \
      "and the +-3 sigma Ds- mass window are stored and cut in ROOT.")
alg_pipi0.note(:pid_correction_method,
      "signal-side pi+ identified from MDC+TOF likelihoods: L(pi)>L(K), EMC/p<0.9, " \
      "no extra pi0 and residual photon energy <0.3 GeV.")
alg_pipi0.note(:background_veto,
      "tau+ -> pi+ pi0 mode: |M(pi+ pi0) - m_rho+| < 0.2 GeV/c^2, residual photon " \
      "energy < 0.1 GeV and missing-mass squared in [3.82,3.98] GeV^2/c^4; " \
      "DeltaE in [-0.2,0.2] GeV - all applied as ROOT cuts on stored quantities.")

alg_pipi0.apply
alg_pipi0.execute_on(all_datasets)

# -----------------------------------------------------------------------------
# Produced ROOT ntuples (tag mBC / DeltaE / missing-mass stored unconditionally)
# for the four tau+ signal channels over the 4.130-4.230 GeV scan points.
# -----------------------------------------------------------------------------