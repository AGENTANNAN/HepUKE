# =============================================================================
#  BOSS-side HepScript DSL
#  Double-tag analysis of  e+e- -> Lambda_c+ Lambda_c- , Lambda_c+ -> p eta'
#  (single tag of Lambda_c- on the opposite side)
#  Seven c.m. energy points (4.5 fb^-1):
#    4.600 / 4.612 / 4.628 / 4.641 / 4.661 / 4.682 / 4.699 GeV
# =============================================================================

### ---------------------------------------------------------------- datasets ---
# Real data and inclusive MC at the seven c.m. energy points
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV
  DatasetManager.real_data.find("706_4610"),   # 4.612 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.628 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.641 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.661 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.682 GeV
  DatasetManager.real_data.find("706_4700"),   # 4.699 GeV
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

### ------------------------------------------------------------- decay cards ---
# Sub-channel A : Lambda_c+ -> p eta' , eta' -> pi+ pi- eta , eta -> gamma gamma
decay_card_etap_eta = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 p+ eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Sub-channel B : Lambda_c+ -> p eta' , eta' -> pi+ pi- gamma
decay_card_etap_gamma = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 p+ eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- gamma PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------------------------------------------ exclusive MC ---
# 500k-event signal MC per sub-channel, one sample per c.m. energy point
exMCs_etap_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_LcToPEtaPrimeToPiPiEta"
  config.events        = 500_000
  config.decay_card    = decay_card_etap_eta
  config.cross_section = :default
end

exMCs_etap_gamma = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_LcToPEtaPrimeToPiPiGamma"
  config.events        = 500_000
  config.decay_card    = decay_card_etap_gamma
  config.cross_section = :default
end

### ============================= sub-channel A =============================== ###
# eta' -> pi+ pi- eta (eta -> gamma gamma) : the eta is treated as missing
alg_etap_eta = TagAnalysis.new("LcToPEtaPrimeToPiPiEta")
alg_etap_eta.set_header(["LcToPEtaPrimeToPiPiEtaAlg/LcToPEtaPrimeToPiPiEta.h"])
            .set_constant({ "ECMS" => [:double, 4.600] })
            .set_alias({ "std::vector<double>" => "Vdouble" })

# Tag side : single tag of Lambda_c- in the hadronic family topologies, charm = -1
alg_etap_eta.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoPKPi,           # family  p K pi
          :LambdacPtoLambdaPi,       # family  Lambda pi
          :LambdacPtoLambdaPiPiPi,   # family  Lambda pi pi pi
          :LambdacPtoPKPiPiPi        # family  p K pi pi pi
  t.charm -1
  t.window :mBC, min: 2.275, max: 2.310   # explicitly requested tag-side M_BC window
end

# Signal side : p, pi+, pi- recoiling against the ST Lambda_c- ; eta is missing
alg_etap_eta.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 1)
  s.missing :eta
  s.require_charge 1
end

# Four-momentum constrained kinematic fit
alg_etap_eta.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap_eta
  .note(:tag_window_deltae_per_mode,
        "Tag-side DeltaE windows are asymmetric and mode-dependent (different bounds for " \
        "each of the hadronic Lambda_c tag topologies); they are applied to the stored tag " \
        "candidates in the ROOT analysis, with DTagAlg's own mBC/DeltaE selector cuts kept " \
        "disabled so that every candidate is stored.")
  .note(:background_veto,
        "Signal-side vetoes: Lambda suppressed by M(p pi-) > 1.125 GeV/c2 and K_S0 suppressed " \
        "by M(pi+ pi-) outside (0.490,0.510) GeV/c2, evaluated on the selected signal tracks.")
  .note(:etap_to_pipim_eta_recoil,
        "eta' -> pi+ pi- eta (eta -> gamma gamma) sub-channel: the eta is treated as missing " \
        "and the recoil mass M_rec(Lambda_c- p pi+ pi-) must lie in the eta window.  The " \
        "result is extracted, together with the eta' -> pi+ pi- gamma sub-channel, from a " \
        "simultaneous unbinned maximum-likelihood fit to M_rec(Lambda_c- p) and " \
        "M(pi+ pi- gamma) at the ROOT level.")

alg_etap_eta.with_decay_card(decay_card_etap_eta)
alg_etap_eta.apply
alg_etap_eta.execute_on(data_points + incMC_points + exMCs_etap_eta)

### ============================= sub-channel B =============================== ###
# eta' -> pi+ pi- gamma : the photon is taken from the showers left unused by the tag
alg_etap_gamma = TagAnalysis.new("LcToPEtaPrimeToPiPiGamma")
alg_etap_gamma.set_header(["LcToPEtaPrimeToPiPiGammaAlg/LcToPEtaPrimeToPiPiGamma.h"])
              .set_constant({ "ECMS" => [:double, 4.600] })
              .set_alias({ "std::vector<double>" => "Vdouble" })

alg_etap_gamma.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoPKPi,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoPKPiPiPi
  t.charm -1
  t.window :mBC, min: 2.275, max: 2.310
end

# Signal side : p, pi+, pi- plus the photons of the signal side
alg_etap_gamma.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 1)
  s.photons 2
  s.require_charge 1
end

alg_etap_gamma.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap_gamma
  .note(:tag_window_deltae_per_mode,
        "Tag-side DeltaE windows are asymmetric and mode-dependent (one set of bounds per " \
        "hadronic Lambda_c tag topology); applied on the stored tag candidates in ROOT with " \
        "DTagAlg's own selector cuts disabled.")
  .note(:background_veto,
        "Signal-side vetoes: Lambda suppressed by M(p pi-) > 1.125 GeV/c2 and K_S0 suppressed " \
        "by M(pi+ pi-) outside (0.490,0.510) GeV/c2.")
  .note(:etap_to_pipim_gamma_deltae,
        "eta' -> pi+ pi- gamma sub-channel: the photon is chosen among the signal-side " \
        "showers left unused by the tag and selected with a DeltaE window; M(pi+ pi- gamma) " \
        "is stored and enters the simultaneous unbinned maximum-likelihood fit.")

alg_etap_gamma.with_decay_card(decay_card_etap_gamma)
alg_etap_gamma.apply
alg_etap_gamma.execute_on(data_points + incMC_points + exMCs_etap_gamma)