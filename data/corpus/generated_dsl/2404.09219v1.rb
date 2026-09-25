# ===== Datasets =====
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# ===== Decay cards (EvtGen format) =====
# Signal: D0 -> pi+ pi- eta ; opposite-side tag D0bar -> K+ pi- ; eta -> gamma gamma
decay_card_D0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal: D+ -> pi+ pi0 eta ; opposite-side tag D- -> K+ pi- pi- ; eta -> gamma gamma ; pi0 -> gamma gamma
decay_card_Dp = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ pi0 eta PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ===== Exclusive MC: 200k events per signal mode =====
exMC_D0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0topipieta"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_D0
  config.cross_section   = :default
end

exMC_Dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dptopipi0eta"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_Dp
  config.cross_section   = :default
end

### ===== Tag analysis: D0 -> pi+ pi- eta ===== ###
alg_D0 = TagAnalysis.new("D0ToPiPiEta")
alg_D0.set_header(["D0ToPiPiEtaAlg/D0ToPiPiEta.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: opposite-side D0bar reconstructed in Kpi, Kpipi0, Kpipi0pi0 and Kpipipi
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPi0Pi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: one pi+ and one pi- plus exactly two photons (eta -> gamma gamma)
alg_D0.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.min_photon_angle 10.0
end

# 4C kinematic fit; one gamma-gamma pair constrained to the nominal eta mass
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

# Inexpressible BOSS-side steps captured as notes (applied later in ROOT)
alg_D0.note(:background_veto,
            "veto events with |M(pi+pi-) - m_K_S0| < 0.03 GeV/c^2 to suppress K_S0 -> pi+pi- contamination")
     .note(:bdtg_classifier,
            "BDTG eta classifier using M(gamma gamma), the eta fit chi^2 and the photon helicity angle; retains 83% signal and rejects 78% background")

alg_D0.apply
alg_D0.execute_on([data_3773, incMC_3773, exMC_D0])

### ===== Tag analysis: D+ -> pi+ pi0 eta ===== ###
alg_Dp = TagAnalysis.new("DpToPiPi0Eta")
alg_Dp.set_header(["DpToPiPi0EtaAlg/DpToPiPi0Eta.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: opposite-side D- in Kpipi, Kpipipi0, K_S pi, K_S pi pi0, K_S pi pi pi and KKpi
alg_Dp.tag_side(:Dp) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: one pi+ and 4-48 photons (pi0 -> gamma gamma and eta -> gamma gamma)
alg_Dp.signal_side do |s|
  s.charged(pip: 1)
  s.photons 4..48
  s.min_photon_angle 10.0
end

# 4C kinematic fit; one gamma-gamma pair constrained to the nominal pi0 mass and another to the nominal eta mass
alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

# Inexpressible BOSS-side step captured as a note (applied later in ROOT)
alg_Dp.note(:bdtg_classifier,
            "BDTG eta classifier using M(gamma gamma), the eta fit chi^2 and the photon helicity angle; retains 77% signal and rejects 84% background")

alg_Dp.apply
alg_Dp.execute_on([data_3773, incMC_3773, exMC_Dp])