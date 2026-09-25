# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# 713 R-scan point, sqrt(s) = 2.12655 GeV (nominal 2.125 GeV)
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")     # real data point
incMC_2125 = DatasetManager.inclusive_mc.find("713_Rscan_2125")  # matching inclusive MC

# --- Decay cards (EvtGen) ---
# Both continuum QED channels radiate an ISR/FSR photon, modelled by PHOTOS.
# No intermediate meson is produced, so the BESIII KKMC convention top mother
# psi(4260) is used.

# Bhabha: e+e- -> (gamma) e+e-
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.0  e+  e-  PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# gamma gamma: e+e- -> gamma gamma
decay_card_gg = <<~DECAYCARD
    Decay psi(4260)
    1.0  gamma  gamma  PHSP;
    Enddecay
    End
DECAYCARD

# dimuon: e+e- -> (gamma) mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0  mu+  mu-  PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# generic light-quark continuum background: e+e- -> u ubar
decay_card_uubar = <<~DECAYCARD
    Decay psi(4260)
    1.0  u  anti-u  PHSP;
    Enddecay
    End
DECAYCARD

# --- Exclusive MC samples (1M events each, all at the 2.12655 GeV point) ---
exMC_bhabha = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2125_bhabha"
  config.related_dataset = data_2125
  config.events          = 1_000_000
  config.decay_card      = decay_card_bhabha
  config.cross_section   = :default
end

exMC_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2125_gg"
  config.related_dataset = data_2125
  config.events          = 1_000_000
  config.decay_card      = decay_card_gg
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2125_mumu"
  config.related_dataset = data_2125
  config.events          = 1_000_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

exMC_uubar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2125_uubar"
  config.related_dataset = data_2125
  config.events          = 1_000_000
  config.decay_card      = decay_card_uubar
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Three independent final states -> three separate Algorithm objects (Rule T1).

# ---------- Channel 1: Bhabha  e+e- -> (gamma) e+e- ----------
alg_bhabha = Algorithm.new("BhabhaLumi2125")
alg_bhabha.set_header(["BhabhaLumi2125Alg/BhabhaLumi2125.h"])
          .set_constant({"ECMS" => [:double, 2.12655]})   # sqrt(s) = 2.12655 GeV
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:bhabha_emc_energy, "each charged track required to have EMC energy
            > 0.69 GeV (0.65 * E_beam); electron shower-energy requirement, not
            expressible with the current DSL track-selection keys")

bhabha_selection = Selection.new
bhabha_selection
  .select_track {
    cos_theta 0.8     # |cos(theta)| < 0.8
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # exactly one positive track
    nChrn     "==1"   # exactly one negative track
    nNet      "==0"   # net charge zero
  }
  .assign({:chrgp => :ep, :chrgn => :em})   # Bhabha tracks taken as e+/e-
  .kinematic_fit([:ep, :em]) {              # 4C fit of the e+e- pair
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_bhabha.with_decay_card(decay_card_bhabha).apply(bhabha_selection)
alg_bhabha.execute_on([data_2125, incMC_2125, exMC_bhabha, exMC_uubar])

# ---------- Channel 2: gamma gamma  e+e- -> gamma gamma ----------
alg_gg = Algorithm.new("GG2125")
alg_gg.set_header(["GG2125Alg/GG2125.h"])
      .set_constant({"ECMS" => [:double, 2.12655]})
      .note(:gg_cluster_cos_theta, "each photon cluster required to have
        |cos(theta)| < 0.8; no polar-angle key exists in select_photon")

gg_selection = Selection.new
gg_selection
  .select_track {
    nChrp "==0"   # no charged tracks
    nChrn "==0"
    nNet  "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025   # 25 MeV (barrel)
    energyThreshold_e 0.050   # 50 MeV (endcap)
    nGam              "==2"   # exactly two photons
  }
  .kinematic_fit([:gamma, :gamma]) {   # 4C fit of the two photons
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_gg.with_decay_card(decay_card_gg).apply(gg_selection)
alg_gg.execute_on([data_2125, incMC_2125, exMC_gg, exMC_uubar])

# ---------- Channel 3: dimuon  e+e- -> (gamma) mu+ mu- ----------
alg_mumu = Algorithm.new("MuMu2125")
alg_mumu.set_header(["MuMu2125Alg/MuMu2125.h"])
        .set_constant({"ECMS" => [:double, 2.12655]})
        .note(:background_veto, "charged tracks with E/p > 0.4 vetoed (before PID)
          to suppress Bhabha e+e- contamination in the dimuon selection")

mumu_selection = Selection.new
mumu_selection
  .select_track {
    cos_theta 0.8     # |cos(theta)| < 0.8
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # exactly one positive track
    nChrn     "==1"   # exactly one negative track
    nNet      "==0"   # net charge zero
  }
  .pid(method: :probability) {   # probability PID; high-momentum tracks treated as leptons
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"   # one positive lepton (mu+ after the E/p <= 0.4 electron veto)
    nlm "==1"   # one negative lepton (mu-)
  }
  .kinematic_fit([:lp, :lm]) {   # 4C fit of the mu+ mu- pair
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mumu.with_decay_card(decay_card_mumu).apply(mumu_selection)
alg_mumu.execute_on([data_2125, incMC_2125, exMC_mumu, exMC_uubar])