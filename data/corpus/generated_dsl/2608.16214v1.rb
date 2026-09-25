# =====================================================================
# BESIII analysis: e+e- -> J/psi / psi(2S) -> Sigma0 anti-Sigma0 eta
#   Sigma0      -> gamma Lambda0,        anti-Sigma0 -> gamma anti-Lambda0
#   Lambda0     -> p+ pi-,               anti-Lambda0 -> anti-p- pi+
#   eta         -> gamma gamma
# BOSS part: datasets + exclusive MC + event selection up to the 4C fit.
# ROOT-level items (pi0 veto, Sigma0/eta joint-chi2 selection, J/psi veto,
# tightened chi2 cuts) are intentionally not expressed here.
# =====================================================================

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # psi(3686) inclusive MC

### Decay cards (EvtGen format) ###
decay_card_jpsi = <<~DECAYCARD
  Decay J/psi
  1.000 Sigma0 anti-Sigma0 eta PHSP;
  Enddecay

  Decay Sigma0
  1.000 gamma Lambda0 PHSP;
  Enddecay

  Decay anti-Sigma0
  1.000 gamma anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_psip = <<~DECAYCARD
  Decay psi(2S)
  1.000 Sigma0 anti-Sigma0 eta PHSP;
  Enddecay

  Decay Sigma0
  1.000 gamma Lambda0 PHSP;
  Enddecay

  Decay anti-Sigma0
  1.000 gamma anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC: 200k events for each resonance, same decay chain ###
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Sigma0Sigma0barEta"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_Sigma0Sigma0barEta"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common selection chain shared by both resonances
event_selection_common = Selection.new
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        20.0    # |Vz| < 20 cm
    Vr        10.0    # Vr < 10 cm
    nChrp     ">=2"   # at least two positive tracks
    nChrn     ">=2"   # at least two negative tracks
  }
  .select_photon {
    tdc_emc_start     0      # EMC TDC lower bound
    tdc_emc_end       14     # EMC TDC upper bound
    energyThreshold_b 0.025  # 25 MeV barrel energy threshold
    energyThreshold_e 0.050  # 50 MeV endcap energy threshold
    nGam              ">=4"  # at least four good photons
  }
  .pid(method: :probability) {
    prob_cut 0.001                                 # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]      # p / p-bar separated from K, pi
    identify :pion,   against: [:kaon, :proton]    # pi+ / pi- separated from K, p
    nprp ">=1"                                     # at least one proton
    nprm ">=1"                                     # at least one anti-proton
    npip ">=1"                                     # at least one pi+
    npim ">=1"                                     # at least one pi-
  }
  # Secondary vertex fit: Lambda -> p pi-
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Secondary vertex fit: anti-Lambda -> anti-p pi+
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit to Lambda anti-Lambda gamma gamma gamma gamma
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200   # loose BOSS-level cut (tightened later in ROOT: 60 J/psi, 40 psi(2S))
  }
  # Alternative 4C fit to Lambda anti-Lambda gamma gamma gamma:
  # competing one-fewer-photon hypothesis; chi2 stored for a ROOT-level veto,
  # no chi2_cut and no nominal here.
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

### Algorithms: one per resonance (different decay-card mother),
### sharing the same selection chain ###
alg_jpsi = Algorithm.new("Sigma0Sigma0barEtaJpsi")
alg_jpsi.set_header(["Sigma0Sigma0barEtaJpsiAlg/Sigma0Sigma0barEtaJpsi.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})
alg_jpsi.with_decay_card(decay_card_jpsi).apply(event_selection_common.dup)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

alg_psip = Algorithm.new("Sigma0Sigma0barEtaPsip")
alg_psip.set_header(["Sigma0Sigma0barEtaPsipAlg/Sigma0Sigma0barEtaPsip.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})
alg_psip.with_decay_card(decay_card_psip).apply(event_selection_common.dup)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])