# ============================================================================
# BOSS DSL — e+e- -> eta J/psi, J/psi -> l+l- (l = e, mu)
# 44-point BESIII energy scan (3.808 - 4.951 GeV, 22.42 fb^-1),
# represented here by the psi(4260) dataset (703_4260).
# Two independent signal modes (eta -> g g  and  eta -> pi0 pi+ pi-).
# ============================================================================

### Dataset description ###
data_4260  = DatasetManager.real_data.find("703_4260")       # real data (representative scan point)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")   # matching inclusive MC

# Signal decay card — Mode I: eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  0.5 e+  e-  PHOTOS VLL;
  0.5 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal decay card — Mode II: eta -> pi0 pi+ pi-, pi0 -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0 pi0 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  0.5 e+  e-  PHOTOS VLL;
  0.5 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive MC for each of the two eta decay modes
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etajpsi_eta2gg"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etajpsi_eta2pi0pipim"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================================
### Event selection — Mode I: e+e- -> eta J/psi, eta -> gamma gamma ###
# ============================================================================
alg_modeI = Algorithm.new("EtaJpsiGG")
alg_modeI.set_header(["EtaJpsiGGAlg/EtaJpsiGG.h"])
         .set_constant({"ECMS" => [:double, 4.26]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                 # exactly two charged tracks, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                # >=2 photons, timing 0-14, >20 deg to tracks, 25/50 MeV
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {    # high-momentum leptons: one l+ and one l-
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.8
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {   # 4C fit to gamma gamma l+ l-
    nominal
    constrain_four_momentum
    chi2_cut 200   # loose in BOSS; tight 4C chi2 < 40 applied in ROOT
  }

# BOSS-side details not expressible in DSL syntax
alg_modeI.note(:pid_correction_method,
               "lepton ID: track with p > 1.0 GeV/c treated as lepton; electron if E/p > 0.8 " \
               "or EMC energy > 0.8 GeV; muon if EMC energy <= 0.4 GeV")
         .note(:isr_correction,
               "Born cross section extracted over the 44 scan points with initial-state-radiation " \
               "and vacuum-polarization corrections applied iteratively")
         .note(:background_veto,
               "J/psi signal window M(l+l-) in [3.067, 3.127] GeV/c^2 with sidebands " \
               "[3.027, 3.057] and [3.137, 3.167] GeV/c^2; photon E > 0.08 GeV after fit " \
               "-> applied at ROOT level post 4C fit")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_4260, incMC_4260, exMC_modeI])

# ============================================================================
### Event selection — Mode II: e+e- -> eta J/psi, eta -> pi0 pi+ pi- ###
# ============================================================================
alg_modeII = Algorithm.new("EtaJpsiPi0PiPi")
alg_modeII.set_header(["EtaJpsiPi0PiPiAlg/EtaJpsiPi0PiPi.h"])
          .set_constant({"ECMS" => [:double, 4.26]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                 # exactly four charged tracks, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {    # one l+, one l-, one pi+, one pi-
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
    nlp  "==1"
    nlm  "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct pi0 from photon pairs (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :pip, :pim, :lp, :lm]) { # 5C: 4C + M(pi0 pi+ pi-) -> m_eta
    nominal
    constrain_four_momentum
    invariant_mass_of(:pi0, :pip, :pim).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200   # loose in BOSS; tight 5C chi2 < 80 applied in ROOT
  }

alg_modeII.note(:pid_correction_method,
                "lepton ID: track with p > 1.0 GeV/c treated as lepton; electron if E/p > 0.8 " \
                "or EMC energy > 0.8 GeV; muon if EMC energy <= 0.4 GeV")
          .note(:isr_correction,
                "Born cross section extracted over the 44 scan points with initial-state-radiation " \
                "and vacuum-polarization corrections applied iteratively")
          .note(:background_veto,
                "J/psi signal window M(l+l-) in [3.067, 3.127] GeV/c^2 with sidebands " \
                "[3.027, 3.057] and [3.137, 3.167] GeV/c^2 -> applied at ROOT level post 5C fit")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_4260, incMC_4260, exMC_modeII])