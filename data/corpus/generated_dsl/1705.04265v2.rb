### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data (2.93 fb^-1) at sqrt(s)=3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample at 3.773 GeV

# Decay card for the ISR mu+mu- channel: psi(4260) -> mu+ mu- gamma (phase space, photon unreconstructed)
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 mu+ mu- gamma PHSP;
    Enddecay
    End
DECAYCARD

# Decay card for the ISR e+e- channel: psi(4260) -> e+ e- gamma (phase space, photon unreconstructed)
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 e+ e- gamma PHSP;
    Enddecay
    End
DECAYCARD

# 100k-event exclusive MC for each of the two ISR modes
exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_isr_mumu_gamma"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_isr_ee_gamma"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ===== Channel I: e+e- -> gamma_ISR mu+ mu- =====
alg_name_mumu = "ISRGammaMuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})          # sqrt(s) = 3.773 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_mumu = Selection.new
sel_mumu.select_track {                                      # charged-track quality + multiplicity
      cos_theta 0.921                                        # |cos(theta)| < 0.921
      Vz        10.0                                         # |Vz| < 10 cm
      Vr        1.0                                          # Vr < 1 cm
      nChrp    "==1"                                         # exactly one positive track
      nChrn    "==1"                                         # exactly one negative track
      nNet     "==0"                                         # net charge zero
    }
    .remove(:chrgp) { condition "pt_of(:chrgp) < 0.3" }      # reject spiralling tracks with pT < 0.3 GeV/c
    .remove(:chrgn) { condition "pt_of(:chrgn) < 0.3" }
    .pid(method: :probability) {                             # probability-based PID
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6  # p>1.0 GeV/c -> lepton; EMC>0.6 GeV -> e, else mu (P(mu)>P(e))
      identify :pion, against: [:kaon]                       # pi/K separation
      nlp "==1"                                              # exactly one positive lepton
      nlm "==1"                                              # exactly one negative lepton
    }
    .kinematic_fit([:lp, :lm]) {                             # nominal 1C fit: ISR photon is the missing track
      nominal
      miss_track_of :gamma                                   # photon not reconstructed, inferred from 4-momentum conservation
      constrain_four_momentum                                # four-momentum constraint
      chi2_cut 20                                            # chi2 < 20 for mu+mu-gamma
    }

alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)

# ===== Channel II: e+e- -> gamma_ISR e+ e- =====
alg_name_ee = "ISRGammaEE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_ee = Selection.new
sel_ee.select_track {                                        # charged-track quality + multiplicity
      cos_theta 0.921
      Vz        10.0
      Vr        1.0
      nChrp    "==1"
      nChrn    "==1"
      nNet     "==0"
    }
    .remove(:chrgp) { condition "pt_of(:chrgp) < 0.3" }      # reject spiralling tracks with pT < 0.3 GeV/c
    .remove(:chrgn) { condition "pt_of(:chrgn) < 0.3" }
    .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon]                       # pi/K separation
      nlp "==1"
      nlm "==1"
    }
    .remove(:lp) { condition "ep_ratio_of(:lp) < 0.8" }      # electrons: large EMC deposit, E/p > 0.8
    .remove(:lm) { condition "ep_ratio_of(:lm) < 0.8" }
    .kinematic_fit([:lp, :lm]) {                             # nominal 1C fit: ISR photon is the missing track
      nominal
      miss_track_of :gamma
      constrain_four_momentum
      chi2_cut 5                                             # chi2 < 5 for e+e-gamma
    }

alg_ee.with_decay_card(decay_card_ee).apply(sel_ee)

### Execute on datasets ###
root_files_mumu = alg_mumu.execute_on([data_3773, incMC_3773, exMC_mumu])
root_files_ee   = alg_ee.execute_on([data_3773, incMC_3773, exMC_ee])