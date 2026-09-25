# Electromagnetic Dalitz decays J/psi -> P e+ e- (P = eta', eta, pi0) at sqrt(s) = 3.097 GeV
# Five independent final states + one peaking-background sample.
# BOSS part only: dataset preparation + event selection up to the final kinematic fit.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data (225.3 x 10^6 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC
cont_data  = DatasetManager.real_data.find("712_3773")      # 3.773 GeV continuum data (~2.9 fb^-1)
cont_incMC = DatasetManager.inclusive_mc.find("712_3773")   # Continuum inclusive MC

### Decay cards (EvtGen format) ###

# Mode I: J/psi -> eta' e+ e-, eta' -> gamma pi+ pi-
decay_card_mode1 = <<~DECAYCARD
  Decay J/psi
  1.000 eta' e+ e- PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: J/psi -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_mode2 = <<~DECAYCARD
  Decay J/psi
  1.000 eta' e+ e- PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: J/psi -> eta e+ e-, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_mode3 = <<~DECAYCARD
  Decay J/psi
  1.000 eta e+ e- PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode IV: J/psi -> eta e+ e-, eta -> gamma gamma
decay_card_mode4 = <<~DECAYCARD
  Decay J/psi
  1.000 eta e+ e- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode V: J/psi -> pi0 e+ e-, pi0 -> gamma gamma
decay_card_mode5 = <<~DECAYCARD
  Decay J/psi
  1.000 pi0 e+ e- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background: J/psi -> eta' gamma, eta' -> gamma pi+ pi- (photon conversion)
decay_card_bkg = <<~DECAYCARD
  Decay J/psi
  1.000 eta' gamma PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples ###
# 100k events for each of the five signal modes ...
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_etap_ee_gammapipi"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_etap_ee_pipieta"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_ee_pipipi0"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_ee_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_mode4
  config.cross_section   = :default
end

exMC_mode5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_pi0_ee_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_mode5
  config.cross_section   = :default
end

# ... and 500k events for the peaking J/psi -> eta' gamma conversion background.
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_etapgamma_gammapipi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###

## ================= Mode I: eta' -> gamma pi+ pi- =================
alg1_name = "JpsiEtapEEGammaPiPi"
alg_mode1 = Algorithm.new(alg1_name)
alg_mode1.set_header(["#{alg1_name}Alg/#{alg1_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .note(:photon_conversion_veto, "photon-conversion veto: events in which a photon converts with reconstructed conversion radius sqrt(Rx^2 + Ry^2) < 2 cm are rejected; the conversion radius is not accessible from the DSL photon-selection block and is applied in the BOSS algorithm code")
         .note(:babayaga_qed_background, "continuum QED background e+e- -> e+e-gamma(gamma) and 3gamma modelled with Babayaga on the 3.773 GeV continuum data sample; Babayaga is not an EvtGen decay card, so the sample is generated separately and analysed with the same selection")

sel_mode1 = Selection.new
sel_mode1.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp "==2"
           nChrn "==2"
           nNet  "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=1"                      # eta' -> gamma pi+ pi- : at least one photon
         }
         .pid(method: :probability) {
           prob_cut 0.001                  # PID probability > 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"                       # exactly one e+
           nlm "==1"                       # exactly one e-
         }
         .remove([:lp <= :chrgp, :lm <= :chrgn])       # drop the identified electrons
         .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks taken as pi+/pi- (no further PID)
         .kinematic_fit([:lp, :lm, :gamma, :pip, :pim]) {
           nominal
           vertex_fit([3, 4])              # common vertex for pi+ (3) and pi- (4)
           constrain_four_momentum         # 4C fit to e+e-gamma pi+ pi-
           invariant_mass_of(:gamma, :lp, :lm).out_of(0.10, 0.16)  # veto 0.10 < M(gamma e+e-) < 0.16 GeV/c^2
           chi2_cut 200                    # loose first-pass cut; tight chi2_4C < 100 applied in ROOT
         }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.execute_on([jpsi_data, jpsi_incMC, cont_data, cont_incMC, exMC_mode1, exMC_bkg])

## ================= Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma =================
alg2_name = "JpsiEtapEEPiPiEta"
alg_mode2 = Algorithm.new(alg2_name)
alg_mode2.set_header(["#{alg2_name}Alg/#{alg2_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .note(:photon_conversion_veto, "photon-conversion veto: events in which a photon converts with reconstructed conversion radius sqrt(Rx^2 + Ry^2) < 2 cm are rejected; the conversion radius is not accessible from the DSL photon-selection block and is applied in the BOSS algorithm code")
         .note(:decay_angle_cut, "|cos(theta_decay)| < 0.9 required for the eta -> gamma gamma decay to suppress combinatorial gamma gamma background; the decay angle of the gamma pair in the eta rest frame is not expressible in the DSL and is applied in the BOSS algorithm code")
         .note(:babayaga_qed_background, "continuum QED background e+e- -> e+e-gamma(gamma) and 3gamma modelled with Babayaga on the 3.773 GeV continuum data sample; Babayaga is not an EvtGen decay card, so the sample is generated separately and analysed with the same selection")

sel_mode2 = Selection.new
sel_mode2.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp "==2"
           nChrn "==2"
           nNet  "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=2"                      # eta -> gamma gamma : at least two photons
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .remove([:lp <= :chrgp, :lm <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct eta from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 25
           neta ">=1"
         }
         .kinematic_fit([:lp, :lm, :pip, :pim, :eta]) {  # 4C fit to e+e- pi+ pi- eta
           nominal
           vertex_fit([2, 3])              # common vertex for pi+ (2) and pi- (3)
           constrain_four_momentum
           chi2_cut 200
         }

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.execute_on([jpsi_data, jpsi_incMC, cont_data, cont_incMC, exMC_mode2, exMC_bkg])

## ================= Mode III: eta -> pi+ pi- pi0, pi0 -> gamma gamma =================
alg3_name = "JpsiEtaEEPiPiPi0"
alg_mode3 = Algorithm.new(alg3_name)
alg_mode3.set_header(["#{alg3_name}Alg/#{alg3_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .note(:photon_conversion_veto, "photon-conversion veto: events in which a photon converts with reconstructed conversion radius sqrt(Rx^2 + Ry^2) < 2 cm are rejected; the conversion radius is not accessible from the DSL photon-selection block and is applied in the BOSS algorithm code")
         .note(:decay_angle_cut, "|cos(theta_decay)| < 0.9 required for the pi0 -> gamma gamma decay to suppress combinatorial gamma gamma background; the decay angle of the gamma pair in the pi0 rest frame is not expressible in the DSL and is applied in the BOSS algorithm code")
         .note(:babayaga_qed_background, "continuum QED background e+e- -> e+e-gamma(gamma) and 3gamma modelled with Babayaga on the 3.773 GeV continuum data sample; Babayaga is not an EvtGen decay card, so the sample is generated separately and analysed with the same selection")

sel_mode3 = Selection.new
sel_mode3.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp "==2"
           nChrn "==2"
           nNet  "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=2"                      # pi0 -> gamma gamma : at least two photons
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .remove([:lp <= :chrgp, :lm <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct pi0 from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=1"
         }
         .kinematic_fit([:lp, :lm, :pip, :pim, :pi0]) {  # 4C fit to e+e- pi+ pi- pi0
           nominal
           vertex_fit([2, 3])              # common vertex for pi+ (2) and pi- (3)
           constrain_four_momentum
           chi2_cut 200
         }

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)
alg_mode3.execute_on([jpsi_data, jpsi_incMC, cont_data, cont_incMC, exMC_mode3, exMC_bkg])

## ================= Mode IV: eta -> gamma gamma =================
alg4_name = "JpsiEtaEEGammaGamma"
alg_mode4 = Algorithm.new(alg4_name)
alg_mode4.set_header(["#{alg4_name}Alg/#{alg4_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .note(:photon_conversion_veto, "photon-conversion veto: events in which a photon converts with reconstructed conversion radius sqrt(Rx^2 + Ry^2) < 2 cm are rejected; the conversion radius is not accessible from the DSL photon-selection block and is applied in the BOSS algorithm code")
         .note(:decay_angle_cut, "|cos(theta_decay)| < 0.9 required for the eta -> gamma gamma decay to suppress combinatorial gamma gamma background; the decay angle of the gamma pair in the eta rest frame is not expressible in the DSL and is applied in the BOSS algorithm code")
         .note(:babayaga_qed_background, "continuum QED background e+e- -> e+e-gamma(gamma) and 3gamma modelled with Babayaga on the 3.773 GeV continuum data sample; Babayaga is not an EvtGen decay card, so the sample is generated separately and analysed with the same selection")

sel_mode4 = Selection.new
sel_mode4.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp "==1"                     # e+ e- eta final state: one positive track
           nChrn "==1"                     # one negative track
           nNet  "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=2"                      # eta -> gamma gamma
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct eta from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 25
           neta ">=1"
         }
         .kinematic_fit([:lp, :lm, :eta]) {              # 4C fit to e+e- eta
           nominal
           constrain_four_momentum
           chi2_cut 200
         }

alg_mode4.with_decay_card(decay_card_mode4).apply(sel_mode4)
alg_mode4.execute_on([jpsi_data, jpsi_incMC, cont_data, cont_incMC, exMC_mode4, exMC_bkg])

## ================= Mode V: pi0 -> gamma gamma =================
alg5_name = "JpsiPi0EEGammaGamma"
alg_mode5 = Algorithm.new(alg5_name)
alg_mode5.set_header(["#{alg5_name}Alg/#{alg5_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .note(:photon_conversion_veto, "photon-conversion veto: events in which a photon converts with reconstructed conversion radius sqrt(Rx^2 + Ry^2) < 2 cm are rejected; the conversion radius is not accessible from the DSL photon-selection block and is applied in the BOSS algorithm code")
         .note(:decay_angle_cut, "|cos(theta_decay)| < 0.9 required for the pi0 -> gamma gamma decay to suppress combinatorial gamma gamma background; the decay angle of the gamma pair in the pi0 rest frame is not expressible in the DSL and is applied in the BOSS algorithm code")
         .note(:babayaga_qed_background, "continuum QED background e+e- -> e+e-gamma(gamma) and 3gamma modelled with Babayaga on the 3.773 GeV continuum data sample; Babayaga is not an EvtGen decay card, so the sample is generated separately and analysed with the same selection")

sel_mode5 = Selection.new
sel_mode5.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp "==1"                     # e+ e- pi0 final state: one positive track
           nChrn "==1"                     # one negative track
           nNet  "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=2"                      # pi0 -> gamma gamma
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct pi0 from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=1"
         }
         .kinematic_fit([:lp, :lm, :pi0]) {             # 4C fit to e+e- pi0
           nominal
           constrain_four_momentum
           invariant_mass_of(:lp, :lm).within(0.0, 0.4) # M(e+e-) <= 0.4 GeV/c^2 for J/psi -> pi0 e+e-
           chi2_cut 200
         }

alg_mode5.with_decay_card(decay_card_mode5).apply(sel_mode5)
alg_mode5.execute_on([jpsi_data, jpsi_incMC, cont_data, cont_incMC, exMC_mode5, exMC_bkg])