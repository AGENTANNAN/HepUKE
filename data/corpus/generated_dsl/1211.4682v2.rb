# =====================================================================
# BOSS spec: J/psi -> Lambda Lambdabar pi0/eta  and  psi(2S) -> Lambda Lambdabar pi0/eta
# Lambda -> p pi-, Lambdabar -> anti-p pi+, pi0/eta -> gamma gamma
# Covers dataset prep (4 decay cards + 4 exclusive MC) and the shared
# event-selection chain up to and including the 4C kinematic fit.
# Mode-specific vetoes (mass windows / recoil mass) are applied downstream
# in the ROOT analysis and are therefore NOT encoded here.
# =====================================================================

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")        # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # psi(2S) inclusive MC

### Decay cards (EvtGen format) ###
# J/psi -> Lambda Lambdabar pi0
decay_card_jpsi_pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 pi0   PHSP;
  Enddecay
  Decay Lambda0
  1.000 p+ pi-                      HypWK;
  Enddecay
  Decay anti-Lambda0
  1.000 anti-p- pi+                 HypWK;
  Enddecay
  Decay pi0
  1.000 gamma gamma                 PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> Lambda Lambdabar eta
decay_card_jpsi_eta = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 eta   PHSP;
  Enddecay
  Decay Lambda0
  1.000 p+ pi-                      HypWK;
  Enddecay
  Decay anti-Lambda0
  1.000 anti-p- pi+                 HypWK;
  Enddecay
  Decay eta
  1.000 gamma gamma                 PHSP;
  Enddecay
  End
DECAYCARD

# psi(2S) -> Lambda Lambdabar pi0
decay_card_psip_pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda0 anti-Lambda0 pi0   PHSP;
  Enddecay
  Decay Lambda0
  1.000 p+ pi-                      HypWK;
  Enddecay
  Decay anti-Lambda0
  1.000 anti-p- pi+                 HypWK;
  Enddecay
  Decay pi0
  1.000 gamma gamma                 PHSP;
  Enddecay
  End
DECAYCARD

# psi(2S) -> Lambda Lambdabar eta
decay_card_psip_eta = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda0 anti-Lambda0 eta   PHSP;
  Enddecay
  Decay Lambda0
  1.000 p+ pi-                      HypWK;
  Enddecay
  Decay anti-Lambda0
  1.000 anti-p- pi+                 HypWK;
  Enddecay
  Decay eta
  1.000 gamma gamma                 PHSP;
  Enddecay
  End
DECAYCARD

### Exclusive MC (200k events per mode) ###
exMC_jpsi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_lambdalambdabar_pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_jpsi_pi0
  config.cross_section   = :default
end

exMC_jpsi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_lambdalambdabar_eta"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_jpsi_eta
  config.cross_section   = :default
end

exMC_psip_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_lambdalambdabar_pi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_psip_pi0
  config.cross_section   = :default
end

exMC_psip_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_lambdalambdabar_eta"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_psip_eta
  config.cross_section   = :default
end

# =====================================================================
# Common event-selection chain (identical for all four modes)
# =====================================================================
event_selection_common = Selection.new
  .select_track {                     # charged-track quality + multiplicity
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        100.0                   # |Vz| < 100 cm
    Vr        10.0                    # |Vr| < 10 cm
    nChrp     ">=2"                   # at least two positive tracks
    nChrn     ">=2"                   # at least two negative tracks
  }
  .select_photon {                    # photon selection
    tdc_emc_start     0               # EMC TDC window start
    tdc_emc_end       14              # EMC TDC window end
    energyThreshold_b 0.025           # 25 MeV barrel threshold
    energyThreshold_e 0.050           # 50 MeV endcap threshold
    nGam              ">=2"           # at least two photons (pi0/eta -> gamma gamma)
  }
  .pid(method: :probability) {        # PID by the probability method
    prob_cut 0.001                    # probability > 0.001
    identify :proton, against: [:kaon, :pion]  # p+ and anti-p- vs K and pi
    nprp ">=1"                        # at least one proton
    nprm ">=1"                        # at least one antiproton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])    # pull (anti)protons out of the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})    # remaining positive/negative -> pi+/pi-
  .for_each(:prp) { where { pt < 0.2 }; remove }   # proton pT > 0.2 GeV/c
  .for_each(:prm) { where { pt < 0.2 }; remove }   # antiproton pT > 0.2 GeV/c
  .secondary_vertex_fit([:prp, :pim]) {          # Lambda: p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # Lambdabar: anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {   # 4C fit to Lambda Lambdabar gamma gamma
    nominal                          # nominal fit (only fitted momenta retained)
    constrain_four_momentum          # 4C energy-momentum conservation
    chi2_cut 200                     # loose chi2 cut; tighter cut applied downstream
  }

# =====================================================================
# One algorithm per mode; all four share the common selection chain
# =====================================================================

# ---- J/psi -> Lambda Lambdabar pi0 ----
alg_jpsi_pi0 = Algorithm.new("JpsiLLbarPi0")
alg_jpsi_pi0.set_header(["JpsiLLbarPi0Alg/JpsiLLbarPi0.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
alg_jpsi_pi0.with_decay_card(decay_card_jpsi_pi0).apply(event_selection_common.dup)
alg_jpsi_pi0.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_pi0])

# ---- J/psi -> Lambda Lambdabar eta ----
alg_jpsi_eta = Algorithm.new("JpsiLLbarEta")
alg_jpsi_eta.set_header(["JpsiLLbarEtaAlg/JpsiLLbarEta.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
alg_jpsi_eta.with_decay_card(decay_card_jpsi_eta).apply(event_selection_common.dup)
alg_jpsi_eta.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_eta])

# ---- psi(2S) -> Lambda Lambdabar pi0 ----
alg_psip_pi0 = Algorithm.new("PsipLLbarPi0")
alg_psip_pi0.set_header(["PsipLLbarPi0Alg/PsipLLbarPi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
alg_psip_pi0.with_decay_card(decay_card_psip_pi0).apply(event_selection_common.dup)
alg_psip_pi0.execute_on([psip_data, psip_incMC, exMC_psip_pi0])

# ---- psi(2S) -> Lambda Lambdabar eta ----
alg_psip_eta = Algorithm.new("PsipLLbarEta")
alg_psip_eta.set_header(["PsipLLbarEtaAlg/PsipLLbarEta.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
alg_psip_eta.with_decay_card(decay_card_psip_eta).apply(event_selection_common.dup)
alg_psip_eta.execute_on([psip_data, psip_incMC, exMC_psip_eta])