# psi -> pi+ pi- eta', eta' -> eta pi+ pi-, eta -> gamma gamma
# Partial wave analyses at J/psi (1.31e9) and psi(3686) (4.48e8).

### Dataset preparation ###
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")

# Signal decay card for J/psi -> pi+ pi- eta' (eta' -> eta pi+ pi-, eta -> gg)
decay_card_jpsi = <<~DECAYCARD
    Decay J/psi
    1.0000 pi+ pi- eta'    PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-     PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma     PHSP;
    Enddecay

    End
DECAYCARD

decay_card_psip = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- eta'    PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-     PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma     PHSP;
    Enddecay

    End
DECAYCARD

exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_pipietap"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Psip_pipietap"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

### Event selection ###

# ---- J/psi -> pi+ pi- eta' ----
alg_jpsi = Algorithm.new("JpsiPiPiEtap")
alg_jpsi.set_header(["JpsiPiPiEtapAlg/JpsiPiPiEtap.h"])
        .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi = Selection.new
sel_jpsi.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     "==2"
          nChrn     "==2"
          nNet      "==0"
        }
        .select_photon {
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          angle_to_track    10.0
          tdc_emc_start     0
          tdc_emc_end       14
          nGam              ">=2"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :pion, against: [:kaon]
          npip ">=2"
          npim ">=2"
        }
        .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
          nominal
          constrain_four_momentum
          chi2_cut 40
        }

alg_jpsi.with_decay_card(decay_card_jpsi).apply(sel_jpsi)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])


# ---- psi(3686) -> pi+ pi- eta' ----
alg_psip = Algorithm.new("PsipPiPiEtap")
alg_psip.set_header(["PsipPiPiEtapAlg/PsipPiPiEtap.h"])
        .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip = Selection.new
sel_psip.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     "==2"
          nChrn     "==2"
          nNet      "==0"
        }
        .select_photon {
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          angle_to_track    10.0
          tdc_emc_start     0
          tdc_emc_end       14
          nGam              ">=2"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :pion, against: [:kaon]
          npip ">=2"
          npim ">=2"
        }
        .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
          nominal
          constrain_four_momentum
          chi2_cut 40
        }

alg_psip.with_decay_card(decay_card_psip).apply(sel_psip)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])
