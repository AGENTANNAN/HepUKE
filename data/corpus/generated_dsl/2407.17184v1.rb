### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data @ 3.686 GeV (2.712e9 psi(3686) events)
cont_data  = DatasetManager.real_data.find("709_3650")     # continuum data @ 3.650 GeV (401 pb^-1)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # inclusive MC @ 3.686 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")  # inclusive MC @ 3.650 GeV

# Decay card (Mode I): psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K+ K- eta', eta' -> pi+ pi- gamma
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c(2S)  PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000 K+ K- eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card (Mode II): psi(3686) -> gamma chi_c1, chi_c1 -> K+ K- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1  P2GC1;
    Enddecay

    Decay chi_c1
    1.0000 K+ K- eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC sample for each of the two modes
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2s_etap_to_pipigam"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic1_etap_to_pipieta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### ------------------------------------------------------------------ ###
### Mode I: eta' -> pi+ pi- gamma                                      ###
### ------------------------------------------------------------------ ###
alg_name_modeI = "PsipToGammaEtac2S"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                    # charged track quality + multiplicity
     cos_theta   0.93                # |cos(theta)| < 0.93
     Vz          10.0                # |Vz| < 10 cm
     Vr          1.0                 # Vr < 1 cm in the transverse plane
     nChrp       ">=2"               # at least two positive tracks
     nChrn       ">=2"               # at least two negative tracks
  }
  .select_photon {                   # photon selection
     tdc_emc_start     0            # EMC timing window 0 - 700 ns
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025         # 25 MeV in the barrel
     energyThreshold_e 0.050         # 50 MeV in the endcap
     nGam              ">=2"         # eta' -> pi+pi-gamma : at least two photons
  }
  .pid(method: :chi2_sum) {          # dE/dx + TOF combined: minimise chi2_PID over track species
     chi_min_cut 4
     identify :kaon, :pion           # K+ K- pi+ pi- final state
  }
  # Nominal 4C fit to K+ K- pi+ pi- gamma with the eta' mass pre-window
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
     nominal
     invariant_mass_of(:pip, :pim, :gamma).within(0.94, 0.97)   # M(pi+pi-gamma) = M(eta') window
     constrain_four_momentum
     chi2_cut 25
  }
  # Final 3C fit (constrain three-momentum only -> radiative photon energy left free),
  # performed on the candidate already selected by the nominal fit; gives M(K+K-eta')
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
     use_track_index_from_nominal_kmfit
     constrain_three_momentum
  }

alg_modeI.note(:background_veto,
               "Mode I: vetoes against pi0, eta, J/psi, chi_cJ, pi+pi-J/psi and phi backgrounds " \
               "on the K+K-pi+pi-gamma final state. The veto windows are not specified in the " \
               "analysis description and must be applied on the stored invariant masses.")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([psip_data, cont_data, psip_incMC, cont_incMC, exMC_modeI])

### ------------------------------------------------------------------ ###
### Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma                   ###
### ------------------------------------------------------------------ ###
alg_name_modeII = "PsipToGammaChicJ"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                    # charged track quality + multiplicity
     cos_theta   0.93                # |cos(theta)| < 0.93
     Vz          10.0                # |Vz| < 10 cm
     Vr          1.0                 # Vr < 1 cm in the transverse plane
     nChrp       ">=2"               # at least two positive tracks
     nChrn       ">=2"               # at least two negative tracks
  }
  .select_photon {                   # photon selection
     tdc_emc_start     0            # EMC timing window 0 - 700 ns
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025         # 25 MeV in the barrel
     energyThreshold_e 0.050         # 50 MeV in the endcap
     nGam              ">=3"         # eta' -> pi+pi-eta, eta -> gamma gamma : at least three photons
  }
  .pid(method: :chi2_sum) {          # dE/dx + TOF combined: minimise chi2_PID over track species
     chi_min_cut 4
     identify :kaon, :pion           # K+ K- pi+ pi- final state
  }
  # 1C Kalman fit over all photon pairs: constrain M(gamma gamma) to the nominal eta mass;
  # the best eta candidate (smallest chi2) is kept and required to exist
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).within(0.52, 0.57)                 # M(gamma gamma) eta window
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta) # 1C mass constraint
     chi2_cut 25
     neta ">=1"
  }
  # Nominal 4C fit to K+ K- pi+ pi- eta with the eta' mass pre-window
  .kinematic_fit([:kp, :km, :pip, :pim, :eta]) {
     nominal
     invariant_mass_of(:pip, :pim, :eta).within(0.93, 0.98)   # M(pi+pi-eta) = M(eta') window
     constrain_four_momentum
     chi2_cut 20
  }
  # Final 3C fit on the candidate selected by the nominal fit; gives M(K+K-eta')
  .kinematic_fit([:kp, :km, :pip, :pim, :eta]) {
     use_track_index_from_nominal_kmfit
     constrain_three_momentum
  }

alg_modeII.note(:background_veto,
                "Mode II: vetoes against pi0, eta-recoil J/psi and phi backgrounds on the " \
                "K+K-pi+pi-eta final state. The veto windows are not specified in the analysis " \
                "description and must be applied on the stored invariant masses.")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([psip_data, cont_data, psip_incMC, cont_incMC, exMC_modeII])