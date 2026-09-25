# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3686) signal point (real data + inclusive MC)
data_3686  = DatasetManager.real_data.find("709_3686")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
# psi(3770) data + inclusive MC at 3.773 GeV, used for QED continuum estimation
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for Mode I: psi(3686) -> Lambda anti-Lambda eta', eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: psi(3686) -> Lambda anti-Lambda eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# 1M-event exclusive MC for each eta' decay mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_LLbar_etap_gammapipim"
  config.related_dataset = data_3686
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_LLbar_etap_etapipim"
  config.related_dataset = data_3686
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Mode I: eta' -> gamma pi+ pi- ###
alg_name_I = "LambdaLambdabarEtapGammaPiPi"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {                 # Charged track selection
            cos_theta 0.93               # |cos(theta)| < 0.93
            Vz 10.0                      # |Vz| < 10 cm
            Vr 1.0                       # Vr < 1 cm (transverse)
            nChrp ">=3"                  # at least 3 positive tracks
            nChrn ">=3"                  # at least 3 negative tracks
            nTot  ">5"                   # more than 5 tracks total
         }
         .select_photon {                # Photon selection
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0          # at least 10 deg from any charged track
            energyThreshold_b 0.025      # 25 MeV barrel
            energyThreshold_e 0.050      # 50 MeV endcap
            nGam ">=1"                   # at least one photon (Mode I)
         }
         .pid(method: :probability) {    # Probability PID, 0.001 cut
            prob_cut 0.001
            identify :proton, against: [:pion, :kaon]   # p+ and anti-p- vs pi/K
            nprp ">=1"                                  # at least one proton
            nprm ">=1"                                  # at least one anti-proton
            identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K/proton
         }
         .remove([:prp <= :chrgp, :prm <= :chrgn])      # remove identified (anti)protons
         .secondary_vertex_fit([:prp, :pim]) {          # Lambda -> p pi-
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         # 4C nominal kinematic fit over Lambda anti-Lambda gamma pi+ pi-; picks best Lambda anti-Lambda pair
         .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :pip, :pim]) {
            nominal                        # nominal fit (corrected four-momenta used)
            vertex_fit([3, 4])             # primary-vertex fit of the prompt eta' pion pair
            constrain_four_momentum        # 4C energy-momentum constraint
            chi2_cut 30
         }

alg_modeI
  .note(:background_veto, "veto events consistent with Sigma0 (M(Lambda gamma)), J/psi and chi_c (M(Lambda anti-Lambda)), Xi and Sigma(1385) (M(Lambda pi)) peaking backgrounds; windows derived from signal vs inclusive/continuum MC and applied as mass-window vetoes")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([data_3686, incMC_3686, data_3773, incMC_3773, exMC_modeI])

### Mode II: eta' -> eta pi+ pi-, eta -> gamma gamma ###
alg_name_II = "LambdaLambdabarEtapEtaPiPi"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {                # Charged track selection (same as Mode I)
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp ">=3"
            nChrn ">=3"
            nTot  ">5"
         }
         .select_photon {                # Photon selection
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"                   # at least two photons (Mode II, from eta -> gamma gamma)
         }
         .pid(method: :probability) {    # Probability PID, 0.001 cut
            prob_cut 0.001
            identify :proton, against: [:pion, :kaon]   # p+ and anti-p- vs pi/K
            nprp ">=1"
            nprm ">=1"
            identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K/proton
         }
         .remove([:prp <= :chrgp, :prm <= :chrgn])      # remove identified (anti)protons
         .secondary_vertex_fit([:prp, :pim]) {          # Lambda -> p pi-
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct eta from gamma gamma (1-C)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 40
            neta ">=1"                                  # at least one eta candidate
         }
         # 4C nominal fit over Lambda anti-Lambda eta pi+ pi-, with (pi+ pi- eta) constrained to eta'
         .kinematic_fit([:Lambda, :Lambda_bar, :eta, :pip, :pim]) {
            nominal                                     # nominal fit
            vertex_fit([3, 4])                          # primary-vertex fit of the prompt eta' pion pair
            constrain_four_momentum                     # 4C energy-momentum constraint
            invariant_mass_of(:pip, :pim, :eta).constrain_to_nominal_mass_of(:etap)
            chi2_cut 40
         }

alg_modeII
  .note(:background_veto, "veto events consistent with Sigma0 (M(Lambda gamma)), J/psi and chi_c (M(Lambda anti-Lambda)), Xi and Sigma(1385) (M(Lambda pi)) peaking backgrounds; windows derived from signal vs inclusive/continuum MC and applied as mass-window vetoes")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([data_3686, incMC_3686, data_3773, incMC_3773, exMC_modeII])