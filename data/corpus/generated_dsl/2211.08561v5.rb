# =============================================================================
#  e+e- -> KS0 KS0 J/psi   (J/psi -> e+e- or mu+mu-,  KS0 -> pi+pi-)
#  36 energy points in 4.130 - 4.946 GeV (real data + inclusive MC)
# =============================================================================

### ---------------------------- Datasets ---------------------------------- ###
# 36 scan points between 4.130 and 4.946 GeV: real data and matching inclusive MC
data_points  = DatasetManager.real_data.where(cms_energy: { value: 4130.0..4946.0 })
incMC_points = DatasetManager.inclusive_mc.where(cms_energy: { value: 4130.0..4946.0 })

# Reference data set used for the exclusive signal MC (produced on the psi(4260))
data_4260 = DatasetManager.real_data.find("703_4260")

### ------------------------- Decay cards (EvtGen) ------------------------- ###
# Signal: e+e- -> KS0 KS0 J/psi, J/psi -> e+e-
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K_S0 K_S0 J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Signal: e+e- -> KS0 KS0 J/psi, J/psi -> mu+mu-
decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K_S0 K_S0 J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### ------------ Exclusive signal MC (100k events, J/psi -> e+e-) ----------- ###
exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ksksjpsi_ee"
  config.related_dataset = data_4260
  config.events          = 100000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

### ================== Chain I : two reconstructed K_S0 ===================== ###
alg_name_2ks = "KSKSJpsiTwoKS"
alg_2ks = Algorithm.new(alg_name_2ks)
alg_2ks.set_header(["#{alg_name_2ks}Alg/#{alg_name_2ks}.h"])
       .set_constant({"ECMS" => [:double, 4.26]})

sel_2ks = Selection.new
sel_2ks.select_track {                       # charged-track selection
          cos_theta 0.93                     # |cos theta| < 0.93
          Vz        10.0                     # |Vz| < 10 cm
          Vr        1.0                      # Vr  < 1 cm
          nChrp     ">=3"                    # >= 6 tracks (2 KS0 -> 4 pi, J/psi -> 2 l)
          nChrn     ">=3"
        }
        .select_photon {                     # photon acceptance
          tdc_emc_start     0                # EMC timing 0 - 700 ns
          tdc_emc_end       14
          angle_to_track    10.0             # >= 10 deg from any charged track
          energyThreshold_b 0.025            # 25 MeV (barrel)
          energyThreshold_e 0.050            # 50 MeV (endcap)
        }
        .pid(method: :probability) {         # probability PID
          prob_cut 0.001
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,   # p > 0.95 GeV -> lepton
                                         treat_as_electron_if_energy_above: 0.95    # EMC E > 0.95 GeV -> e, else mu
        }
        .remove([:lp <= :chrgp, :lm <= :chrgn])    # keep leptons out of the pion pool
        .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks treated as pions
        # first K_S0 from an opposite-charge pion pair (best mass difference)
        .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        # second K_S0 from the remaining pion pair
        .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        # 4C kinematic fit of KS0 KS0 l+l-
        .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) {
          nominal
          constrain_four_momentum
          chi2_cut 200                       # chi2 < 200
        }

alg_2ks.with_decay_card(decay_card_ee).apply(sel_2ks)

### ================== Chain II : one reconstructed K_S0 ==================== ###
alg_name_1ks = "KSKSJpsiOneKS"
alg_1ks = Algorithm.new(alg_name_1ks)
alg_1ks.set_header(["#{alg_name_1ks}Alg/#{alg_name_1ks}.h"])
       .set_constant({"ECMS" => [:double, 4.26]})

sel_1ks = Selection.new
sel_1ks.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     ">=2"                    # >= 4 tracks (1 KS0 -> 2 pi, J/psi -> 2 l)
          nChrn     ">=2"
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=2"            # at least two photons
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                         treat_as_electron_if_energy_above: 0.95
        }
        .remove([:lp <= :chrgp, :lm <= :chrgn])
        .assign({:chrgp => :pip, :chrgn => :pim})
        # one reconstructed K_S0 from an opposite-charge pion pair
        .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        # 1C fit: the second (undetected) K_S0 enters with its nominal mass, momentum floated
        .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) {
          nominal
          miss_track_of :K_S0                # missing K_S0 constrained to the K_S0 mass
          constrain_four_momentum
          chi2_cut 20                        # chi2 < 20
        }

alg_1ks.with_decay_card(decay_card_ee).apply(sel_1ks)

### ------------------------ Execute on all datasets ----------------------- ###
datasets = data_points + incMC_points + [exMC_signal_ee]
alg_2ks.execute_on(datasets)
alg_1ks.execute_on(datasets)