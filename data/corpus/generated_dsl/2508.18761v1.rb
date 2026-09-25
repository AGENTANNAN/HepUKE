# Core DSL classes and dependencies are loaded automatically at execution

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# ---- Decay cards (EvtGen) ----
# Mode I : psi(2S) -> gamma chi_c0, chi_c0 -> Lambda anti-Lambda eta', eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
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

# Mode II: psi(2S) -> gamma chi_c0, chi_c0 -> Lambda anti-Lambda eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
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

# ---- Exclusive MC samples (100k events each, both through chi_c0) ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic0_LLbareta_etap_gammapipi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic0_LLbareta_etap_etapipi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection -- Mode I : eta' -> gamma pi+ pi- ###
alg_name_modeI = "LLbarEtaPrimeGam"
alg_modeI = Algorithm.new(alg_name_modeI)                     # independent decay mode -> its own Algorithm
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})          # 3.686 GeV psi(2S)
         .note(:background_veto, "Sigma0 background suppressed with M(Lambda gamma) > 1.2 GeV/c^2; J/psi and pi0 mass windows vetoed; Lambda/anti-Lambda mass window +/-5 MeV/c^2 and eta' window +/-12 MeV/c^2 applied in the ROOT stage")

sel_modeI = Selection.new
  .select_track {                                            # charged track selection
      cos_theta 0.93                                         # |cos(theta)| < 0.93
      Vz        10.0                                         # |Vz| < 10 cm
      Vr        1.0                                          # Vr < 1 cm
      nChrp     ">=3"                                        # at least 3 positive tracks
      nChrn     ">=3"                                        # at least 3 negative tracks
  }
  .select_photon {                                           # photon selection
      energyThreshold_b 0.025                                # EMC barrel energy > 25 MeV
      energyThreshold_e 0.050                                # EMC endcap energy > 50 MeV
      angle_to_track    10.0                                 # angle to nearest charged track > 10 deg
      tdc_emc_start     0                                    # EMC timing window
      tdc_emc_end       14
      nGam              ">=2"                                # at least 2 photons (radiative + eta' -> gamma)
  }
  .pid(method: :probability) {                               # PID, probability method
      prob_cut 0.001                                         # probability > 0.001
      identify :proton, against: [:kaon, :pion]              # p and p-bar vs K/pi
      nprp     ">=1"
      nprm     ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])                  # take protons out of the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})                  # remaining tracks treated as pi+/pi-
  .secondary_vertex_fit([:prp, :pim]) {                      # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {                      # anti-Lambda -> p-bar pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :pip, :pim]) {   # 4C fit to Lambda anti-Lambda gamma pi+ pi-
      nominal                                                # nominal fit (corrected four-momenta saved)
      constrain_four_momentum                                # 4-momentum conservation
      chi2_cut 18                                            # chi^2 < 18
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

### Event selection -- Mode II : eta' -> eta pi+ pi-, eta -> gamma gamma ###
alg_name_modeII = "LLbarEtaPrimeEta"
alg_modeII = Algorithm.new(alg_name_modeII)                   # independent decay mode -> its own Algorithm
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .note(:background_veto, "Sigma0 background suppressed with M(Lambda gamma) > 1.2 GeV/c^2; J/psi and pi0 mass windows vetoed; Lambda/anti-Lambda mass window +/-5 MeV/c^2 and eta' window +/-10 MeV/c^2 applied in the ROOT stage")

sel_modeII = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=3"
      nChrn     ">=3"
  }
  .select_photon {
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
      tdc_emc_start     0
      tdc_emc_end       14
      nGam              ">=3"                                # at least 3 photons (radiative + eta -> gamma gamma)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp     ">=1"
      nprm     ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) {                      # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {                      # anti-Lambda -> p-bar pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                  # reconstruct eta from a gamma-gamma pair (1C mass)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 200
      neta     ">=1"
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :eta, :pip, :pim]) { # 5C fit (4C + eta mass) to Lambda anti-Lambda eta pi+ pi-
      nominal
      constrain_four_momentum
      chi2_cut 53                                            # chi^2 < 53
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])