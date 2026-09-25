# ============================================================
# Datasets: psi(3686) real data + inclusive MC
# ============================================================
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay cards (EvtGen syntax)
# ============================================================
# Mode I: psi(2S) -> omega K_S0 K+ pi-,  omega -> pi+pi-pi0, K_S0 -> pi+pi-, pi0 -> gamma gamma
decay_card_modeI = <<~DECAYCARD
  Decay psi(2S)
  1.0000 omega K_S0 K+ pi- PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: psi(2S) -> omega K+ K- pi0,  omega -> pi+pi-pi0, pi0 -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay psi(2S)
  1.0000 omega K+ K- pi0 PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ============================================================
# Exclusive MC samples (100k events each)
# ============================================================
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_omegaKsKpi_modeI_exMC"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_omegaKKpi0_modeII_exMC"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Mode I: psi(2S) -> omega K_S0 K+ pi-
#   final state: K_S0 K+ pi+ pi- pi- pi0
# ============================================================
alg_name_I = "OmegaKsKpi"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_I = Selection.new
sel_I.select_track {                     # charged-track pre-selection
        cos_theta 0.93                   # |cos(theta)| < 0.93
        Vz        20.0                   # |Vz| < 20 cm
        Vr        2.0                    # Vr < 2 cm
        nChrp     ">=2"                  # at least 2 positive tracks
        nChrn     ">=2"                  # at least 2 negative tracks
        nNet      "==0"                  # net charge zero
      }
      .pid(method: :probability) {       # kaon identification
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        nkp      ">=1"                   # at least one K+
      }
      .remove([:kp <= :chrgp])           # remove identified kaons from positive charged list
      .assign({:chrgp => :pip, :chrgn => :pim})  # treat remaining tracks as pions
      .select_photon {                   # photon selection
        tdc_emc_start     0              # TDC window 0-700 ns
        tdc_emc_end       14
        energyThreshold_b 0.025          # E > 25 MeV  (barrel)
        energyThreshold_e 0.050          # E > 50 MeV  (endcap)
        angle_to_track    10.0           # > 10 deg from any charged track
        nGam              ">=2"          # at least two photons (for pi0 -> gamma gamma)
      }
      .secondary_vertex_fit([:pip, :pim]) {   # reconstruct K_S0 from an opposite-charge pair
        build_virtual_particle(:K_S0).by_minimizing_mass_difference  # closest to nominal K_S0 mass
        remove_used_particle_from_candidate_list
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct pi0 (1-C mass constraint)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=1"                    # at least one pi0 candidate
      }
      .kinematic_fit([:K_S0, :kp, :pip, :pim, :pim, :pi0]) {   # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 40
      }

alg_I.with_decay_card(decay_card_modeI).apply(sel_I)
root_files_I = alg_I.execute_on([psip_data, psip_incMC, exMC_modeI])

# ============================================================
# Mode II: psi(2S) -> omega K+ K- pi0
#   final state: K+ K- pi+ pi- pi0 pi0
# ============================================================
alg_name_II = "OmegaKKpi0"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_II = Selection.new
sel_II.select_track {                     # charged-track pre-selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        20.0                    # |Vz| < 20 cm
        Vr        2.0                     # Vr < 2 cm
        nChrp     "==2"                   # exactly 2 positive tracks
        nChrn     "==2"                   # exactly 2 negative tracks
        nNet      "==0"                   # net charge zero
      }
      .pid(method: :probability) {        # kaon identification
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        nkp      ">=1"                    # at least one K+
        nkm      ">=1"                    # at least one K-
      }
      .remove([:kp <= :chrgp, :km <= :chrgn])  # remove identified kaons from charged lists
      .assign({:chrgp => :pip, :chrgn => :pim}) # treat remaining tracks as pions
      .select_photon {                    # photon selection
        tdc_emc_start     0               # TDC window 0-700 ns
        tdc_emc_end       14
        energyThreshold_b 0.025           # E > 25 MeV  (barrel)
        energyThreshold_e 0.050           # E > 50 MeV  (endcap)
        angle_to_track    10.0            # > 10 deg from any charged track
        nGam              ">=4"           # at least four photons (two pi0's)
      }
      # Competing hypotheses: store chi2 of 4C fits with 3gamma and 5gamma (no chi2 cut, non-nominal)
      .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
        constrain_four_momentum           # stored for ROOT-level background veto
      }
      .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
        constrain_four_momentum           # stored for ROOT-level background veto
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct pi0's (1-C mass constraint)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=2"                    # at least two pi0 candidates
      }
      .kinematic_fit([:kp, :km, :pip, :pim, :pi0, :pi0]) {  # 6C fit (4C + two pi0 mass constraints)
        nominal
        constrain_four_momentum
        chi2_cut 100
      }

alg_II.with_decay_card(decay_card_modeII).apply(sel_II)
root_files_II = alg_II.execute_on([psip_data, psip_incMC, exMC_modeII])