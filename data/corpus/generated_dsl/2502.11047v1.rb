# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Single-tag search for Lambda_c+ -> Sigma0 K+ pi0 and Lambda_c+ -> Sigma0 K+ pi+ pi-
# Seven energy points from 4599.53 to 4698.82 MeV (~4.5 fb^-1)
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

data_points  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card for Lambda_c+ -> Sigma0 K+ pi0 (Sigma0 -> Lambda gamma, Lambda -> p pi-, pi0 -> gamma gamma)
decay_card_mode_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma0 K+ pi0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 Lambda0 gamma PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Lambda_c+ -> Sigma0 K+ pi+ pi- (Sigma0 -> Lambda gamma, Lambda -> p pi-)
decay_card_mode_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma0 K+ pi+ pi- PHSP;
    Enddecay

    Decay Sigma0
    1.0000 Lambda0 gamma PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 1,000,000 events for each mode at each energy point
exMCs_pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_Sigma0Kpi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card_mode_pi0
  config.cross_section = :default
end

exMCs_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_Sigma0Kpipipi"
  config.events        = 1_000_000
  config.decay_card    = decay_card_mode_pipi
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ------------------------------------------------------------------
# Mode I: Lambda_c+ -> Sigma0 K+ pi0 (Sigma0 -> Lambda gamma, Lambda -> p pi-, pi0 -> gamma gamma)
# ------------------------------------------------------------------
alg_name_pi0 = "LcToSigma0Kpi0"
alg_pi0 = Algorithm.new(alg_name_pi0)
alg_pi0.set_header(["#{alg_name_pi0}Alg/#{alg_name_pi0}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})

sel_pi0 = Selection.new
sel_pi0.select_track {                       # charged track selection
          cos_theta 0.93                     # |cos(theta)| < 0.93
          Vz 10.0                            # |Vz| < 10 cm
          Vr 1.0                             # Vr < 1 cm
          nChrp ">=1"                        # >= 1 positive track
          nChrn ">=1"                        # >= 1 negative track (>= 2 total)
        }
        .pid(method: :probability) {         # identify p, K+, pi- (final state)
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]   # p+/p-
          identify :kaon,   against: [:pion, :proton] # K+/K-
          identify :pion,   against: [:proton, :kaon] # pi+/pi-
          nprp ">=1"
          nkp  ">=1"
          npim ">=1"
        }
        .select_photon {                     # photon selection
          tdc_emc_start 0                    # EMC timing [0,700] ns
          tdc_emc_end 14
          angle_to_track 10.0                # > 10 deg from nearest charged track
          energyThreshold_b 0.025            # 25 MeV (barrel)
          energyThreshold_e 0.050            # 50 MeV (endcap)
          nGam ">=3"                         # >= 3 photons (pi0 mode)
        }
        .secondary_vertex_fit([:prp, :pim]) {   # Lambda -> p pi- secondary vertex
          build_virtual_particle(:Lambda).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C mass-constrained fit -> pi0 (nominal)
          invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)          # M(gamma gamma) window
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 1C: M(gg) = m(pi0)
          chi2_cut 200
          npi0 ">=1"
        }
        .kinematic_fit([:Lambda, :kp, :pi0, :gamma]) {   # nominal fit, corrected four-momenta saved
          nominal
          constrain_four_momentum
          chi2_cut 200
        }

# BOSS-side procedures / criteria that cannot be expressed in the DSL surface
alg_pi0
  .note(:background_veto, "Lambda -> p pi- secondary-vertex fit quality cuts applied: vertex chi2 < 100 and decay length > 2; |M(p pi-) - 1.1157| < 0.005 GeV/c^2 mass window on the p pi- pair")
  .note(:efficiency_curve, "Sigma0 candidate selection: M(Lambda gamma) in [1.179,1.203] GeV/c^2 and cos(angle(p,gamma)) < cos(20 deg) applied before the nominal kinematic fit")
  .note(:background_veto, "Delta E(p pi- K+ gamma gamma) in [-0.160,-0.030] GeV; best Lambda_c+ candidate chosen by minimum |Delta E|; final DeltaE window [-27,6] MeV and M_BC in [2.282,2.291] GeV/c^2 applied in the downstream ROOT analysis")

alg_pi0.with_decay_card(decay_card_mode_pi0).apply(sel_pi0)
alg_pi0.execute_on(data_points + incMC_points + exMCs_pi0)

# ------------------------------------------------------------------
# Mode II: Lambda_c+ -> Sigma0 K+ pi+ pi- (Sigma0 -> Lambda gamma, Lambda -> p pi-)
# ------------------------------------------------------------------
alg_name_pipi = "LcToSigma0Kpipipi"
alg_pipi = Algorithm.new(alg_name_pipi)
alg_pipi.set_header(["#{alg_name_pipi}Alg/#{alg_name_pipi}.h"])
        .set_constant({"ECMS" => [:double, 4.600]})

sel_pipi = Selection.new
sel_pipi.select_track {                      # charged track selection
          cos_theta 0.93                     # |cos(theta)| < 0.93
          Vz 10.0                            # |Vz| < 10 cm
          Vr 1.0                             # Vr < 1 cm
          nChrp ">=2"                        # >= 2 positive tracks
          nChrn ">=2"                        # >= 2 negative tracks (>= 4 total)
        }
        .pid(method: :probability) {         # identify p, K+, pi+, pi-
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]   # p+/p-
          identify :kaon,   against: [:pion, :proton] # K+/K-
          identify :pion,   against: [:proton, :kaon] # pi+/pi-
          nprp ">=1"
          nkp  ">=1"
          npip ">=1"
          npim ">=1"
        }
        .select_photon {                     # photon selection
          tdc_emc_start 0                    # EMC timing [0,700] ns
          tdc_emc_end 14
          angle_to_track 10.0                # > 10 deg from nearest charged track
          energyThreshold_b 0.025            # 25 MeV (barrel)
          energyThreshold_e 0.050            # 50 MeV (endcap)
          nGam ">=1"                         # >= 1 photon (Sigma0 photon)
        }
        .secondary_vertex_fit([:prp, :pim]) {   # Lambda -> p pi- secondary vertex
          build_virtual_particle(:Lambda).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
        .kinematic_fit([:Lambda, :kp, :pip, :pim, :gamma]) {   # nominal fit
          nominal
          constrain_four_momentum
          chi2_cut 200
        }

alg_pipi
  .note(:background_veto, "Lambda -> p pi- secondary-vertex fit quality cuts applied: vertex chi2 < 100 and decay length > 2; |M(p pi-) - 1.1157| < 0.005 GeV/c^2 mass window on the p pi- pair")
  .note(:efficiency_curve, "Sigma0 candidate selection: M(Lambda gamma) in [1.179,1.203] GeV/c^2 applied before the nominal kinematic fit")
  .note(:background_veto, "Delta E(p pi- K+ pi+ pi-) < -0.040 GeV; best Lambda_c+ candidate chosen by minimum |Delta E|; final DeltaE window [-21,7] MeV and M_BC in [2.282,2.291] GeV/c^2 applied in the downstream ROOT analysis")

alg_pipi.with_decay_card(decay_card_mode_pipi).apply(sel_pipi)
alg_pipi.execute_on(data_points + incMC_points + exMCs_pipi)