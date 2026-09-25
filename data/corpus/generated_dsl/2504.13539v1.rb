### Dataset preparation ###
# Representative energy point of the 4.258–4.681 GeV R-scan: 4.681 GeV
# BESIII sample-name convention [BOSS version]_[CMS energy in MeV] -> 706_4680 (Ecms = 4681.92 MeV)
data_4681  = DatasetManager.real_data.find("706_4680")
incMC_4681 = DatasetManager.inclusive_mc.find("706_4680")

# Decay card (representative eta_c hadronic mode: p pbar):
# e+e- -> gamma eta eta_c,  eta -> gamma gamma,  eta_c -> p pbar
decay_card_gamma_eta_etac = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma eta eta_c PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card (representative eta_c hadronic mode: p pbar):
# e+e- -> gamma eta' eta_c,  eta' -> eta pi+ pi-,  eta -> gamma gamma,  eta_c -> p pbar
decay_card_gamma_etap_etac = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma eta' eta_c PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k events for each of the two signal modes
exMC_gamma_eta_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4681_gamma_eta_etac_ppbar"
  config.related_dataset = data_4681
  config.events          = 100000
  config.decay_card      = decay_card_gamma_eta_etac
  config.cross_section   = :default
end
exMC_gamma_eta_etac.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_gamma_etap_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4681_gamma_etap_etac_ppbar"
  config.related_dataset = data_4681
  config.events          = 100000
  config.decay_card      = decay_card_gamma_etap_etac
  config.cross_section   = :default
end
exMC_gamma_etap_etac.save_to_config(format: :yaml, file_path: 'temp_for_test')


### Event selection (BOSS) ###
# ================= Signal mode I: e+e- -> gamma eta eta_c (eta -> gamma gamma, eta_c -> p pbar) =================
alg_name_1 = "GammaEtaEtac"
alg_1 = Algorithm.new(alg_name_1)
alg_1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
     .set_constant({"ECMS" => [:double, 4.681]})   # CMS energy 4.681 GeV
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_1 = Selection.new
  .select_track {                # Charged track selection
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     ">=1"              # at least one positive track (p)
    nChrn     ">=1"              # at least one negative track (pbar)
  }
  .select_photon {               # Photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025      # E_gamma > 25 MeV (barrel)
    energyThreshold_e 0.050
    angle_to_track    10.0       # angle to nearest charged track > 10 deg
    nGam              ">=3"      # at least three photons (radiative + 2 from eta)
  }
  .pid(method: :probability) {   # PID: probability method, prob_cut 0.001
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # p / pbar
    identify :pion,   against: [:kaon, :proton]   # pi+ / pi-
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # eta -> gamma gamma (1C mass-constrained fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:gamma, :eta, :prp, :prm]) {    # nominal 4C fit to gamma eta p pbar
    nominal
    constrain_four_momentum
    chi2_cut 20
  }
  .kinematic_fit([:gamma, :gamma, :prp, :prm]) {  # competing hypothesis "2gamma + hadrons" (store 4C chi2 only)
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prp, :prm]) {  # competing hypothesis "4gamma + hadrons"
    constrain_four_momentum
  }

alg_1.note(:ks0_secondary_vertex,
           "For the eta_c modes containing a K_S0 (16-mode combination), reconstruct K_S0 by a " \
           "secondary_vertex_fit of the pi+pi- pair, requiring |M(pi+pi-) - m(K_S0)| < 20 MeV and " \
           "decay length > 2 sigma; not applicable to the representative eta_c -> p pbar mode.")
     .note(:pi0_mass_constraint,
           "For the eta_c modes containing a pi0, constrain the gamma-gamma invariant mass to the " \
           "nominal pi0 mass, requiring |M(gamma gamma) - m(pi0)| < 15 MeV; not applicable to the " \
           "representative eta_c -> p pbar mode.")

alg_1.with_decay_card(decay_card_gamma_eta_etac).apply(sel_1)
root_files_1 = alg_1.execute_on([data_4681, incMC_4681, exMC_gamma_eta_etac])


# ================= Signal mode II: e+e- -> gamma eta' eta_c (eta' -> eta pi+ pi-, eta -> gamma gamma, eta_c -> p pbar) =================
alg_name_2 = "GammaEtapEtac"
alg_2 = Algorithm.new(alg_name_2)
alg_2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
     .set_constant({"ECMS" => [:double, 4.681]})   # CMS energy 4.681 GeV
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_2 = Selection.new
  .select_track {                # Charged track selection
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     ">=2"              # at least two positive tracks (pi+, p)
    nChrn     ">=2"              # at least two negative tracks (pi-, pbar)
  }
  .select_photon {               # Photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025      # E_gamma > 25 MeV (barrel)
    energyThreshold_e 0.050
    angle_to_track    10.0       # angle to nearest charged track > 10 deg
    nGam              ">=3"      # at least three photons (radiative + 2 from eta)
  }
  .pid(method: :probability) {   # PID: probability method, prob_cut 0.001
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # p / pbar
    identify :pion,   against: [:kaon, :proton]   # pi+ / pi-
    nprp ">=1"
    nprm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # eta -> gamma gamma (1C mass-constrained fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:gamma, :eta, :pip, :pim, :prp, :prm]) {  # nominal 4C fit to gamma eta pi+pi- p pbar
    nominal
    constrain_four_momentum
    chi2_cut 30
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :prp, :prm]) {  # competing hypothesis "2gamma + hadrons + pi+pi-"
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :prp, :prm]) {  # competing hypothesis "4gamma + hadrons + pi+pi-"
    constrain_four_momentum
  }

alg_2.note(:etap_mass_difference,
           "eta' candidates are selected with the mass-difference criterion on M(eta pi+pi-) - M(eta) " \
           "consistent with m(eta') - m(eta), to suppress the eta' -> gamma gamma background; the " \
           "chi2-based comparison against the 2gamma+hadrons and 4gamma+hadrons competing 4C hypotheses " \
           "is performed at the analysis (ROOT) level from the stored chi2 values.")
     .note(:ks0_secondary_vertex,
           "For the eta_c modes containing a K_S0 (16-mode combination), reconstruct K_S0 by a " \
           "secondary_vertex_fit of the pi+pi- pair, requiring |M(pi+pi-) - m(K_S0)| < 20 MeV and " \
           "decay length > 2 sigma; not applicable to the representative eta_c -> p pbar mode.")
     .note(:pi0_mass_constraint,
           "For the eta_c modes containing a pi0, constrain the gamma-gamma invariant mass to the " \
           "nominal pi0 mass, requiring |M(gamma gamma) - m(pi0)| < 15 MeV; not applicable to the " \
           "representative eta_c -> p pbar mode.")

alg_2.with_decay_card(decay_card_gamma_etap_etac).apply(sel_2)
root_files_2 = alg_2.execute_on([data_4681, incMC_4681, exMC_gamma_etap_etac])