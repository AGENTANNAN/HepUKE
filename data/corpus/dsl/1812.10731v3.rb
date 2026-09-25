# arXiv: 1812.10731v3
# Measurement of branching fractions for Lambda_c+ -> Lambda eta pi+ and Sigma(1385)+ eta
# BESIII Collaboration
# Single-tag method at sqrt(s) = 4.600 GeV, 567 pb^-1

# ============================================================================
# Dataset preparation
# ============================================================================

ds_4600 = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000  Lambda_c+ anti-Lambda_c-   PHSP;
  Enddecay

  Decay Lambda_c+
  1.000  Lambda eta pi+             PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000  anti-p K+ pi-              PHSP;
  Enddecay

  Decay Lambda
  1.000  p+ pi-                     HypWK;
  Enddecay

  Decay eta
  1.000  gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Lambdac_LambdaEtapi"
  config.related_dataset = ds_4600
  config.events          = 200_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ============================================================================
# Algorithm: Lambda_c+ -> Lambda eta pi+
# ============================================================================

algorithm = Algorithm.new("LambdacLambdaEtapi")
algorithm.set_header(["LambdacLambdaEtapiAlg/LambdacLambdaEtapi.h"])
          .set_constant("ECMS" => [:double, 4.600])

selection = Selection.new
selection
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon]
    nprp ">=1"
    nprm ">=0"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .assign(chrgp: :pip, chrgn: :pim)
  .kinematic_fit([:pip, :Lambda, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

algorithm
  .note(:lambda_daughters_loose, "Lambda daughter tracks (p, pi-) use looser Vz < 20 cm and no Vr constraint, unlike the generic |cos(theta)| < 0.93, Vz < 10 cm, Vr < 1 cm applied to non-Lambda tracks")
  .note(:lambda_mass_window, "Lambda candidate mass window [1.111, 1.121] GeV/c^2 (3 sigma around nominal)")
  .note(:lambda_vertex, "Lambda vertex fit chi2 < 100; flight distance > 2 sigma of resolution")
  .note(:eta_mass_window, "eta -> gamma gamma mass window [505, 575] MeV/c^2 pre-Kalman-fit")
  .note(:single_tag_method, "Single-tag method: reconstruct one Lambda_c+; the other anti-Lambda_c- assumed in recoil. Delta_E and M_BC cuts/spectra handled in ROOT analysis. Delta_E signal window: [-0.03, 0.03] GeV. Best candidate: minimal |Delta_E|. M_BC fit range [2.25, 2.30] GeV/c^2 for signal yield extraction")
  .note(:sigma_star_intermediate, "Sigma(1385)+ intermediate state studied via M(Lambda pi+) spectrum; Dalitz-plot analysis and resonance fit performed in ROOT")
  .note(:total_lambdac_pairs, "Total Lambda_c+ anti-Lambda_c- pairs in data: (105.9 +/- 4.8(stat) +/- 0.5(syst)) x 10^3, from Ref.[23]")
  .note(:isr_model, "ISR simulation uses observed Lambda_c+ anti-Lambda_c- cross section line shape from Ref.[20]")
  .with_decay_card(decay_card)
  .apply(selection)

# ============================================================================
# Execute
# ============================================================================

root = algorithm.execute_on([ds_4600, incMC_4600, exMC])