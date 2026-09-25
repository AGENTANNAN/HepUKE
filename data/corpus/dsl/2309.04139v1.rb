# Paper: 2309.04139v1 — J/psi -> Lambda Sigma^0 + c.c. hyperon form factors
# Method: Direct reconstruction with secondary vertex fits; 4C kinematic fit
# Data: 10 billion J/psi events at 3.097 GeV
# Result: R = 0.860 +/- 0.029; DeltaPhi measured; CP violation consistent with zero

### Dataset preparation ###
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/psi -> Lambda_bar Sigma^0, Sigma^0 -> gamma Lambda
# Lambda -> p pi-, Lambda_bar -> anti-p pi+
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 anti-Lambda0 Sigma0 HELAMP 0 0 1 0 1 0 0 0 1 0 1 0 0 0;
    Enddecay

    Decay Sigma0
    1.000 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_Lambda_Sigma0"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("JpsiLambdaSigma0")
alg.set_header(["JpsiLambdaSigma0Alg/JpsiLambdaSigma0.h"])
alg.set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new
event_selection.select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
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
    nGam ">=1"
  }
  # No PID applied — paper states "No particle identification is performed to maintain high efficiency"
  .assign({:chrgp => :prp, :chrgn => :pim})   # one positive = proton, others = pion hypothesis
  .remove([:prp <= :chrgp])  # remove identified proton from chrgp
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining: assume pi+, pi-
  # Lambda reconstruction: p pi- secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Lambda_bar reconstruction: anti-p pi+ secondary vertex fit
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit with gamma Lambda Lambda_bar hypothesis
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg
  .note(:vacuum_polarization_method,
    "novel method to extract hyperon form factors by exploiting vacuum polarization enhancement at J/psi peak; " \
    "J/psi -> Lambda Sigma^0 is a purely electromagnetic decay")
  .note(:sigma0_signal_region,
    "Sigma^0 candidates selected from gamma Lambda invariant mass peak; signal region defined by fit to M(gamma Lambda) spectrum")
  .note(:lambda_mass_window,
    "Lambda and Lambda_bar candidates selected with |M(p pi-) - M_Lambda| < 5 MeV/c^2")
  .note(:sigma0_mass_cut,
    "M(gamma Lambda_bar) > 1.135 GeV/c^2 and M(gamma Lambda) > 1.135 GeV/c^2 to suppress Lambda Lambda_bar background")
  .note(:chi2_4c_cut,
    "4C kinematic fit with chi2_4C < 30 applied for final selection")
  .note(:helicity_analysis,
    "6D helicity angle analysis (theta, theta_Lambda, phi_Lambda, theta_p, theta_pbar, phi_pbar) performed in ROOT; " \
    "form factors extracted from angular distribution fit")
  .note(:best_candidate_selection,
    "if multiple Lambda Lambda_bar combinations pass the vertex fit, the one with minimum sqrt((M_p_pim - M_Lambda)^2 + (M_prm_pip - M_Lambda)^2) is chosen")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])