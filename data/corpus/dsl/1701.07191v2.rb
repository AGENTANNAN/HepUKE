# =============================================================================
# 1701.07191v2
# Branching fractions and angular distributions of J/psi and psi(3686)
# -> Lambda Lambda_bar and Sigma0 Sigma0_bar  (BESIII)
#
# BOSS-side spec: dataset preparation + event selection up to the final
# nominal kinematic fit. The unbinned maximum-likelihood mass fits, the
# efficiency-corrected cos(theta) fits and the alpha extraction are ROOT-level
# and are not part of this spec.
#
# The four physics channels (two parent states x two final states) have
# different parent particles, different centre-of-mass energies, different
# photon multiplicities and different mass windows, so one Algorithm object
# per channel is used (Rule T1).
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
jpsi_data  = DatasetManager.real_data.find("708_3097")   # 1310.6e6 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")
psip_data  = DatasetManager.real_data.find("709_3686")   #  447.9e6 psi(3686) events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Continuum QED background samples (e+e- -> B B_bar) used in this analysis:
# 30 pb^-1 at sqrt(s) = 3.08 GeV and 44 pb^-1 at sqrt(s) = 3.65 GeV.
cont3080_data = DatasetManager.real_data.find("708_3080")
cont3650_data = DatasetManager.real_data.find("709_3650")

### ---------------------------------------------------------------------------
### Decay cards — one per channel
### ---------------------------------------------------------------------------
# J/psi -> Lambda Lambda_bar
decay_card_jpsi_ll = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> Lambda Lambda_bar
decay_card_psip_ll = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda0 anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> Sigma0 Sigma0_bar, with Sigma0 -> gamma Lambda, Lambda -> p pi-
decay_card_jpsi_ss = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma0 anti-Sigma0 PHSP;
  Enddecay

  Decay Sigma0
  1.0000 gamma Lambda0 PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 gamma anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> Sigma0 Sigma0_bar
decay_card_psip_ss = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma0 anti-Sigma0 PHSP;
  Enddecay

  Decay Sigma0
  1.0000 gamma Lambda0 PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 gamma anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC — one signal sample per channel
### ---------------------------------------------------------------------------
exMC_jpsi_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambdalambdabar_mc"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_ll
  config.cross_section   = :default
end

exMC_psip_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_lambdalambdabar_mc"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_psip_ll
  config.cross_section   = :default
end

exMC_jpsi_ss = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_sigma0sigma0bar_mc"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_ss
  config.cross_section   = :default
end

exMC_psip_ss = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_sigma0sigma0bar_mc"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_psip_ss
  config.cross_section   = :default
end

# Shared decay-length constraint applied to every Lambda / Lambda_bar candidate:
# the paper requires a decay length > 0.2 cm and keeps the largest-value
# candidate. The decay length is not a DSL property, so it is recorded in the
# notes of each algorithm below.

### ---------------------------------------------------------------------------
### Event selection (BOSS)
### ---------------------------------------------------------------------------
### --- J/psi -> Lambda Lambda_bar ---------------------------------------------
alg1 = "JpsiToLambdaLambda"
algorithm_jpsi_ll = Algorithm.new(alg1)
algorithm_jpsi_ll
  .set_header(["#{alg1}Alg/#{alg1}.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi_ll = Selection.new
sel_jpsi_ll
  .select_track do
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # standard BESIII vertex cut along the beam direction
    Vr        1.0      # standard BESIII vertex cut in the transverse plane
    nChrp     ">=2"    # at least four charged tracks ...
    nChrn     ">=2"    # ... with total charge zero
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0     # >= 10 deg from any charged track
    energyThreshold_b 0.025    # 25 MeV in the barrel  (|cos(theta)| < 0.8)
    energyThreshold_e 0.050    # 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
  end
  .select_isolated_photon do
    angle_to_prm_track 30.0    # >= 30 deg from the anti-proton
    angle_to_prp_track 10.0    # >= 10 deg from other charged tracks
  end
  # Proton / pion separation. The paper uses a momentum criterion
  # (p > 0.5 GeV/c -> proton, p < 0.5 GeV/c -> pion); the momentum-based
  # assignment has no DSL counterpart, so the BOSS probability PID is used and
  # the momentum criterion is recorded in the notes.
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     ">=1"
    nprm     ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # take the protons out of the raw track lists
  .assign({ :chrgp => :pip, :chrgn => :pim }) # the remaining tracks are the pions from Lambda
  # Lambda -> p pi- and Lambda_bar -> p_bar pi+ secondary vertex fits
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Nominal 4C kinematic fit over the fully reconstructed final state
  .kinematic_fit([:Lambda, :Lambda_bar]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm_jpsi_ll
  .note(:kinematic_fit_chi2, "the paper applies no kinematic fit; BOSS runs the loose default 4C fit (chi2_cut 200) over the reconstructed Lambda Lambda_bar system and any tighter requirement is applied in ROOT")
  .note(:proton_pion_assignment, "the paper separates proton from pion by momentum (p > 0.5 GeV/c assumed proton, p < 0.5 GeV/c assumed pion); the probability PID is used at BOSS level and the momentum criterion is applied in ROOT")
  .note(:decay_length_selection, "Lambda(Lambda_bar) candidates are required to have a decay length larger than 0.2 cm and, when more than one candidate is found, the one with the largest decay length is retained; the decay length is not exposed as a DSL property")
  .note(:mass_window, "M(Lambda Lambda_bar) is required to be within [3.05, 3.15] GeV/c2 for J/psi -> Lambda Lambda_bar; this window is applied in ROOT")
  .note(:lambda_mass_window, "the Lambda_bar mass must satisfy |M(p_bar pi+) - M_Lambda_bar| < 3 sigma with sigma = 2.3 MeV/c2 for the J/psi data; applied in ROOT")
  .note(:angular_distribution, "the signal MC is generated with a baryon polar angular distribution 1 + alpha cos^2(theta) using the alpha values measured in this analysis (J/psi -> Lambda Lambda_bar: 0.469); the alpha values are inputs to the exclusive MC generation")
  .note(:efficiency_curve, "the detection efficiency depends on the baryon polar angle cos(theta); data/MC correction factors are determined bin-by-bin from psi -> Lambda Lambda_bar control samples and applied in ROOT")
  .note(:continuum_background, "continuum QED e+e- -> B B_bar background is estimated with the 30 pb^-1 sample at sqrt(s) = 3.08 GeV; no event survives the J/psi -> B B_bar selection")
  .note(:background_veto, "peaking backgrounds (J/psi -> Lambda Sigma0 + c.c., J/psi -> gamma K_S0 K_S0, J/psi -> gamma eta_c with eta_c -> Lambda Lambda_bar) are modelled with exclusive MC samples in the ROOT fit")
  .with_decay_card(decay_card_jpsi_ll)
  .apply(sel_jpsi_ll)

### --- psi(3686) -> Lambda Lambda_bar ------------------------------------------
alg2 = "PsipToLambdaLambda"
algorithm_psip_ll = Algorithm.new(alg2)
algorithm_psip_ll
  .set_header(["#{alg2}Alg/#{alg2}.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip_ll = Selection.new
sel_psip_ll
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  end
  .select_isolated_photon do
    angle_to_prm_track 30.0
    angle_to_prp_track 10.0
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     ">=1"
    nprm     ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({ :chrgp => :pip, :chrgn => :pim })
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :Lambda_bar]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm_psip_ll
  .note(:kinematic_fit_chi2, "the paper applies no kinematic fit; BOSS runs the loose default 4C fit (chi2_cut 200) over the reconstructed Lambda Lambda_bar system and any tighter requirement is applied in ROOT")
  .note(:proton_pion_assignment, "the paper separates proton from pion by momentum (p > 0.5 GeV/c assumed proton, p < 0.5 GeV/c assumed pion); the probability PID is used at BOSS level and the momentum criterion is applied in ROOT")
  .note(:decay_length_selection, "Lambda(Lambda_bar) candidates are required to have a decay length larger than 0.2 cm and the candidate with the largest decay length is retained; the decay length is not exposed as a DSL property")
  .note(:mass_window, "M(Lambda Lambda_bar) is required to be within [3.63, 3.75] GeV/c2 for psi(3686) -> Lambda Lambda_bar; applied in ROOT")
  .note(:lambda_mass_window, "the Lambda_bar mass must satisfy |M(p_bar pi+) - M_Lambda_bar| < 3 sigma with sigma = 4.0 MeV/c2 for the psi(3686) data; applied in ROOT")
  .note(:angular_distribution, "the signal MC is generated with 1 + alpha cos^2(theta) using the alpha value measured in this analysis (psi(3686) -> Lambda Lambda_bar: 0.824)")
  .note(:efficiency_curve, "data/MC efficiency correction factors are determined bin-by-bin in cos(theta) from psi -> Lambda Lambda_bar control samples and applied in ROOT")
  .note(:continuum_background, "continuum QED e+e- -> B B_bar background is estimated with the 44 pb^-1 sample at sqrt(s) = 3.65 GeV; only a few events survive and no peak is observed")
  .note(:background_veto, "backgrounds from psi(3686) -> pi+ pi- J/psi (J/psi -> p p_bar), psi(3686) -> Sigma0 Sigma0_bar and psi(3686) -> Lambda Sigma0 + c.c. are modelled with exclusive MC in the ROOT fit")
  .with_decay_card(decay_card_psip_ll)
  .apply(sel_psip_ll)

### --- J/psi -> Sigma0 Sigma0_bar ----------------------------------------------
alg3 = "JpsiToSigma0Sigma0bar"
algorithm_jpsi_ss = Algorithm.new(alg3)
algorithm_jpsi_ss
  .set_header(["#{alg3}Alg/#{alg3}.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi_ss = Selection.new
sel_jpsi_ss
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"    # the two photons from Sigma0 -> gamma Lambda
  end
  .select_isolated_photon do
    angle_to_prm_track 30.0
    angle_to_prp_track 10.0
    nGam               ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     ">=1"
    nprm     ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({ :chrgp => :pip, :chrgn => :pim })
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # The two photons and the two Lambda candidates enter the nominal 4C fit,
  # which selects the combination with the smallest chi2
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm_jpsi_ss
  .note(:kinematic_fit_chi2, "the paper applies no kinematic fit; BOSS runs the loose default 4C fit (chi2_cut 200) over the Lambda Lambda_bar gamma gamma system")
  .note(:proton_pion_assignment, "the paper separates proton from pion by momentum (p > 0.5 GeV/c assumed proton, p < 0.5 GeV/c assumed pion); the probability PID is used at BOSS level and the momentum criterion is applied in ROOT")
  .note(:decay_length_selection, "Lambda(Lambda_bar) candidates are required to have a decay length larger than 0.2 cm and the candidate with the largest decay length is retained; not exposed as a DSL property")
  .note(:delta_m_photon_selection, "the photon pair from Sigma0 and Sigma0_bar is chosen by minimising Delta_m = sqrt((M(Lambda gamma1) - M_Sigma0)^2 + (M(Lambda_bar gamma2) - M_Sigma0_bar)^2); this two-photon combination selection is applied in ROOT")
  .note(:mass_window, "M(Lambda Lambda_bar) is required to be within [2.82, 3.02] GeV/c2 for J/psi -> Sigma0 Sigma0_bar; applied in ROOT")
  .note(:sigma0_mass_window, "the Sigma0_bar mass must satisfy |M(p_bar pi+ gamma) - M_Sigma0_bar| < 3 sigma with sigma = 4.3 MeV/c2 for the J/psi data; applied in ROOT")
  .note(:opening_angle, "the opening angle between the reconstructed Sigma0 and Sigma0_bar in the c.m. system must exceed 178 deg for the J/psi data; applied in ROOT")
  .note(:angular_distribution, "the signal MC is generated with 1 + alpha cos^2(theta) using the alpha value measured in this analysis (J/psi -> Sigma0 Sigma0_bar: -0.449; the negative alpha confirms previous measurements)")
  .note(:efficiency_curve, "data/MC efficiency correction factors are determined bin-by-bin in cos(theta) for photon detection and Sigma0 reconstruction, and applied in ROOT")
  .note(:continuum_background, "continuum QED e+e- -> B B_bar background is estimated with the 30 pb^-1 sample at sqrt(s) = 3.08 GeV")
  .note(:background_veto, "peaking backgrounds (J/psi -> Lambda Sigma0 + c.c., J/psi -> gamma eta_c with eta_c -> Lambda Lambda_bar / Sigma0 Sigma0_bar / Lambda Sigma0 + c.c., J/psi -> Sigma0 Sigma*0 + c.c.) are modelled with exclusive MC in the ROOT fit")
  .with_decay_card(decay_card_jpsi_ss)
  .apply(sel_jpsi_ss)

### --- psi(3686) -> Sigma0 Sigma0_bar ------------------------------------------
alg4 = "PsipToSigma0Sigma0bar"
algorithm_psip_ss = Algorithm.new(alg4)
algorithm_psip_ss
  .set_header(["#{alg4}Alg/#{alg4}.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip_ss = Selection.new
sel_psip_ss
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  end
  .select_isolated_photon do
    angle_to_prm_track 30.0
    angle_to_prp_track 10.0
    nGam               ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     ">=1"
    nprm     ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({ :chrgp => :pip, :chrgn => :pim })
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm_psip_ss
  .note(:kinematic_fit_chi2, "the paper applies no kinematic fit; BOSS runs the loose default 4C fit (chi2_cut 200) over the Lambda Lambda_bar gamma gamma system")
  .note(:proton_pion_assignment, "the paper separates proton from pion by momentum (p > 0.5 GeV/c assumed proton, p < 0.5 GeV/c assumed pion); the probability PID is used at BOSS level and the momentum criterion is applied in ROOT")
  .note(:decay_length_selection, "Lambda(Lambda_bar) candidates are required to have a decay length larger than 0.2 cm and the candidate with the largest decay length is retained; not exposed as a DSL property")
  .note(:delta_m_photon_selection, "the photon pair from Sigma0 and Sigma0_bar is chosen by minimising Delta_m; this two-photon combination selection is applied in ROOT")
  .note(:mass_window, "M(Lambda Lambda_bar) is required to be within [3.34, 3.61] GeV/c2 for psi(3686) -> Sigma0 Sigma0_bar; applied in ROOT")
  .note(:sigma0_mass_window, "the Sigma0_bar mass must satisfy |M(p_bar pi+ gamma) - M_Sigma0_bar| < 3 sigma with sigma = 6.0 MeV/c2 for the psi(3686) data; applied in ROOT")
  .note(:opening_angle, "the opening angle between the reconstructed Sigma0 and Sigma0_bar in the c.m. system must exceed 178.5 deg for the psi(3686) data; applied in ROOT")
  .note(:angular_distribution, "the signal MC is generated with 1 + alpha cos^2(theta) using the alpha value measured in this analysis (psi(3686) -> Sigma0 Sigma0_bar: 0.71)")
  .note(:efficiency_curve, "data/MC efficiency correction factors are determined bin-by-bin in cos(theta) for photon detection and Sigma0 reconstruction, and applied in ROOT")
  .note(:continuum_background, "continuum QED e+e- -> B B_bar background is estimated with the 44 pb^-1 sample at sqrt(s) = 3.65 GeV")
  .note(:background_veto, "peaking backgrounds from psi(3686) -> gamma chi_cJ (chi_cJ -> Lambda Lambda_bar, J = 0,1,2) and psi(3686) -> Xi0 Xi0_bar (Xi0 -> Lambda pi0) are modelled with exclusive MC in the ROOT fit")
  .with_decay_card(decay_card_psip_ss)
  .apply(sel_psip_ss)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
root_files_jpsi_ll = algorithm_jpsi_ll.execute_on([jpsi_data, jpsi_incMC, cont3080_data, exMC_jpsi_ll])
root_files_psip_ll = algorithm_psip_ll.execute_on([psip_data, psip_incMC, cont3650_data, exMC_psip_ll])
root_files_jpsi_ss = algorithm_jpsi_ss.execute_on([jpsi_data, jpsi_incMC, cont3080_data, exMC_jpsi_ss])
root_files_psip_ss = algorithm_psip_ss.execute_on([psip_data, psip_incMC, cont3650_data, exMC_psip_ss])
