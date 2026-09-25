# 2212.05296v1: First Direct Measurement of B(Sigma+ -> Lambda e+ nu_e) at J/psi
# Double-tag (DT) method: ST anti-Sigma- -> anti-p pi0, DT Sigma+ -> Lambda e+ nu_e
# NOT a D/Lambdac TagAnalysis — hyperon-level ST/DT, ordinary analysis

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/psi -> Sigma+ anti-Sigma-
# Sigma+ -> Lambda e+ nu_e (semileptonic)
# anti-Sigma- -> anti-p- pi0
# Lambda -> p pi-
# pi0 -> gamma gamma
decay_card = <<~DECAY
Decay J/psi
  1.000 Sigma+ anti-Sigma- VSS;
Enddecay
Decay Sigma+
  1.000 Lambda e+ nu_e SLL;
Enddecay
Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
Enddecay
Decay Lambda
  1.000 p pi- PHSP;
Enddecay
Decay pi0
  1.000 gamma gamma PHSP;
Enddecay
End
DECAY

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Sigma_enu"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

algorithm = Algorithm.new("Sigma2LambdaENu", "00-00-01")
algorithm.set_header(["Sigma2LambdaENu/Sigma2LambdaENu.h"])
algorithm.set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new

# Track selection: N_Track = 4 for DT events
# Charged tracks within |cos(theta)| < 0.93
event_selection.select_track do
  nTot "==4"
  cos_theta 0.93
end

# Photon selection: standard BESIII cuts
# pi0 -> gamma gamma
event_selection.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

# pi0 reconstruction via 1C kinematic fit (chi2 < 25)
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
end

# PID for protons, electrons, pions
# Proton: L(p) > L(K) and L(p) > L(pi) and L(p) > 0.001
# e+: combined EMC+TOF likelihood L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8
# pi-: L'(pi) > L'(K) and L'(pi) > L'(e)
event_selection.pid(method: :probability) do
  identify :pim, against: [:kaon, :electron]
  prob_cut :pim, 0.001
  # Proton identification
  identify :proton, against: [:pion, :kaon]
  prob_cut :proton, 0.001
  # Electron identification — high-momentum lepton
  identify_high_momentum_leptons(:ep)
  prob_cut :ep, 0.8
end

# Lambda reconstruction: p pi- secondary vertex
# Mass window |M(p pi-) - m_Lambda| < 10 MeV/c^2
event_selection.secondary_vertex_fit([:proton, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Partial reconstruction: missing neutrino
# The signal Sigma+ -> Lambda e+ nu_e has an undetected nu_e
# Reconstruct tagged particles, infer neutrino from recoil
event_selection.partial_miss([5]) do
  # The neutrino (recID 5) is missed
  # All other particles (anti-p-, gamma, gamma, p, pi-, e+) are tagged
end

algorithm.note(:st_selection,
  "ST anti-Sigma- selection: DeltaE_tag in [-48, 52] MeV window applied post-BOSS. " \
  "ST yield extracted from M_BC^tag fit (MC shape convoluted with Gaussian + 3rd-order Chebyshev bkg)." \
  "If multiple combinations, minimum |DeltaE_tag| retained."
)

algorithm.note(:dt_selection,
  "DT selection: events with |M_BC^tag - m(Sigma-)| < 0.049 GeV/c^2. " \
  "N_Track = 4 required. e+ PID requires L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8. " \
  "pi- requires L'(pi) > L'(K) and L'(pi) > L'(e). Other track assumed proton."
)

algorithm.note(:lambda_selection,
  "Lambda candidate: |M(p pi-) - m_Lambda| < 10 MeV/c^2. " \
  "Decay length > 2x vertex resolution from IP."
)

algorithm.note(:signal_kinematics,
  "Constrained Sigma+ momentum: p_Sigma+ = -p_Sigma-_hat * sqrt(E_beam^2 - m_Sigma+^2). " \
  "M_miss^2 = E_miss^2 - p_miss^2, E_miss = E_beam - E_Lambda - E_e+. " \
  "p_miss = |p_Sigma+ - p_Lambda - p_e+|."
)

algorithm.note(:background_suppression,
  "Post-BOSS cuts: angle(anti-p, pi0) > 170 deg in Sigma- rest frame; " \
  "p_pbar in [0.16, 0.21] GeV/c in Sigma- rest frame; " \
  "M_recoil(Sigma-Lambda) + M(p pi-) - m_Lambda > -60 MeV/c^2; " \
  "M_Lambda(pi+->e+) < 1.32 GeV/c^2 (suppresses Sigma(1385)+ -> Lambda pi+ bkg)."
)

algorithm.note(:yield_extraction,
  "ST yield from binned extended ML fit to M_BC^tag (signal: MC shape * Gaussian, bkg: 3rd-order Chebyshev). " \
  "DT yield from unbinned extended ML fit to M_miss^2 (signal: MC shape * Gaussian, bkg: flat)."
)

algorithm.with_decay_card(decay_card).apply(event_selection)
algorithm.execute_on([jpsi_data, jpsi_incMC, exMC])