#==============================================================================
# Paper: 1808.06946v1
# Title: Amplitude analysis of the K_S K_S system produced in radiative
#        J/psi decays
# Journal: Phys. Rev. D 98, 072003 (2018)
# Data: 1311 million J/psi events, sqrt(s) = 3.097 GeV
# Final state: gamma K_S K_S, K_S -> pi+ pi-
# Method: 6C kinematic fit (4C E-p + 2x K_S mass constraints)
#         Mass-dependent (MD) and mass-independent (MI) amplitude analyses
#         in ROOT stage (GPUPWA framework)
#==============================================================================

decay_card = <<~DECAY
Decay vpho
1.0 J/psi VSS;
Enddecay

Decay J/psi
1.0 gamma K_S0 K_S0 VSP_PWAVE;
Enddecay

Decay K_S0
1.0 pi+ pi- PHSP;
Enddecay
DECAY

# --- Algorithm ---
alg = Algorithm.new("JpsiRadiativeKsKs")

alg.set_header(["JpsiRadiativeKsKs/JpsiRadiativeKsKs.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .with_decay_card(decay_card)

# --- Selection ---
event_selection = Selection.new

# 4 charged tracks: pi+ pi- from each K_S0
event_selection.select_track do
  nTot "==4"
  nChrp "==2"
  nChrn "==2"
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
end

# At least 1 photon
event_selection.select_photon do
  nGam ">=1"
  energyThreshold_b 0.025
  energyThreshold_e 0.050
end

# No PID for pions — all tracks assumed to be pions

# Reconstruct K_S0 #1: pi+pi- -> K_S0 via secondary vertex fit
event_selection.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Reconstruct K_S0 #2 from remaining pi+ pi-
event_selection.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# 6C kinematic fit: 4C (E-p conservation) + 2x 1C (K_S0 mass constraints)
event_selection.kinematic_fit([:gamma, :K_S0, :K_S0]) do
  constrain_four_momentum
  invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
  invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
  chi2_cut 60
  nominal
end

alg.note(:ks0_flight_distance,
  "K_S0 flight distance requirements: L/sigma_L > 0 for each K_S0 and " \
  "combined sqrt((L1/sigma1)^2 + (L2/sigma2)^2) > 2.2. " \
  "Suppresses combinatorial background without K_S0.")

alg.note(:background_level,
  "Background from inclusive MC: ~0.5% of data. Continuum background " \
  "(e+e- -> gamma K_S K_S without J/psi) scaled to ~0.7%. " \
  "Both negligible for amplitude analysis; ignored in fit.")

alg.note(:amplitude_analysis,
  "Mass-dependent (MD) analysis: coherent sum of Breit-Wigner amplitudes " \
  "with covariant tensor formalism. Dominant amplitudes: f0(1710), " \
  "f0(2200), f2'(1525). Additional: f0(1370), f0(1500), f0(1790), " \
  "f0(2330), f2(1270), f2(2340), plus K1(1270) and K*(892) -> gamma K_S. " \
  "Mass-independent (MI) analysis: piecewise 0++ and 2++ amplitudes in " \
  "bins of K_S K_S invariant mass. " \
  "Fits performed with GPUPWA framework in ROOT stage. " \
  "B(J/psi -> gamma K_S K_S) = (8.1 +/- 0.4) x 10^-4.")

alg.note(:vertex_fit,
  "Each pi+pi- pair fitted to a common vertex before kinematic fit. " \
  "Charged track momenta used in kinematic fit are updated values " \
  "after vertex fit. No events have >1 combination surviving selection.")

alg.apply(event_selection)

# --- Datasets ---
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_inc_mc = DatasetManager.inclusive_mc.find("708_3097")

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_rad_ksks"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg.execute_on([jpsi_data, jpsi_inc_mc, sig_mc])