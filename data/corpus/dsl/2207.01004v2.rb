# Paper: 2207.01004v2 — Evidence for cusp effect in η' → η π0 π0
# Type: Ordinary (J/ψ → γ η', η' → η π0 π0, all-neutral)
# Dataset: 708_3097 (J/ψ, 10 billion events)
# Selection: 7+ photons, no charged tracks, 1C Kalman + 8C kinematic fit

algorithm = Algorithm.new("EtapToEtaPi0Pi0Cusp")

algorithm.set_header([
  "EventModel/Event.h",
  "EvtRecEvent/EvtRecTrack.h",
  "EventModel/EventModel.h"
])

algorithm.set_constant(ECMS: 3.097)

decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma etap   VSP_PWAVE;
  Enddecay
  Decay etap
  1.0000 eta pi0 pi0   PHSP;
  Enddecay
  End
DECAYCARD

event_selection = Selection.new

# All-neutral: no charged tracks
event_selection.select_track do
  nTot 0
end

# At least 7 photons: radiative photon + 3×γγ for η + 2π0
event_selection.select_photon do
  min_energy 0.025
  min_angle 10.0
  cos_theta_barrel(0.80)
  cos_theta_endcap(0.86, 0.92)
end

# Reconstruct eta from gamma gamma (1C Kalman fit)
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# Reconstruct two pi0 from the remaining gamma gamma pairs
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=2"
end

# 8C kinematic fit: gamma + eta + pi0 + pi0 with full constraints
# 4C (four-momentum) + 1C (eta mass) + 1C (pi0_1 mass) + 1C (pi0_2 mass) + 1C (etap mass) = 8C
# Participants: radiative gamma + eta + 2 pi0
event_selection.kinematic_fit([:gamma, :eta, :pi0, :pi0]) do
  constrain_four_momentum
  invariant_mass_of(:eta, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
  chi2_cut 100
  nominal
end

# Select best combination by smallest chi2 (handled automatically by kinematic_fit)

algorithm.with_decay_card(decay_card).apply(event_selection)

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "708_3097"
  c.decay_card = decay_card
  c.event_count = 500_000
end

algorithm.note(:all_neutral,
  "All-neutral final state: 7+ photons, no charged tracks. " \
  "Radiative photon identified as the most energetic photon. " \
  "1C Kalman kinematic fits for eta and pi0 reconstruction with chi2 < 25. " \
  "|cos(theta_pi0)| < 0.95 for photon angle in pi0 rest frame. " \
  "Photon timing: -500 < T < 500 ns relative to most energetic photon.")

algorithm.note(:best_combination,
  "If multiple gamma eta pi0 pi0 combinations exist, " \
  "the one with smallest chi2_8C is retained (handled by kinematic_fit block).")

algorithm.note(:dalitz_plot_fit,
  "Dalitz plot analysis performed in ROOT using NREFT framework. " \
  "Tree-level amplitude parameters fitted to data; cusp effect " \
  "evaluated via 1-loop and 2-loop amplitude contributions. " \
  "Scattering length a0-a2 determined from unbinned ML fits to " \
  "M^2(pi0 pi0) vs M^2(eta pi0) Dalitz plot.")

algorithm.note(:background,
  "Background from eta' -> 3pi0 (peaking) and J/psi -> omega eta " \
  "(omega -> gamma pi0, eta -> 3pi0) estimated at 0.82% contamination. " \
  "Neglected in NREFT amplitude fit.")

algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])