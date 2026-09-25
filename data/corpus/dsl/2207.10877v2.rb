# Paper: 2207.10877v2 — Search for ψ(3686) → Λ_c^+ Σ^- + c.c.
# Type: Ordinary (Λ_c^+ → p K^- π^+, Σ^- → p̄ π^0, 5C kinematic fit)
# Dataset: 709_3686 (ψ(3686), 448.1×10^6 events)

algorithm = Algorithm.new("PsipToLambdacSigma")

algorithm.set_header([
  "EventModel/Event.h",
  "EvtRecEvent/EvtRecTrack.h",
  "EventModel/EventModel.h"
])

algorithm.set_constant(ECMS: 3.686)

decay_card = <<~DECAYCARD
  Decay psi(3686)
  1.0000 Lambda_c+ Sigma-   PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 p+ K- pi+   PHSP;
  Enddecay
  Decay Sigma-
  1.0000 anti-p- pi0   PHSP;
  Enddecay
  End
DECAYCARD

event_selection = Selection.new

# At least 4 charged tracks: p, K^-, π^+, p̄
# Σ^- has long lifetime → looser vertex cuts for its daughter p̄
event_selection.select_track do
  nChrp 2, 10        # at least 2 positive (p, π^+)
  nChrn 2, 10        # at least 2 negative (K^-, p̄)
  nTot 4, 10
  cos_theta(-0.93, 0.93)
  Vr 0.0, 1.0
  Vz(-10.0, 10.0)
end

# PID: identify p, K^-, π^+, p̄
event_selection.pid(method: :probability) do
  identify :prp, against: [:km, :pip, :prm]   # proton
  identify :km,  against: [:prp, :pip, :prm]  # kaon
  identify :pip, against: [:prp, :km, :prm]   # pion
  identify :prm, against: [:prp, :km, :pip]   # antiproton
  prob_cut 0.0                                 # highest-probability assignment
end

# Photon selection: at least 2 photons for π^0 from Σ^- decay
event_selection.select_photon do
  min_energy 0.025
  min_angle 10.0
  cos_theta_barrel(0.80)
  cos_theta_endcap(0.86, 0.92)
end

# Reconstruct π^0 via Kalman 1C fit
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
end

# 5C kinematic fit: e+e- → p K^- π^+ p̄ π^0
# 4C (4-momentum conservation) + 1C (π^0 mass)
event_selection.kinematic_fit([:prp, :km, :pip, :prm, :pi0]) do
  constrain_four_momentum
  chi2_cut 60
  nominal
end

algorithm.note(:sigma_track_looser,
  "Charged tracks from Σ^- decay have looser vertex cuts: " \
  "|V_xy| < 10 cm, |V_z| < 20 cm (due to Σ^- long lifetime). " \
  "The select_track block applies standard cuts; " \
  "looser cuts for Σ^- daughter handled in generated C++.")

algorithm.note(:lambda_bar_veto,
  "Background from ψ(3686) → K*(892)^- p anti-Λ suppressed by requiring " \
  "M(π^+ p̄) not in [1.090, 1.130] GeV/c² (Λ̄ mass window). Applied in ROOT stage.")

algorithm.note(:kstar_veto,
  "Background from ψ(3686) → K*^0(892) p Σ^- suppressed by requiring " \
  "M(K^- π^+) not in [0.756, 1.036] GeV/c² (K*(892)^0 mass window). Applied in ROOT stage.")

algorithm.note(:sigma_mass_window,
  "Σ^- candidates required to have M(p̄ π^0) in [1.150, 1.230] GeV/c² " \
  "(3σ around Σ^- mass). Applied in ROOT stage.")

algorithm.note(:best_candidate,
  "If multiple candidates in an event, the one with smallest χ²_5C is retained.")

algorithm.note(:helix_correction,
  "Helix parameter correction applied to MC charged tracks " \
  "to improve data-MC consistency before the 5C kinematic fit.")

algorithm.note(:signal_extraction,
  "Signal yield extracted from unbinned ML fit to M(p K^- π^+) distribution. " \
  "Signal shape from MC, background from 1st-order Chebyshev polynomial. " \
  "No significant signal observed. Upper limit: B < 1.4×10^{-5} at 90% CL.")

algorithm.with_decay_card(decay_card).apply(event_selection)

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "709_3686"
  c.decay_card = decay_card
  c.event_count = 500_000
end

algorithm.execute_on([psip_data, psip_incMC, exMC_signal])