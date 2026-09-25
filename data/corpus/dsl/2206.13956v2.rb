# Paper: 2206.13956v2 — Search for J/psi -> e mu (LFV)
# Type: Ordinary (counting experiment, no kinematic fit)
# Dataset: 708_3097 (J/psi, 8.998e9 events)

algorithm = Algorithm.new("JpsiEMuLFV")

algorithm.set_header([
  "EventModel/Event.h",
  "EvtRecEvent/EvtRecTrack.h",
  "EventModel/EventModel.h"
])

algorithm.set_constant(ECMS: 3.097)

decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 e+ mu-   VLL;
  Enddecay
  End
DECAYCARD

event_selection = Selection.new

event_selection.select_track do
  nTot 2
  cos_theta(-0.93, 0.93)
  Vr 0.0, 1.0
  Vz(-10.0, 10.0)
end

# Assign one track as electron, one as muon
event_selection.assign(:ep, from: :charged)
event_selection.assign(:mum, from: :charged)

event_selection.pid(method: :probability) do
  identify :ep, against: [:mum, :pip, :kp]
  prob_cut 0.8
  identify :mum, against: [:ep, :pip, :kp]
  prob_cut 0.001
end

# Photon veto: reject events with any good photon
event_selection.select_photon do
  min_energy 0.025
  min_angle 20.0
  cos_theta_barrel(0.80)
  cos_theta_endcap(0.86, 0.92)
end

# No kinematic fit needed — this is a counting experiment
# Signal region selection done in ROOT stage (|Sigma p|/sqrt(s) and E_vis/sqrt(s))

algorithm.note(:no_kinematic_fit, "This is a counting experiment with signal region selection in ROOT. No BOSS-side kinematic fit is performed. Selection: two back-to-back oppositely charged tracks (|Delta_theta| < 1.2 deg, |Delta_phi| < 1.5 deg), e/mu PID with E/p > 0.96 for e, MUC depth > 40 cm for mu, photon veto, TOF timing < 1.0 ns. Signal region: |Sigma p|/sqrt(s) <= 0.02 and 0.95 <= E_vis/sqrt(s) <= 1.04.")

algorithm.note(:background_estimation, "Background estimated from inclusive MC (J/psi decay backgrounds) and continuum data at 3.773, 3.510, 3.080 GeV. Total expected background in signal region: 36.8 +/- 4.0 events.")

algorithm.with_decay_card(decay_card).apply(event_selection)

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "708_3097"
  c.decay_card = decay_card
  c.event_count = 300_000
end

algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])