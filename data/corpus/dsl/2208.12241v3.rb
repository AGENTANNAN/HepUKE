# 2208.12241v3: Observation of hindered electromagnetic Dalitz decay
# ψ(3686) → e⁺e⁻ η_c at ψ(2S) (√s = 3.686 GeV)
# Recoil mass technique, inclusive η_c decays
# Ordinary analysis (not tag-based)

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.000 e+ e- eta_c VSS;
  Enddecay
  Decay eta_c
  1.000 inclusive PHSP;
  Enddecay
  End
DECAYCARD

algorithm = Algorithm.new("Psi2StoEEeta_c", version: '00-00-01')
algorithm.set_header(["Psi2StoEEeta_c/Psi2StoEEeta_c.h"])
algorithm.set_constant({ "ECMS" => [:double, 3.686] })
algorithm.with_decay_card(decay_card)

# Dataset: ψ(2S) at 3.686 GeV
data_psi2s = DatasetManager.load_real_data.find("709_3686")
inc_mc_psi2s = DatasetManager.load_inclusive_mc.find("709_3686")

# Exclusive signal MC
signal_mc = DatasetManager.create_exclusive_mc do |c|
  c.name = "Psi2StoEEeta_c_signal"
  c.decay_card = decay_card
  c.related_dataset = data_psi2s
end

# Event selection
event_selection = Selection.new

# Track selection: |cosθ| < 0.93, Vz < 10 cm, Vr < 1 cm
event_selection.select_track do |t|
  t.cos_theta 0.93
  t.nTot { |n| n >= 2 }
  t.nNet 0
end

# Photon selection for π⁰/η veto and γ conversion veto
event_selection.select_photon do |p|
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
end

# Note: EMC time window [0, 700] ns
algorithm.note(:emc_time_window, "EMC shower time required within [0, 700] ns of event start time")

# Electron/positron PID
event_selection.pid(method: :probability) do |pid|
  pid.identify(:ep, :em, against: [:pip, :pim, :kp, :km])
  pid.prob_cut 0.001
end

# Note: Electron PID: L(e) > 0.001 and L(e)/(L(e)+L(π)+L(K)) > 0.8
# At least one e⁺e⁻ pair selected, keep all candidates
algorithm.note(:electron_pid, "L(e) > 0.001 and L(e)/(L(e)+L(π)+L(K)) > 0.8; at least one e⁺e⁻ pair; keep all e⁺e⁻ candidates")

# Note: e± momentum < 0.8 GeV/c
algorithm.note(:electron_momentum, "p(e±) < 0.8 GeV/c")

# J/ψ veto: recoil mass of π⁺π⁻ (treating each track as pion)
# outside [3.090, 3.104] GeV/c²
# Note: J/ψ veto via recoil mass
algorithm.note(:jpsi_veto, "Each track treated as pion; recoil mass of each oppositely charged pair outside [3.090, 3.104] GeV/c² to veto ψ(3686)→π⁺π⁻J/ψ")

# π⁰ veto: M_γe⁺e⁻ outside [0.115, 0.150] GeV/c²
# η veto: M_γe⁺e⁻ outside [0.505, 0.570] GeV/c²
algorithm.note(:pi0_eta_veto, "M(γe⁺e⁻) outside [0.115, 0.150] GeV/c² (π⁰ veto) and [0.505, 0.570] GeV/c² (η veto) for each γe⁺e⁻ combination")

# γ conversion veto: R_xy < 2 cm and θ_e⁺e⁻ < 40°
algorithm.note(:gamma_conversion_veto, "Veto events with e⁺e⁻ vertex R_xy > 2 cm and opening angle θ_e⁺e⁻ > 40°")

# Note: Signal extracted from recoil mass RM(e⁺e⁻) distribution
# η_c signal peak in RM(e⁺e⁻); inclusive η_c decays
algorithm.note(:recoil_mass_analysis, "Signal yield from unbinned maximum likelihood fit to RM(e⁺e⁻) distribution; η_c signal MC shape convolved with Gaussian; background: 2nd order Chebyshev polynomial + peaking background (γ*γ*→η_c, γ conversions)")

algorithm.apply(event_selection)
algorithm.execute_on([data_psi2s, inc_mc_psi2s, signal_mc])