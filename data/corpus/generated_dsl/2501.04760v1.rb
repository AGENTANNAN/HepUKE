# frozen_string_literal: true

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at √s = 3.773 GeV (20.3 fb⁻¹)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # inclusive MC at the same energy

# Decay card for the signal process D+ → e+ ν_e (EvtGen format, top mother ψ(3770))
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 e+ nu_e  PHOTOS VLL;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for D+ → e+ ν_e
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dplus_to_enue_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: "temp_for_test")

### Event selection (BOSS) — tag-based (D-tag) analysis ###
alg_name = "DplusToENu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})            # √s = 3.773 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})
   .with_decay_card(decay_card_signal)

# Tag side: single tag D^- through six hadronic modes (charge pinned to -1)
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :mBC, min: 1.863, max: 1.877                 # M_BC ∈ (1.863, 1.877) GeV/c² (explicit opt-in window)
end

# Signal side: one e+ of charge +1, plus one missing ν_e
alg.signal_side do |s|
  s.charged(ep: 1)                                       # exactly one positron
  s.require_charge 1                                     # total signal-side charge = +1
  s.missing :nu_e                                        # massless missing neutrino
end

# Kinematic fit: 4-momentum conservation + tag D^- and signal D+ mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)       # tag  D^- → known D± mass
  f.invariant_mass_of(:ep, :nu_e).constrain_to_nominal_mass_of(:Dplus)  # signal D+ → known D± mass (ν floated)
  f.chi2_cut 50
end

# BOSS-side procedures that the tag DSL cannot express — captured for the systematic-uncertainty step
alg.note(:track_selection,
  "Charged tracks required |cosθ| < 0.93, V_r < 1 cm, |V_z| < 10 cm; K_S0 daughters relaxed to |V_z| < 20 cm.")
alg.note(:kaon_pion_pid,
  "Kaon/pion identification by likelihood: L(K) > L(π) for K, L(π) > L(K) for π.")
alg.note(:k_s0_reconstruction,
  "K_S0 → π+π- candidates required a secondary-vertex fit χ² < 100, decay length L/σ_L > 2 and "
  "M(π+π-) ∈ (0.487, 0.511) GeV/c².")
alg.note(:pi0_reconstruction,
  "π0 → γγ candidates required M(γγ) ∈ (0.115, 0.150) GeV/c² and a 1C mass-constrained fit with χ² < 50.")
alg.note(:photon_selection,
  "Photons required E > 25 MeV (barrel, |cosθ| < 0.80) or E > 50 MeV (endcap, 0.86 < |cosθ| < 0.92), "
  "shower-track angle > 10°, and EMC time in [0, 700] ns.")
alg.note(:positron_pid,
  "Signal e+ identified by L(e) > 0.001, L(e)/(L(e)+L(K)+L(π)) > 0.8 and E/p > 0.8; the v1 signal-side "
  "lepton PID thresholds are fixed defaults and not DSL-tunable.")
alg.note(:bremsstrahlung_recovery,
  "Bremsstrahlung / FSR photons within 5° of the e+ direction added to the e+ four-momentum.")
alg.note(:delta_e_cut,
  "Mode-dependent ΔE cuts (from [-25, 24] to [-62, 49] MeV) applied in ROOT; when several single-tag "
  "candidates pass, the one with minimum |ΔE| is kept.")
alg.note(:background_veto,
  "Extra-particle vetoes: no extra π0, no extra good charged track, and maximum extra photon energy < 0.2 GeV.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])