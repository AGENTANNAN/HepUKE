# Paper: 2112.15076v3 — Measurement of e+e- → ωπ0π0 cross section at √s = 2.0-3.08 GeV
# Type: ORDINARY + ConExc (continuum Born cross section, multi-energy scan, 19 energy points)
# arXiv: https://arxiv.org/abs/2112.15076

# === Datasets (19 R-scan energy points, 713) ===
energy_points = [
  "713_Rscan_2000", "713_Rscan_2050", "713_Rscan_2100",
  "713_Rscan_2125", "713_Rscan_2150", "713_Rscan_2175",
  "713_Rscan_2200", "713_Rscan_2232", "713_Rscan_2309",
  "713_Rscan_2386", "713_Rscan_2396", "713_Rscan_2644",
  "713_Rscan_2646", "713_Rscan_2900", "713_Rscan_2950",
  "713_Rscan_2981", "713_Rscan_3000", "713_Rscan_3020",
  "713_Rscan_3080"
]

data_samples = energy_points.map { |ep| DatasetManager.real_data.find(ep) }
incMC_samples = energy_points.map { |ep| DatasetManager.inclusive_mc.find(ep) }

# === ConExc decay card — e+e- → ωπ0π0 ===
# Form D (-2): DIY mode with user-supplied xs_user.txt
# Final state: ω(223)π0(111)π0(111), ω → π+π-π0, π0 → γγ
# → observable particles: π+π- + 6γ
decay_card_omega_pi0pi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc -2 223 111 111;
  Enddecay
  Decay vhdr
  1 pi+ pi- pi0 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# === Exclusive MC for all scan points ===
sig_mc_omega_pi0pi0 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "sig_conexc_omega_pi0pi0"
  config.events        = 100_000
  config.decay_card    = decay_card_omega_pi0pi0
  config.cross_section = :default   # inert for ConExc; required field
end

# === Event Selection (BOSS) ===
# Final state: π+π- + 6γ → reconstruct 3π0 → ω(π+π-π0)π0π0
alg = Algorithm.new("OmegaPi0Pi0")
alg.set_header(["OmegaPi0Pi0Alg/OmegaPi0Pi0.h"])
   .set_constant({ "ECMS" => [:double, 2.125] })  # placeholder; overridden per energy at execute_on
   .note(:conexc_user_xs, "ConExc DIY mode -2 requires xs_user.txt at simulation time;
     the DSL does not auto-generate xs_user.txt for ConExc. Supply xs_user.txt with
     one record per line: mass(GeV) cross_section(pb) error.")
   .note(:best_pi0_combination, "three π0 from γγ pairs: all possible γγ pairings are tried;
     the combination with smallest sum of |M(γγ)-m(π0)| is selected. Best π0 is assigned
     to ω(π+π-π0). Post-fit quality cuts applied in ROOT.")

event_selection = Selection.new

# Charged track selection: exactly 2 pions with opposite charge
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
end

# Photon selection: ≥6 good photons
event_selection.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=6"
end

# PID: identify pions against kaons and protons
event_selection.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip "==1"
  npim "==1"
end

# Reconstruct π0 from photon pairs (1C Kalman fit)
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=3"
end

# 4C kinematic fit: π+π-π0π0π0 hypothesis (ω → π+π-π0 + 2π0)
# Loose chi2 cut; optimal tight cut applied in ROOT (Rule T3)
event_selection.kinematic_fit([:pip, :pim, :pi0, :pi0, :pi0]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg.with_decay_card(decay_card_omega_pi0pi0).apply(event_selection)
alg.execute_on(data_samples + incMC_samples + sig_mc_omega_pi0pi0)