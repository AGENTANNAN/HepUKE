# DSL for arxiv:1811.08742v1
# e+e- -> K+K- cross section measurement at sqrt(s) = 2.00-3.08 GeV
# 22 energy points, ConExc mode 45 (K+K-), R-scan data

# ConExc decay card for e+e- -> K+K- (mode 45)
# Particle vpho omitted per multi-energy rule — DSL auto-injects per-point ECMS
decay_card_kk = <<~DECAYCARD
  Decay vpho
  1 ConExc 45;
  Enddecay
  Decay vhdr
  1 K+ K- PHSP;
  Enddecay
  End
DECAYCARD

# 22 R-scan energy points under BOSS 713
rscan_points = [
  "713_Rscan_2000", "713_Rscan_2050", "713_Rscan_2100", "713_Rscan_2125",
  "713_Rscan_2150", "713_Rscan_2175", "713_Rscan_2200", "713_Rscan_2232",
  "713_Rscan_2309", "713_Rscan_2386", "713_Rscan_2396", "713_Rscan_2500",
  "713_Rscan_2644", "713_Rscan_2646", "713_Rscan_2700", "713_Rscan_2800",
  "713_Rscan_2900", "713_Rscan_2950", "713_Rscan_2981", "713_Rscan_3000",
  "713_Rscan_3020", "713_Rscan_3080"
]

data_points   = rscan_points.map { |n| DatasetManager.real_data.find(n) }
incMC_points  = rscan_points.map { |n| DatasetManager.inclusive_mc.find(n) }

# Exclusive signal MC: ConExc mode 45 at all 22 scan points
sig_kk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name     = "sig_conexc_kk"
  config.events          = 100_000
  config.decay_card      = decay_card_kk
  config.cross_section   = :default  # inert for ConExc
end

# Algorithm: K+K- selection with kaon PID, kinematic fit
alg = Algorithm.new("KKCrossSection")
alg.set_header(["KKCrossSectionAlg/KKCrossSection.h"])

sel = Selection.new
  .select_track do
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .pid(has_lepton: false) do
    prob_cut 0.001
    identify :kp, "from_pion_and_proton"
    identify :km, "from_pion_and_proton"
  end
  .kinematic_fit([:kp, :km]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

# Inexpressible BOSS procedures documented as notes
alg.note(:ep_ratio_cut, "E/p < 0.7-0.8 (energy-dependent, optimised per sqrt(s))")
alg.note(:cos_theta_asymmetry, "cos_theta < 0.8 for positive track, cos_theta > -0.8 for negative track (suppress e+e- bkg)")
alg.note(:opening_angle, "Opening angle between two tracks > 179 degrees in CM frame")
alg.note(:tof_difference, "|Delta_TOF| < 3 ns between the two tracks (reject cosmic rays)")
alg.note(:momentum_window, "Momentum of negative track within (p_exp - 3*sigma_p, p_exp + 3*sigma_p); signal yield from unbinned ML fit to positive track momentum spectrum")
alg.note(:mumu_background, "Dominant background e+e- -> (gamma) mu+ mu- subtracted via MC shape in fit")
alg.note(:isr_iterative, "ISR correction factor (1+delta) determined iteratively with ConExc; converges when Born cross section changes < 0.5%")
alg.note(:jpsi_interference, "Near J/psi (sqrt(s) >= 3.0 GeV): interference between continuum and J/psi -> K+K- corrected via additional J/psi data sample")
alg.note(:form_factor, "Kaon electromagnetic form factor |F_K|^2 extracted from measured Born cross section — ROOT-level analysis")
alg.note(:resonance_fit, "Breit-Wigner fit to cross section line shape for structure at ~2.24 GeV — ROOT-level analysis")

alg.with_decay_card(decay_card_kk).apply(sel)
alg.execute_on(data_points + incMC_points + sig_kk)