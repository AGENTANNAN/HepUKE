# Paper 2005.05072v4: Ds+ -> PP' decays — partial reconstruction at 6 energy points
# e+e- -> Ds*+ Ds- + c.c. -> gamma Ds+ Ds-
# Only Ds+ reconstructed along with soft photon from Ds*+; Ds- inferred from recoil
# 7 signal modes + 1 normalization (Ds+ -> K+K-pi+), BF relative to normalization
# Multi-energy: 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV
# Uses partial_miss (Ds- as missing particle) replacing kinematic_fit

# ===== Decay card for normalization mode Ds+ -> K+ K- pi+ =====
decay_card_norm = <<~DECAYCARD
Decay vpho
1.0 D_s*+ D_s- PHSP;
Enddecay
Decay D_s*+
1.0 gamma D_s+ PHSP;
Enddecay
Decay D_s+
1.0 K+ K- pi+ PHSP;
Enddecay
End
DECAYCARD

# 6 energy points: 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV (BOSS 703)
scan_datasets = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
]

scan_incMC = scan_datasets.map { |ds| DatasetManager.inclusive_mc.find(ds) }

exMC_norm = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name = "Ds_KpKmPip_partial"
  config.events = 100000
  config.decay_card = decay_card_norm
  config.cross_section = :default
end

# ===== Algorithm for normalization mode =====
alg_norm = Algorithm.new("Ds_to_KpKmPip_partial")
alg_norm.set_header(["Ds_to_KpKmPip_partialAlg/Ds_to_KpKmPip_partial.h"])
# Multi-energy scan: no single ECMS

# recID mapping from decay_card_norm:
#   0: vpho (CMS mother — NEVER in partial_rec/miss)
#   1: D_s*+ (reconstruct)
#   2: D_s- (MISS — not reconstructed)
#   3: gamma (from Ds*+, reconstruct)
#   4: D_s+ (reconstruct)
#   5: K+ (daughter of Ds+, auto-included)
#   6: K- (daughter of Ds+, auto-included)
#   7: pi+ (daughter of Ds+, auto-included)

sel_norm = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
  }
  # Partial reconstruction: Ds- is the missing (not reconstructed) particle
  # recID 2 = Ds-; all other particles (D_s*+, gamma, D_s+, daughters) are reconstructed
  .partial_miss([2]) {
    require_recoil_mass 1.95, 1.99
  }

alg_norm.note(:signal_modes,
  "7 signal modes (Ds+ -> PP'): K+ eta', eta' pi+, K+ eta, eta pi+, " \
  "K+ K_S0, K_S0 pi+, K+ pi0. Each mode requires its own Algorithm + Selection + decay card. " \
  "This spec covers only the normalization mode Ds+ -> K+ K- pi+. " \
  "6 additional specs needed (Rule T1: separate Algorithm per independent decay mode)."
)

alg_norm.note(:partial_reconstruction_details,
  "e+e- -> Ds*+ Ds- + c.c. Only Ds+ is reconstructed along with soft photon from Ds*+ -> gamma Ds+. " \
  "Ds- is inferred from recoil (recID 2 in partial_miss). " \
  "The combination with minimum |DeltaE| is selected, where " \
  "DeltaE = (E_Ds+ + E_gamma + E_rec) - E0. " \
  "After deltaE cut, a kinematic fit constrains M_rec(Ds+ gamma) to nominal Ds- mass. " \
  "This kinematic fit is inexpressible in DSL and must be implemented in C++ code."
)

alg_norm.note(:deltaE_requirements,
  "DeltaE cuts per mode (Table 2): " \
  "K+K-pi+: (-0.030, 0.020) GeV. " \
  "Other modes: K+eta': (-0.040,0.025); eta'pi+: (-0.040,0.025); " \
  "K+eta: (-0.045,0.025); eta pi+: (-0.045,0.025); " \
  "K+K_S0: (-0.040,0.020); K_S0 pi+: (-0.040,0.020); " \
  "K+pi0: (-0.050,0.020) GeV."
)

alg_norm.note(:band_selections,
  "After kinematic fit, candidates must satisfy: " \
  "M_rec(Ds+) in (2.100, 2.130) GeV, M(Ds+ gamma) in (2.095, 2.130) GeV (varies per mode)."
)

alg_norm.note(:intermediate_particle_selection,
  "pi0 -> gamma gamma: mass window (0.120, 0.145) GeV, 1C Kalman fit to nominal mass. " \
  "eta -> gamma gamma: mass window (0.510, 0.560) GeV, 1C Kalman fit to nominal mass. " \
  "eta' -> pi+ pi- eta: mass window (0.945, 0.970) GeV. " \
  "K_S0 -> pi+ pi-: L/sigma_L > 2, mass window (0.491, 0.505) GeV."
)

alg_norm.note(:photon_selection,
  "Soft photon from Ds*+ decay: energy > 25 MeV barrel, > 50 MeV endcap. " \
  "Angle to nearest charged track > 10 deg. " \
  "Shower time: 0 < t < 700 ns."
)

alg_norm.note(:conexc_generator,
  "Signal MC samples generated with ConExc for e+e- -> Ds*+ Ds- with Ds+ -> signal mode " \
  "and Ds- decaying inclusively. Cross section line shape from sigma(e+e- -> Ds*+ Ds-)."
)

alg_norm.with_decay_card(decay_card_norm).apply(sel_norm)
alg_norm.execute_on(scan_datasets + scan_incMC + exMC_norm)