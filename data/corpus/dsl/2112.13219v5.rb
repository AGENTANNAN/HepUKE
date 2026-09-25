# 2112.13219v5: Measurement of e+e- → φπ+π- cross sections
# at 22 center-of-mass energies from 2.00 to 3.08 GeV.
# φ → K+K- reconstruction; one kaon allowed missing (1C kinematic fit).
# ISR correction via ConExc generator (mode 35, form B).
# Cross section: Born/dressed, iterative procedure.

# ConExc decay card for e+e- → φπ+π- (mode 35).
# Multi-energy: omit Particle vpho line; DSL auto-injects per-point √s.
decay_card_phipipi = <<~DECAYCARD
  Decay vpho
  1 gamma gamma* ConExc 35;
  Enddecay
  Decay gamma*
  1 pi+ pi- K+ K- PHSP;
  Enddecay
  End
DECAYCARD

# ── Datasets ──────────────────────────────────────────────────────
# 22 energy points from 2.00 to 3.08 GeV (R-scan, BOSS 713).
# Full energy list: 2.0000, 2.0500, 2.1000, 2.1250, 2.1500, 2.1750,
#   2.2000, 2.2324, 2.3094, 2.3864, 2.3960, 2.5000, 2.6444, 2.6464,
#   2.7000, 2.8000, 2.9000, 2.9500, 2.9810, 3.0000, 3.0200, 3.0800.
# Represented below by 3 example energies; repeat for all 22.

data_2125  = DatasetManager.real_data.find("713_Rscan_2125")
incMC_2125 = DatasetManager.inclusive_mc.find("713_Rscan_2125")

# Multi-energy scan: same ConExc card at all 22 points
data_points = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_3080"),
  # ... all 22 energies; look up exact sample names in BES3_dataset.md
]

sig_phipipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name     = "sig_phipipi"
  config.events          = 1_000_000   # 1M per energy point
  config.decay_card      = decay_card_phipipi
  config.cross_section   = :default    # required; inert for ConExc
end

# ── Algorithm ────────────────────────────────────────────────────
alg = Algorithm.new("PhiPiPiCrossSection")
alg.set_header(["PhiPiPiAlg/PhiPiPi.h"])
   .set_constant({ "ECMS" => [:double, 2.1250] })
   .with_decay_card(decay_card_phipipi)

# Event selection: π+π-K± reconstructed, one kaon missing.
# Vertex fit + 1C kinematic fit; both K+ and K- missing hypotheses tested.
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot ">=3"        # π+π- + at least 1 kaon
}
.select_photon {
  nGam ">=0"        # no photon selection in this analysis
}
.pid(method: :probability) {
  prob_cut 0.001
  # Particle type with largest probability assigned per track
  identify :pion, against: [:kaon, :proton]
  identify :kaon, against: [:pion, :proton]
  npip ">=1"
  npim ">=1"
  nkp  ">=1"        # K+ reconstructed, K- missing hypothesis
}
# Vertex fit on π+π-K+ combination, then 1C kinematic fit (missing K-).
# χ²_1C cut at 10 (paper's BOSS-level cut).
.kinematic_fit([:pip, :pim, :kp]) {
  vertex_fit([0, 1, 2])
  miss_track_of :km
  constrain_four_momentum
  chi2_cut 10
  nominal
}

alg
  .note(:missing_kaon, "One kaon allowed missing to increase efficiency (factor ~3 at 2.00 GeV, ~30% at 3.08 GeV). Both π+π-K+ (missing K-) and π+π-K- (missing K+) combinations tested; combination with smallest χ²_1C retained per event. The DSL block above shows the missing-K- hypothesis.")
  .note(:phi_window, "φ signal region: |M_K+K- − m_φ| < 0.01 GeV/c² (post-fit ROOT cut). Sidebands: [0.995, 1.005] and [1.035, 1.045] GeV/c² used for non-φ background estimation.")
  .note(:multi_energy, "Analysis at 22 c.m. energies from 2.00 to 3.08 GeV. Separate jobOptions per energy point. Dataset sample names follow 713_Rscan_[Energy_in_MeV]; look up exact names in BES3_dataset.md.")
  .note(:conexc_isr, "ISR correction factor (1+δ) from ConExc generator (mode 35). Iterative procedure: BABAR result [6] as initial input → MC → updated cross section → repeat until convergence. VP correction factor 1/|1−Π|² from F. Jegerlehner.")
  .note(:cross_section_formula, "Dressed cross section: σ = N_obs / [L · (1+δ) · ε · B(φ→K+K−)]. Born cross section: σ_B = σ_D / (1/|1−Π|²).")
  .note(:amplitude_reweighting, "Signal MC: uniform PHSP reweighted by amplitude analysis with 4 subprocesses (φf0(980), φσ, φf0(1370), φf2(1270)). Covariant tensor amplitudes; unbinned maximum likelihood fit with MINUIT. For low-statistics points, parameters from adjacent high-statistics energies used.")
  .note(:signal_extraction, "Signal yields from unbinned ML fits to M_K+K- in [2m_K±, 1.08] GeV/c². φ peak: P-wave BW ⊗ Gaussian. Background: reversed ARGUS function. No peaking background in φ signal region.")
  .note(:dominant_background, "Dominant background: e+e- → K*(892)K±π±. Non-φ contribution from φ sideband. ρK+K- accumulation also observed.")
  .note(:line_shape_fit, "Cross section line-shape fit with BW for φ(2170) + continuum term. φ(2170) mass: 2178±20±5 MeV/c², width: 140±36±16 MeV. Phase space and φf0(980)-enhanced fits both performed.")
  .apply(event_selection)

# Example execution for one energy point; repeat for all 22 energies
# alg.execute_on([data_2125, incMC_2125, sig_phipipi[1]])