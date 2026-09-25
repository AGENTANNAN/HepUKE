# Paper 1: 2312.12719v1 — Measurements of Σ EM Form Factors using untagged ISR technique
# e+e- → Σ+Σ- via ISR, at √s = 3.773 to 4.258 GeV (12 energy points)

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# Primary ISR scan energy points from BOSS 703
energy_scan_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
]

# Dedicated points for J/ψ and ψ(3686) → Σ+Σ- BF measurement
# Note: 3.773 GeV data is from BOSS 712
psi3773_data = DatasetManager.real_data.find("712_3773")

# Decay card: KKMC + psi(4260) → Σ+ Σ-, Σ+ → p π0, Σ- → pbar π0
decay_card_sigma_isr = <<~DECAYCARD
  Decay psi(4260)
  1 Sigma+ anti-Sigma- PHSP;
  Enddecay
  Decay Sigma+
  1 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1 anti-p- pi0 PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC for the ISR scan (multi-energy)
sig_sigma_scan = DatasetManager.create_exclusive_mc_for(energy_scan_points) do |config|
  config.sample_name   = "sig_sigma_isr"
  config.events        = 100_000
  config.decay_card    = decay_card_sigma_isr
  config.cross_section = :default
end

# Inclusive MC for all scan points
incMC_scan = energy_scan_points.map { |pt|
  DatasetManager.inclusive_mc.find("703_#{pt.sample_name.split('_').last}")
}

# ========== Algorithm — ISR Scan ==========
alg_sigma_isr = Algorithm.new("SigmaISR", version: '00-00-01')
alg_sigma_isr.set_header(["SigmaISRAlg/SigmaISR.h"])
# No ECMS constant — multi-energy scan, ECMS injected per point

# ========== Selection ==========
event_selection = Selection.new
  # Track selection: 2 charged tracks, net zero charge (p and pbar)
  .select_track do
    nChrp ">=1"
    nChrn ">=1"
    nTot "==2"
    nNet "==0"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  # Photon selection: >=4 photons (for two pi0 decays)
  .select_photon do
    nGam ">=4"
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
  end
  # Proton/anti-proton PID
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pip, :kp]
  end
  # pi0 reconstruction from γγ pairs (1C Kalman kinematic fit, chi2 < 25)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
  end
  # Kinematic fit: p pbar pi0 pi0 with 4C constraint
  # Use loose chi2_cut 200 (Rule T3 — tight cut in ROOT)
  .kinematic_fit([:proton_p, :proton_m, :pi0, :pi0]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

# Post-selection: best Σ+ combination by minimizing mass difference from nominal Σ+ mass
# Applied at ROOT level. Store relevant variables for offline analysis.
alg_sigma_isr.with_decay_card(decay_card_sigma_isr).apply(event_selection)
alg_sigma_isr.execute_on(energy_scan_points + [psi3773_data] + incMC_scan + sig_sigma_scan)

# Note: U_miss requirement ([-0.14, 0.06] GeV), ISR photon polar angle cut (<0.25 or >2.90 rad),
# best Σ+ combination selection, and M_gammagamma window ([-60, 40] MeV around pi0 mass)
# are applied at ROOT level. The untagged ISR technique requires ROOT-level computation
# of M(Σ+Σ-) and ISR photon kinematics.