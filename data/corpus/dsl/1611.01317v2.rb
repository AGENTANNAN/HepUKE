# Paper: 1611.01317v2 — Precise measurement of e+e- → π+π-J/ψ cross section at √s=3.77-4.60 GeV
# Analysis type: Multi-energy cross-section scan, ordinary (Algorithm + Selection)
# J/ψ reconstructed via leptonic decays (μ+μ-, e+e-)
# Uses ConExc generator (mode 90 = J/ψ π+π-) for continuum production

# Load dataset tables
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# ConExc decay card for e+e- → π+π-J/ψ (mode 90)
# ============================================================
decay_card_conexc = <<~DECAYCARD
  Decay vpho
  1 ConExc 90;
  Enddecay
  Decay vhdr
  1 pi+ pi- J/psi PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Data energy points — XYZ data (high luminosity, ~500+ pb-1 each)
# ============================================================
xyz_data_points = [
  DatasetManager.real_data.find("712_3773"),   # 3.7730 GeV
  DatasetManager.real_data.find("703_3810"),   # 3.8077 GeV (3810)
  DatasetManager.real_data.find("703_3900"),   # 3.8962 GeV
  DatasetManager.real_data.find("703_4009"),   # 4.0076 GeV
  DatasetManager.real_data.find("703_4090"),   # 4.0855 GeV
  DatasetManager.real_data.find("703_4190"),   # 4.1886/4.2077/4.2171/4.2263 GeV (high-lumi round10: 526.7)
  DatasetManager.real_data.find("703_4210"),   # 4.2077 GeV (high-lumi round10: 517.1)
  DatasetManager.real_data.find("703_4220"),   # 4.2263 GeV (high-lumi round10: 514.6)
  DatasetManager.real_data.find("703_4230"),   # 4.2263/4.2417 GeV (high-lumi: 1056.4)
  DatasetManager.real_data.find("703_4260"),   # 4.2580 GeV (828.4)
  DatasetManager.real_data.find("703_4360"),   # 4.3583 GeV (543.9)
  DatasetManager.real_data.find("703_4420"),   # 4.4156 GeV (round07 high-lumi: 1043.9)
  DatasetManager.real_data.find("703_4470"),   # 4.4671 GeV (111.09)
  DatasetManager.real_data.find("703_4600"),   # 4.5995 GeV (586.9)
]

# Also include scan data points (low luminosity, ~7-9 pb-1 each)
scan_data_points = [
  DatasetManager.real_data.find("703_4190"),   # scan round06: 43.33 (but duplicate name!)
]

# Note: multi-round datasets at same energy need special handling.
# For simplicity, use the high-luminosity entries where available.
# The full list of 130+ energy points is extensive; representative points shown.

# ============================================================
# Exclusive MC for signal (multi-energy)
# ============================================================
all_data_points = xyz_data_points
sig_mc_samples = DatasetManager.create_exclusive_mc_for(all_data_points) do |config|
  config.sample_name   = "sig_pipi_jpsi_conexc"
  config.events        = 60_000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

# ============================================================
# Algorithm — shared selection for both J/ψ→μ+μ- and e+e- modes
# The e/μ separation is done in ROOT via stored EMC energies
# ============================================================
alg = Algorithm.new("PipiJpsiScan")
alg.set_header(["PipiJpsiScanAlg/PipiJpsiScan.h"])
   .set_constant({ "ECMS" => [:double, 4.260] })  # nominal; per-point √s from MeasuredEcmsSvc

alg.with_decay_card(decay_card_conexc)

# ============================================================
# Selection — tracks, photons, PID, kinematic fit
# ============================================================
sel = Selection.new

# --- Charged tracks ---
# 4 charged tracks, zero net charge, |cosθ|<0.93, Vz<±10cm, Vr<1cm
sel.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot "==4"
  nNet "==0"
end

# --- Photon selection (not used in fit, but needed for EMC access) ---
sel.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
end

# --- PID: high-momentum leptons + pions ---
# Tracks with p > 1.06 GeV/c are leptons; EMC energy separates e/μ
# NOTE: paper applies EMC<0.35 for μ mode and EMC>1.1 for e mode;
# full separation is done in ROOT
sel.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.06,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  npip "==1"; npim "==1"
  nlp "==1"; nlm "==1"
end

alg.note(:lepton_separation,
  "Paper separates μ+μ- and e+e- modes via EMC energy: " \
  "μ candidates require EMC < 0.35 GeV, e candidates require EMC > 1.1 GeV. " \
  "This separation is applied in ROOT analysis using stored EMC energies.")

# --- 4C kinematic fit: π+π-ℓ+ℓ- ---
# χ²/ndf < 60/4
sel.kinematic_fit([:pip, :pim, :lp, :lm]) do
  constrain_four_momentum
  chi2_cut 60
  nominal
end

alg.note(:opening_angle_cuts,
  "Applied in ROOT: cos(θ_π+π-) < 0.98 to suppress radiative Bhabha/dimuon backgrounds. " \
  "For e+e- mode also: cos(θ_π±e∓) < 0.98 for photon conversion veto.")

alg.note(:jpsi_mass_window,
  "J/ψ mass window 3.08 < M(ℓ+ℓ-) < 3.12 GeV/c² applied in ROOT. " \
  "Sideband regions: 3.00-3.06 and 3.14-3.20 GeV/c² for background estimation.")

alg.note(:isr_correction,
  "ISR correction factor (1+δ) computed iteratively using KKMC program. " \
  "Detection efficiency corrected iteratively. Cross section formula: " \
  "σ = N_sig / (L_int * (1+δ) * ε * B(J/ψ→ℓ+ℓ-)).")

# Apply selection
alg.apply(sel)

# Execute on all data points
alg.execute_on(xyz_data_points + sig_mc_samples)