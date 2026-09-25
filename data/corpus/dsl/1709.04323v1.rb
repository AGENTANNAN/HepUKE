# Paper: 1709.04323v1 — Observation of Y(2175) in e+e- → η Y(2175) at √s=3.686–4.600 GeV
# Analysis type: Ordinary multi-energy scan
# Y(2175) → φ f0(980) → φ π+π-, φ → K+K-, η → γγ
# ConExc: mode 36 = φ f0(980) (but this is not exactly η Y(2175) → η φ f0(980))
# Use KKMC + psi(4260) for continuum at higher energies;
# ψ(3686) resonant production also studied

# Load dataset tables
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# Decay card — KKMC for continuum production up to √s=4.600
# η Y(2175) with Y(2175) → φ f0(980) → φ π+π-
# ============================================================
decay_card_eta_y2175 = <<~DECAYCARD
  Decay psi(4260)
  1.0 eta Y(2175) PHSP;
  Enddecay
  Decay Y(2175)
  1.0 phi f_0 PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHOTOS;
  Enddecay
  Decay phi
  1.0 K+ K- VSS;
  Enddecay
  Decay f_0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Data energy points (from Table I of paper)
# ============================================================
data_points = [
  DatasetManager.real_data.find("709_3686"),   # 3.686 GeV — ψ(3686) data
  DatasetManager.real_data.find("712_3773"),   # 3.773 GeV
  DatasetManager.real_data.find("703_4009"),   # 4.008 GeV
  DatasetManager.real_data.find("703_4230"),   # 4.226 GeV
  DatasetManager.real_data.find("703_4260"),   # 4.258 GeV
  DatasetManager.real_data.find("703_4360"),   # 4.358 GeV
  DatasetManager.real_data.find("703_4420"),   # 4.416 GeV (round07 high-lumi)
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV
]

incMC_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") rescue nil }.compact

# ============================================================
# Exclusive MC for signal (multi-energy)
# ============================================================
sig_mc_samples = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_eta_y2175"
  config.events        = 100_000
  config.decay_card    = decay_card_eta_y2175
  config.cross_section = :default
end

# ============================================================
# Algorithm
# ============================================================
alg = Algorithm.new("EtaY2175Scan")
alg.set_header(["EtaY2175ScanAlg/EtaY2175Scan.h"])
   .set_constant({ "ECMS" => [:double, 4.260] })  # nominal

alg.with_decay_card(decay_card_eta_y2175)

# ============================================================
# Selection
# ============================================================
sel = Selection.new

# --- Charged tracks: 4 tracks (π+π-K+K-), zero net charge ---
sel.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot "==4"
  nNet "==0"
end

# --- Photon selection: at least 2 photons for η → γγ ---
sel.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

# --- PID: kaons identified via probability ---
sel.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"; nkm "==1"
  identify :pion, against: [:kaon, :proton]
  npip "==1"; npim "==1"
end

# --- Reconstruct η via Kalman kinematic fit ---
sel.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# --- 4C kinematic fit: η K+K- π+π- (η Y(2175) → η φ f0(980) → η K+K- π+π-) ---
sel.kinematic_fit([:eta, :kp, :km, :pip, :pim]) do
  constrain_four_momentum
  chi2_cut 60
  nominal
end

alg.note(:mass_windows,
  "Mass windows applied in ROOT: η [0.513, 0.578] GeV/c², " \
  "φ [1.009, 1.031] GeV/c², f0(980) [0.868, 1.089] GeV/c². " \
  "Sidebands: [μ-5W, μ-2W] and [μ+2W, μ+5W].")

alg.note(:y2175_fit,
  "Y(2175) signal extracted via unbinned maximum likelihood fit to M(φ f0(980)). " \
  "Simultaneous fit over all energy points with shared mass/width. " \
  "Born cross section: σB = Nobs / (Lint * B * ε * (1+δ) * (1+δvac)). " \
  "Line shape parameterized as σ ∝ 1/s^n with iterative ISR correction.")

alg.note(:psi3686_search,
  "Search for ψ(3686) → η Y(2175) at √s=3.686 GeV. " \
  "Continuum contribution estimated from fit to higher-energy data extrapolated down. " \
  "Upper limit on B(ψ(3686)→ηY(2175))×B(Y(2175)→φf0(980)→φπ+π-) < 2.2×10⁻⁶ at 90% CL.")

alg.note(:etap_y2175,
  "Search for e+e- → η' Y(2175) with η' → γπ+π-. " \
  "Ratio R = σ(η'Y(2175))/σ(ηY(2175)) < 0.43 at 90% CL.")

# Apply selection
alg.apply(sel)

# Execute on all data points
alg.execute_on(data_points + sig_mc_samples)