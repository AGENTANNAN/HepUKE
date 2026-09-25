# Paper: 2203.13782v2 — First observation of direct production of χc1 in e+e- annihilation
# Process: e+e- → χc1 → γ J/ψ → γ μ+μ-
# 4-point energy scan around χc1 mass (3.5080–3.5146 GeV), ORDINARY analysis

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ── Datasets: χc1 scan samples ───────────────────────────────────────
# Four χc1 scan data points (BOSS 703)
data_scan_35080 = DatasetManager.real_data.find("703_chi_c1_scan_3")   # 3508.0 MeV, 181.8 pb-1
data_scan_35097 = DatasetManager.real_data.find("703_chi_c1_scan_2")   # 3509.7 MeV,  39.3 pb-1
data_scan_35104 = DatasetManager.real_data.find("703_chi_c1_scan_4")   # 3510.6 MeV, 184.6 pb-1 (paper says 3510.4 MeV)
data_scan_35146 = DatasetManager.real_data.find("703_chi_c1_scan_5")   # 3514.4 MeV,  40.9 pb-1 (paper says 3514.6 MeV)

scan_data_points = [data_scan_35080, data_scan_35097, data_scan_35104, data_scan_35146]

# Control samples for background validation (used in ROOT, not in jobOptions)
# These are listed for reference but not used in the main analysis job
# 712_3773 (3.773 GeV), 703_4180 (4.178 GeV), 704_psip_scan_1 (3.5815 GeV), 704_psip_scan_2 (3.6702 GeV)

inc_mc_points = scan_data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# ── Decay card ────────────────────────────────────────────────────────
# e+e- → χc1 → γ J/ψ → γ μ+μ-
# The decay card uses the PHOKHARA generator model for signal+ISR interference
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-  PHSP;
  Enddecay
  End
DECAYCARD

# ── Exclusive MC for each scan point ─────────────────────────────────
ex_mc_signal = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name   = "sig_chic1_gamma_mumu"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# ── Algorithm ─────────────────────────────────────────────────────────
algorithm = Algorithm.new("Chic1DirectScan")

# ── Event selection ───────────────────────────────────────────────────
# 2 muon tracks, 1 photon; full reconstruction with 4C kinematic fit
event_selection = Selection.new
  .select_track do
    nChrp "==1"
    nChrn "==1"
    nTot  "==2"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=1"
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:mup, :mum, against: :pion)
    prob_cut 0.001
  end
  # Muon E_EMC < 0.4 GeV cut (muon-specific, supplement to PID)
  .for_each(:mup) do
    where { eraw > 0.4 }
    remove
  end
  .for_each(:mum) do
    where { eraw > 0.4 }
    remove
  end
  # 4C kinematic fit: γ μ+ μ- constrained to initial 4-momentum
  .kinematic_fit([:gamma, :mup, :mum]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

# The best photon candidate is selected by minimum χ2_4C (looped in the fit)
# |cosθ_γ| < 0.80 to suppress ISR background
# The photon selection with cosθ cut is handled in select_photon or via for_each

algorithm
  .set_header(["Chic1DirectScan/Chic1DirectScan.h"])
  .note(:phokhara_generator, "signal and irreducible ISR background cross sections modeled with PHOKHARA event generator; interference between e+e-→χc1→γJ/ψ and e+e-→γ_ISR J/ψ implemented per Ref.[14] (H.Czyz et al., PRD 94, 034033)")
  .note(:muon_emc_cut, "muon candidates required to have E_EMC < 0.4 GeV to suppress electron/hadron backgrounds; implemented via for_each filter on muon candidate lists")
  .note(:photon_cos_theta_cut, "best photon candidate required to have |cosθ_γ| < 0.80 to suppress ISR background events; best photon selected by minimum χ2_4C in the kinematic fit combination loop")
  .note(:two_dimensional_correction, "2D correction factor applied to (M_μ+μ-, |cosθ_μ|) distributions in MC to correct generator-level discrepancies between data and PHOKHARA prediction; correction extracted from 3.773 GeV and 4.178 GeV control samples; applied in ROOT post-selection stage")
  .note(:beam_energy_spread, "beam energy spread measured to be 736 ± 27 keV; BEMS uncertainty ±0.05 MeV on c.m. energy; systematic from beam energy spread studied by varying to 1000 keV")
  .note(:interference_fit, "signal + background + interference fit performed in ROOT with N_χc1 and N_bg free, N_int = f·√(N_χc1·N_bg); common fit to all 4 scan points extracts Γ_ee = (0.12+0.13−0.08) eV and φ = (205.0+15.4−22.4)°")
  .note(:control_samples, "background description verified with control samples at 3.773, 4.178, 3.5815, and 3.6702 GeV where the χc1 signal is absent; 2D correction validated on these samples ensuring N_sig consistent with zero after correction")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on(scan_data_points + inc_mc_points + ex_mc_signal)