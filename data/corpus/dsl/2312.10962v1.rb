# Paper 2312.10962v1: Form factors of neutral kaons and ψ(3770)→Ks0KL0 observation
# Ordinary analysis, energy scan 3.51-4.95 GeV
# e+e- → K_S0 K_L0, K_S0 → π+π-, K_L0 identified via missing momentum
# Uses ConExc generator for continuum + ψ(3770) resonance

# Data: massive energy scan — representative datasets shown
data_3773 = DatasetManager.load_real_data.find("712_3773")
data_4180 = DatasetManager.load_real_data.find("703_4180")
data_4230 = DatasetManager.load_real_data.find("703_4230")
data_4914 = DatasetManager.load_real_data.find("707_4914")
data_4946 = DatasetManager.load_real_data.find("707_4946")

inc_mc_3773 = DatasetManager.load_inclusive_mc.find("712_3773")
inc_mc_4180 = DatasetManager.load_inclusive_mc.find("703_4180")
inc_mc_4230 = DatasetManager.load_inclusive_mc.find("703_4230")
inc_mc_4914 = DatasetManager.load_inclusive_mc.find("707_4914")
inc_mc_4946 = DatasetManager.load_inclusive_mc.find("707_4946")

# Exclusive MC: e+e- → K_S0 K_L0 via continuum (VSS model)
# ConExc generator for multi-energy scan
excl_mc = DatasetManager.create_exclusive_mc_for([data_3773, data_4180, data_4230, data_4914, data_4946]) do |c|
  c.decay_card = <<~DECAY
    Decay vpho
    1.0 K_S0 K_L0 VSS;
    Decay K_S0
    1.0 pi+ pi- PHSP;
    End
  DECAY
end

# Single Algorithm + Selection applied across all energies
# Signal extraction via X = E(K_S0)/E_beam is done in ROOT
alg = Algorithm.new("Ks0KL0CrossSection", version: "00-00-01")
alg.set_header(["Ks0KL0Alg/Ks0KL0.h"])

# Note: ECMS varies by energy point — the constant here is representative for ψ(3770)
# In production, each execute_on call would use the appropriate ECMS
alg.set_constant({ "ECMS" => [:double, 3.773] })

event_selection = Selection.new("Ks0KL0Sel")

# Charged tracks: exactly 2 with net charge 0 (π+ π- from K_S0)
event_selection.select_track do |t|
  t.nTot 2
  t.nChrp 0
end

# Photon selection: variable count; used for π0 veto and K_L0 cluster check
event_selection.select_photon do |p|
  p.min_energy 0.025
  p.min_angle 20.0  # angle > 20° from charged tracks to suppress track-induced showers
end

# PID: both tracks identified as pions; electron rejection important for Bhabha suppression
event_selection.pid(method: :probability) do |pid|
  pid.identify(:pip, :pim, against: [:kp, :km, :prp, :prm, :ep, :em, :mup, :mum])
  pid.prob_cut 0.001
end

# Note: e+ e- PID rejection is critical for Bhabha suppression
event_selection.note(:electron_rejection, "L(π) > L(e) for both tracks; E/p < 1.2 GeV/(GeV/c) to suppress Bhabha events")

# K_S0 secondary vertex fit
event_selection.secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) do |v|
  v.mass_window 0.478, 0.518
  v.decay_length_gt 2.0  # > 2 cm
  v.chi2_cut 15
end

# 1C kinematic fit constraining π+π- mass to K_S0 mass
event_selection.kalman_kinematic_fit([:pip, :pim]) do |kkf|
  kkf.constrain_to_nominal_mass_of(:K_S0)
  kkf.chi2_cut 12
end

# K_L0 identification: no explicit reconstruction — missing momentum computed in ROOT
# K_L0 interaction clusters used for validation
event_selection.note(:K_L0_identification, "K_L0 identified from missing momentum. X = E(K_S0)/E_beam used as discriminating variable in ROOT")
event_selection.note(:K_L0_cluster, "K_L0 may interact with detector: search for neutral clusters in a 20° cone opposite K_S0 direction, with second moment > 20 cm2")

# Background suppression
event_selection.note(:pi0_veto, "Reject event if any γγ combination satisfies M(γγ) ∈ [0.123, 0.144] GeV/c2")
event_selection.note(:extra_energy_cut, "Total energy of neutral clusters outside K_L0 cone < 0.2 GeV to suppress K*0(892)K0+c.c. background")
event_selection.note(:EMC_timing, "EMC time difference in [0, 700] ns to suppress electronic noise")

# X = E(K_S0)/E_beam signal extraction (ROOT)
event_selection.note(:X_signal_region, "X ∈ [0.98, 1.02] for signal optimization")
event_selection.note(:X_fit, "Unbinned ML fit to X distribution: signal + e+e-→γISR ψ(3686) + K*0(892)K0+c.c. + exponential continuum background; shapes convolved with Gaussian for data-MC resolution differences")

# Cross section extraction (ROOT level)
event_selection.note(:cross_section_formula, "σ_dressed = N_obs / (ε × L × (1+δ) × B(K_S0→π+π-))")
event_selection.note(:form_factor, "|F_{K_S0 K_L0}|^2 = σ_dressed × |1-Π|^2 × 3s / (πα^2 β^3)")
event_selection.note(:psi3770_resonance, "ψ(3770)→K_S0 K_L0 extracted via coherent sum fit: σ = |BW×e^{iφ} + a/(√s)^n × √Φ|^2. Branching fraction B and relative phase φ determined from likelihood scan")

alg.with_decay_card(<<~DECAY)
  Decay vpho
  1.0 K_S0 K_L0 VSS;
  Decay K_S0
  1.0 pi+ pi- PHSP;
  End
DECAY

alg.apply(event_selection)
alg.execute_on([data_3773, data_4180, data_4230, data_4914, data_4946])

alg.note(:full_energy_scan, "Analysis covers 3.51-4.95 GeV with ~50 energy points. Only representative datasets listed here; full analysis requires all scan points from BOSS 703, 705, 706, 707, 712, 713 data")
alg.note(:iterative_isr, "ISR correction factor (1+δ) and efficiencies obtained iteratively — cross section line shape fed back into KKMC generator")
alg.note(:supplemental_material, "See Supplemental Material of paper for full list of c.m. energies, luminosities, signal yields, efficiencies, ISR factors, VP factors, and Born cross sections")