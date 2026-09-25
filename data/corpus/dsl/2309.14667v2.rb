# BESIII Analysis: Investigation of ΔI=1/2 rule and CP symmetry in Xi- decays
# Paper: 2309.14667v2
# Data: J/psi at 3.097 GeV, BOSS 708
# J/psi → Xi- Xibar+ → Lambda(p pi-) pi- Lambdabar(nbar pi0) pi+

# ============================================================
# Datasets
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Decay Card
# ============================================================
decay_card_Xi = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi- anti-Xi+  PHSP;
  Enddecay
  Decay Xi-
  1.0000 Lambda pi-  PHSP;
  Enddecay
  Decay anti-Xi+
  1.0000 anti-Lambda- pi+  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay anti-Lambda-
  1.0000 anti-n0 pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Jpsi_XiXibar"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_Xi
  config.cross_section   = :default
end

# ============================================================
# Analysis Algorithm
# Xi- → Lambda pi-, Lambda → p pi- (reconstructed via vertex fit)
# Xibar+ → Lambdabar pi+, Lambdabar → nbar pi0 (nbar missing, pi0 → gamma gamma)
# ============================================================
alg = Algorithm.new("XiDecayAnalysis")
alg
  .set_header(["XiDecayAnalysis/XiDecayAnalysis.h"])

sel = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=3"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 20.0
    nGam ">=2"
  end
  # Proton assignment: positive track with momentum > 0.32 GeV → proton
  # Pions: positive/negative tracks with momentum < 0.30 GeV → pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
  end
  .remove([:prp <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  # Lambda → p pi- (secondary vertex)
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Kalman fit for pi0 → gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # Full kinematic fit: Xi- (Lambda pi-) + Xibar+ (pi+, pi0) with missing anti-n0
  # Anti-n0 treated as missing particle with unknown mass
  .kinematic_fit([:Lambda, :pim, :pip, :pi0]) do
    miss_track_of :anti_n0
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg
  .note(:momentum_based_pid, "Proton: positive track p>0.32 GeV/c. Pions: positive tracks p<0.30 GeV/c, negative tracks assigned as pi-")
  .note(:xi_vertex_fit, "Xi- reconstructed via sequential decay vertex fit: Xi-→Lambda pi-→p pi- pi-; best combination by (M_pi-pi- - m_Xi)^2 + (M_p pi- - m_Lambda)^2")
  .note(:lambda_mass_window, "|M(p pi-) - m_Lambda| < 11 MeV")
  .note(:xi_mass_window, "|M(p pi- pi-) - m_Xi| < 11 MeV")
  .note(:positive_decay_length, "Both Xi- and Lambda decay lengths required positive")
  .note(:xi_cos_theta_cut, "|cos theta_Xi| < 0.84 in e+e- CM frame")
  .note(:missing_antineutron, "Antineutron treated as missing particle with unknown mass in kinematic fit; mass reconstructed from fit")
  .note(:xibar_mass_window, "|M(nbar gamma gamma pi+) - m_Xibar| < 11 MeV after kinematic fit")
  .note(:bdt_photon, "BDT classifier applied to photon candidates for signal/noise discrimination (90% signal eff, 55% bg rejection)")
  .note(:antineutron_veto, "Photon separated from Xibar direction with opening angle > 15 degrees to veto antineutron EMC showers")
  .note(:angular_analysis, "Joint angular distribution fit performed in ROOT to extract decay asymmetry parameters")
  .note(:blind_analysis, "Central values blinded using hidden-answer technique until all selections finalized")
  .with_decay_card(decay_card_Xi)
  .apply(sel)

alg.execute_on([jpsi_data, jpsi_incMC, signal_mc])